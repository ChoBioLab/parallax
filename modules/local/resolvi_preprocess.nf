process RESOLVI_PREPROCESS {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/resolvi", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/scverse/scvi-tools:py3.12-cu12-1.3.2-dev' :
        'ghcr.io/scverse/scvi-tools:py3.12-cu12-1.3.2-dev' }"

    containerOptions '--writable-tmpfs'

    input:
    tuple val(meta), path(zarr_file)
    val annotation_label

    output:
    tuple val(meta), path("*_preprocessed.h5ad"), emit: adata
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/usr/bin/env python3
    
    # Disable numba caching to avoid container issues
    import os
    os.environ['NUMBA_CACHE_DIR'] = '/tmp'
    os.environ['NUMBA_DISABLE_JIT'] = '1'
    
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import anndata as ad
    import zarr
    import sys
    import warnings
    warnings.filterwarnings('ignore')
    
    # Set scanpy settings to avoid caching issues
    sc.settings.cache_compression = None
    
    print("Loading zarr file directly...")
    
    try:
        # Try to load as AnnData directly from zarr
        adata = ad.read_zarr("${zarr_file}")
        print(f"Loaded as AnnData from zarr: {adata.shape}")
    except Exception as e:
        print(f"Could not load as AnnData directly: {e}")
    
        # Try to load zarr store and extract data manually
        try:
            store = zarr.open("${zarr_file}", mode='r')
            print(f"Zarr store contents: {list(store.keys())}")
    
            # Look for common zarr structure patterns
            if 'tables' in store:
                tables = store['tables']
                print(f"Found tables: {list(tables.keys())}")
    
                # Get the main table (usually 'table' or first available)
                table_key = 'table' if 'table' in tables else list(tables.keys())[0]
                table_data = tables[table_key]
    
                # Load as AnnData
                adata = ad.read_zarr(f"${zarr_file}/tables/{table_key}")
    
            elif 'X' in store:
                # Direct zarr format
                adata = ad.read_zarr("${zarr_file}")
            else:
                print(f"ERROR: Unknown zarr structure. Contents: {list(store.keys())}")
                sys.exit(1)
    
        except Exception as e2:
            print(f"ERROR: Could not load zarr file: {e2}")
            sys.exit(1)
    
    print(f"Loaded data shape: {adata.shape}")
    print(f"Original obs columns: {list(adata.obs.columns)}")
    print(f"Original obsm keys: {list(adata.obsm.keys())}")
    
    # Check for annotation column
    if "${annotation_label}" not in adata.obs.columns:
        print(f"ERROR: Annotation column '${annotation_label}' not found!")
        print(f"Available columns: {list(adata.obs.columns)}")
        sys.exit(1)
    
    # Make variable names unique
    adata.var_names_make_unique()
    
    # Basic preprocessing
    print("Running basic preprocessing...")
    
    # Calculate QC metrics (handle sparse matrices properly)
    if hasattr(adata.X, 'toarray'):
        # For sparse matrices
        adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1
        adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    else:
        # For dense matrices
        adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
        adata.obs['n_counts'] = adata.X.sum(axis=1)
    
    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=10)
    sc.pp.filter_genes(adata, min_cells=3)
    
    print(f"After filtering: {adata.shape}")
    
    # Store raw data before normalization
    adata.raw = adata
    
    # More robust normalization to avoid the bug
    print("Normalizing data...")
    try:
        # Try the standard scanpy normalization
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    except Exception as e:
        print(f"Standard normalization failed: {e}")
        print("Using manual normalization...")
    
        # Manual normalization as fallback
        # Convert to dense if sparse
        if hasattr(adata.X, 'toarray'):
            X = adata.X.toarray()
        else:
            X = adata.X.copy()
    
        # Calculate size factors
        size_factors = X.sum(axis=1)
        size_factors[size_factors == 0] = 1  # Avoid division by zero
    
        # Normalize
        X_norm = X / size_factors[:, np.newaxis] * 1e4
    
        # Log transform
        X_log = np.log1p(X_norm)
    
        # Store back
        adata.X = X_log
    
        print("Manual normalization completed")
    
    # Find highly variable genes
    print("Finding highly variable genes...")
    try:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata.var.highly_variable]
    except Exception as e:
        print(f"HVG selection failed: {e}")
        print("Skipping HVG selection, using all genes")
    
    print(f"After HVG selection: {adata.shape}")
    
    # Get spatial coordinates
    if 'spatial' in adata.obsm:
        adata.obsm['X_spatial'] = adata.obsm['spatial']
        print("Using 'spatial' coordinates")
    elif 'X_spatial' in adata.obsm:
        print("Using existing 'X_spatial' coordinates")
    else:
        print("ERROR: No spatial coordinates found!")
        print(f"Available obsm keys: {list(adata.obsm.keys())}")
        sys.exit(1)
    
    print(f"Spatial coordinates shape: {adata.obsm['X_spatial'].shape}")
    
    # Ensure we have the required QC metrics
    if 'n_genes' not in adata.obs.columns:
        if hasattr(adata.X, 'toarray'):
            adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1
        else:
            adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
    
    if 'n_counts' not in adata.obs.columns:
        if hasattr(adata.X, 'toarray'):
            adata.obs['n_counts'] = adata.X.sum(axis=1).A1
        else:
            adata.obs['n_counts'] = adata.X.sum(axis=1)
    
    print(f"Final data shape: {adata.shape}")
    print(f"Final obs columns: {list(adata.obs.columns)}")
    print(f"Final obsm keys: {list(adata.obsm.keys())}")
    
    # Save preprocessed data
    output_file = f"${prefix}_preprocessed.h5ad"
    adata.write_h5ad(output_file)
    print(f"Preprocessed data saved to: {output_file}")
    
    # Create versions file
    versions = {
        'RESOLVI_PREPROCESS': {
            'python': sys.version.split()[0],
            'scanpy': sc.__version__,
            'pandas': pd.__version__,
            'numpy': np.__version__
        }
    }
    
    import yaml
    with open('versions.yml', 'w') as f:
        yaml.dump(versions, f)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_preprocessed.h5ad
    touch versions.yml
    """
}

