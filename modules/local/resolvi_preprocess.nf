process RESOLVI_PREPROCESS {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'parallax.sif' :
        'parallax:latest' }"

    input:
    tuple val(meta), path(zarr_path)
    val annotation_label

    output:
    tuple val(meta), path("*_preprocessed.h5ad"), emit: adata
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    
    # Import shared environment setup
    import sys
    sys.path.insert(0, '${projectDir}/bin')

    import spatialdata as sd
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import anndata as ad
    import logging
    import warnings
    import os
    import sys
    import os
    import tempfile
    from pathlib import Path
    import re

    # Set up container-safe numba caching
    if os.path.exists('/.singularity.d') or os.environ.get('SINGULARITY_CONTAINER'):
        # Inside Singularity - use writable temp location
        numba_cache = tempfile.mkdtemp(prefix='numba_cache_')
        os.environ['NUMBA_CACHE_DIR'] = numba_cache
        # Keep JIT enabled but use safe cache location
    else:
        # Regular environment - use default caching
        pass
    
    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    
    # Suppress warnings
    warnings.filterwarnings("ignore")
    
    def sanitize_h5_keys(adata):
        '''Replace forward slashes in AnnData keys and DataFrame columns with underscores'''
    
        # Sanitize obsm keys and DataFrame columns within obsm
        if hasattr(adata, 'obsm') and adata.obsm is not None:
            obsm_dict = dict(adata.obsm)
            for old_key in list(obsm_dict.keys()):
                new_key = re.sub(r'/', '_', old_key)
                if new_key != old_key:
                    logger.info(f"Sanitizing obsm key: '{old_key}' -> '{new_key}'")
                    obsm_dict[new_key] = obsm_dict.pop(old_key)
    
                # Check if the value is a DataFrame and sanitize column names
                if hasattr(obsm_dict[new_key], 'columns'):
                    df = obsm_dict[new_key]
                    new_columns = [re.sub(r'/', '_', str(col)) for col in df.columns]
                    if list(new_columns) != list(df.columns):
                        logger.info(f"Sanitizing columns in obsm['{new_key}']: {list(df.columns)} -> {new_columns}")
                        df.columns = new_columns
                        obsm_dict[new_key] = df
    
            adata.obsm = obsm_dict
    
        # Sanitize varm keys and DataFrame columns within varm
        if hasattr(adata, 'varm') and adata.varm is not None:
            varm_dict = dict(adata.varm)
            for old_key in list(varm_dict.keys()):
                new_key = re.sub(r'/', '_', old_key)
                if new_key != old_key:
                    logger.info(f"Sanitizing varm key: '{old_key}' -> '{new_key}'")
                    varm_dict[new_key] = varm_dict.pop(old_key)
    
                # Check if the value is a DataFrame and sanitize column names
                if hasattr(varm_dict[new_key], 'columns'):
                    df = varm_dict[new_key]
                    new_columns = [re.sub(r'/', '_', str(col)) for col in df.columns]
                    if list(new_columns) != list(df.columns):
                        logger.info(f"Sanitizing columns in varm['{new_key}']: {list(df.columns)} -> {new_columns}")
                        df.columns = new_columns
                        varm_dict[new_key] = df
    
            adata.varm = varm_dict
    
        # Sanitize obs and var column names
        if hasattr(adata, 'obs') and adata.obs is not None:
            new_obs_columns = [re.sub(r'/', '_', str(col)) for col in adata.obs.columns]
            if list(new_obs_columns) != list(adata.obs.columns):
                logger.info(f"Sanitizing obs columns: {list(adata.obs.columns)} -> {new_obs_columns}")
                adata.obs.columns = new_obs_columns
    
        if hasattr(adata, 'var') and adata.var is not None:
            new_var_columns = [re.sub(r'/', '_', str(col)) for col in adata.var.columns]
            if list(new_var_columns) != list(adata.var.columns):
                logger.info(f"Sanitizing var columns: {list(adata.var.columns)} -> {new_var_columns}")
                adata.var.columns = new_var_columns
    
        return adata

    logger.info("=== Starting ResolVI Preprocessing ===")

    # Define the annotation label to use
    annotation_label = "${annotation_label}"
    logger.info(f"Using annotation label: '{annotation_label}'")
    
    try:
        # Load spatialdata zarr store
        logger.info("Loading spatialdata zarr store...")
        sdata = sd.read_zarr("${zarr_path}")
        logger.info(f"Successfully loaded spatialdata with keys: {list(sdata.tables.keys())}")
    
        # Extract the main table (usually 'table' or the first available table)
        table_keys = list(sdata.tables.keys())
        if not table_keys:
            raise ValueError("No tables found in spatialdata object")
    
        # Try to find the main table
        main_table_key = None
        for key in ['table', 'counts', 'expression', 'adata']:
            if key in table_keys:
                main_table_key = key
                break
    
        if main_table_key is None:
            main_table_key = table_keys[0]
            logger.warning(f"Using first available table: {main_table_key}")
    
        logger.info(f"Using table: {main_table_key}")
        adata = sdata.tables[main_table_key]
    
        # Convert to AnnData if needed
        if not isinstance(adata, ad.AnnData):
            logger.info("Converting table to AnnData format...")
            adata = ad.AnnData(adata)
    
        logger.info(f"AnnData shape: {adata.shape}")
        logger.info(f"Available obs columns: {list(adata.obs.columns)}")
        logger.info(f"Available var columns: {list(adata.var.columns)}")
    
        # Check for annotation column
        annotation_found = False
        original_annotation_col = None
    
        # List of possible annotation column names to check
        possible_annotation_cols = [
            annotation_label,
            'cell_type',
            'celltype', 
            'cell_types',
            'annotation',
            'annotations',
            'cluster',
            'clusters',
            'leiden',
            'seurat_clusters'
        ]
    
        for col_name in possible_annotation_cols:
            if col_name in adata.obs.columns:
                original_annotation_col = col_name
                annotation_found = True
                logger.info(f"Found annotation column: '{col_name}'")
                break
    
        if not annotation_found:
            logger.warning(f"No annotation column found. Available columns: {list(adata.obs.columns)}")
            # Create a dummy annotation column
            adata.obs['annotation'] = 'unknown'
            logger.info(f"Created dummy annotation column: 'annotation'")
        else:
            # Standardize the annotation column name to 'annotation'
            adata.obs['annotation'] = adata.obs[original_annotation_col]
            logger.info(f"Standardized annotation column from '{original_annotation_col}' to 'annotation'")
    
        # Log annotation statistics
        if 'annotation' in adata.obs.columns:
            unique_annotations = adata.obs['annotation'].unique()
            logger.info(f"Found {len(unique_annotations)} unique cell types/annotations:")
            for i, ann in enumerate(unique_annotations[:10]):  # Show first 10
                count = sum(adata.obs['annotation'] == ann)
                logger.info(f"  {ann}: {count} cells")
            if len(unique_annotations) > 10:
                logger.info(f"  ... and {len(unique_annotations) - 10} more")
    
        # Ensure gene names are properly formatted
        logger.info("Processing gene names...")
        adata.var_names_make_unique()
    
        # Check if gene names need to be set
        if adata.var.index.name != 'gene_name' and 'gene_name' in adata.var.columns:
            logger.info("Setting gene names from 'gene_name' column")
            adata.var_names = adata.var['gene_name']
        elif adata.var.index.name != 'gene_symbol' and 'gene_symbol' in adata.var.columns:
            logger.info("Setting gene names from 'gene_symbol' column")
            adata.var_names = adata.var['gene_symbol']
    
        # Make gene names unique
        adata.var_names_make_unique()
        logger.info(f"Gene names shape: {adata.var_names.shape}")
    
        # Basic quality control metrics - FIXED for sparse matrices
        logger.info("Computing basic QC metrics...")
    
        # Calculate QC metrics if not already present
        if 'n_genes' not in adata.obs.columns:
            n_genes = (adata.X > 0).sum(axis=1)
            # Handle sparse matrix result
            if hasattr(n_genes, 'A1'):
                n_genes = n_genes.A1
            elif hasattr(n_genes, 'toarray'):
                n_genes = n_genes.toarray().flatten()
            else:
                n_genes = np.array(n_genes).flatten()
            adata.obs['n_genes'] = n_genes
    
        if 'n_counts' not in adata.obs.columns:
            n_counts = adata.X.sum(axis=1)
            # Handle sparse matrix result
            if hasattr(n_counts, 'A1'):
                n_counts = n_counts.A1
            elif hasattr(n_counts, 'toarray'):
                n_counts = n_counts.toarray().flatten()
            else:
                n_counts = np.array(n_counts).flatten()
            adata.obs['n_counts'] = n_counts
    
        # Calculate per-gene metrics
        if 'n_cells' not in adata.var.columns:
            n_cells = (adata.X > 0).sum(axis=0)
            # Handle sparse matrix result
            if hasattr(n_cells, 'A1'):
                n_cells = n_cells.A1
            elif hasattr(n_cells, 'toarray'):
                n_cells = n_cells.toarray().flatten()
            else:
                n_cells = np.array(n_cells).flatten()
            adata.var['n_cells'] = n_cells
    
        if 'total_counts' not in adata.var.columns:
            total_counts = adata.X.sum(axis=0)
            # Handle sparse matrix result
            if hasattr(total_counts, 'A1'):
                total_counts = total_counts.A1
            elif hasattr(total_counts, 'toarray'):
                total_counts = total_counts.toarray().flatten()
            else:
                total_counts = np.array(total_counts).flatten()
            adata.var['total_counts'] = total_counts
    
        # Log QC statistics
        logger.info(f"Median genes per cell: {np.median(adata.obs['n_genes'])}")
        logger.info(f"Median counts per cell: {np.median(adata.obs['n_counts'])}")
        logger.info(f"Median cells per gene: {np.median(adata.var['n_cells'])}")
    
        # Extract spatial coordinates if available
        spatial_found = False
        if 'spatial' in adata.obsm.keys():
            logger.info("Found spatial coordinates in obsm['spatial']")
            spatial_found = True
        elif 'X_spatial' in adata.obsm.keys():
            logger.info("Found spatial coordinates in obsm['X_spatial']")
            adata.obsm['spatial'] = adata.obsm['X_spatial']
            spatial_found = True
        else:
            # Check if spatial coordinates are in obs columns
            spatial_cols = []
            for col in ['x', 'y', 'X', 'Y', 'x_coord', 'y_coord', 'X_coord', 'Y_coord']:
                if col in adata.obs.columns:
                    spatial_cols.append(col)
    
            if len(spatial_cols) >= 2:
                logger.info(f"Found spatial coordinates in obs columns: {spatial_cols[:2]}")
                adata.obsm['spatial'] = adata.obs[spatial_cols[:2]].values
                spatial_found = True
    
        if not spatial_found:
            logger.warning("No spatial coordinates found. This may limit spatial analysis capabilities.")
        else:
            logger.info(f"Spatial coordinates shape: {adata.obsm['spatial'].shape}")
    
        # Ensure data is in the correct format for ResolVI
        logger.info("Preparing data for ResolVI...")
    
        # Ensure X is not sparse if it's very small, or convert to dense if needed
        if hasattr(adata.X, 'toarray'):
            if adata.X.shape[0] * adata.X.shape[1] < 1e6:  # Only convert small matrices
                logger.info("Converting sparse matrix to dense for compatibility")
                adata.X = adata.X.toarray()
    
        # Ensure proper data types
        adata.X = adata.X.astype(np.float32)
    
        # Add sample information
        adata.obs['sample_id'] = "${meta.id}"
        adata.obs['batch'] = "${meta.id}"  # For batch correction if needed
    
        # Log final statistics
        logger.info(f"Final AnnData shape: {adata.shape}")
        logger.info(f"Final obs columns: {list(adata.obs.columns)}")
        logger.info(f"Final var columns: {list(adata.var.columns)}")
        logger.info(f"Final obsm keys: {list(adata.obsm.keys())}")
    
        # Save preprocessed data
        output_file = "${prefix}_preprocessed.h5ad"
        logger.info(f"Saving preprocessed data to: {output_file}")
    
        # Sanitize keys before saving
        adata = sanitize_h5_keys(adata)
    
        adata.write_h5ad(output_file)
    
        logger.info("=== ResolVI Preprocessing Completed Successfully ===")
    
    except Exception as e:
        logger.error(f"Error during preprocessing: {str(e)}")
        logger.error(f"Error type: {type(e).__name__}")
        import traceback
        logger.error(f"Traceback: {traceback.format_exc()}")
    
        # Create a minimal error report
        error_msg = f"Preprocessing failed for sample ${meta.id}:\\n"
        error_msg += f"Error: {str(e)}\\n"
        error_msg += f"Error type: {type(e).__name__}\\n"
    
        with open("${prefix}_preprocessing_error.txt", "w") as f:
            f.write(error_msg)
    
        # Re-raise the exception to fail the process
        raise e
    
    # Write versions file
    try:
        with open("versions.yml", "w") as f:
            f.write('"${task.process}":\\n')
            f.write(f'    spatialdata: {sd.__version__}\\n')
            f.write(f'    scanpy: {sc.__version__}\\n')
            f.write(f'    pandas: {pd.__version__}\\n')
            f.write(f'    numpy: {np.__version__}\\n')
            f.write(f'    anndata: {ad.__version__}\\n')
    except Exception as e:
        logger.warning(f"Could not write versions file: {e}")
        # Create a minimal versions file
        with open("versions.yml", "w") as f:
            f.write('"${task.process}":\\n')
            f.write('    status: completed\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_preprocessed.h5ad
    touch versions.yml
    """
}

