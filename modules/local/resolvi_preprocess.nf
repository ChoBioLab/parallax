process RESOLVI_PREPROCESS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::scanpy bioconda::spatialdata"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.9.3--pyhd8ed1ab_0' :
        'quay.io/biocontainers/scanpy:1.9.3--pyhd8ed1ab_0' }"

    input:
    tuple val(meta), path(zarr_store)

    output:
    tuple val(meta), path("*.h5ad"), emit: adata
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3

    import spatialdata as sd
    import scanpy as sc
    import pandas as pd
    import numpy as np

    # Load spatialdata zarr store
    sdata = sd.read_zarr("${zarr_store}")

    # Convert to AnnData (adjust based on your sopa output structure)
    adata = sdata.tables['table']  # Adjust key as needed

    # Setup spatial coordinates
    if 'X_spatial' not in adata.obsm:
        # Extract coordinates from spatialdata
        coords = sdata.points['transcripts'][['x', 'y']].values  # Adjust as needed
        adata.obsm['X_spatial'] = coords

    # Ensure required columns exist
    if 'annotation' not in adata.obs and 'cell_type' in adata.obs:
        adata.obs['annotation'] = adata.obs['cell_type']

    # Setup counts layer
    if 'counts' not in adata.layers:
        if 'raw_counts' in adata.layers:
            adata.layers['counts'] = adata.layers['raw_counts']
        else:
            adata.layers['counts'] = adata.X.copy()

    # Save preprocessed data
    adata.write_h5ad("${meta.id}_preprocessed.h5ad")

    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    scanpy: {sc.__version__}\\n')
        f.write(f'    spatialdata: {sd.__version__}\\n')
    """
}
