process RESOLVI_PREPROCESS {
    tag "$meta.id"
    label 'process_medium'

    conda "/sc/arion/projects/untreatedIBD/ctastad/conda/envs/scvi"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("*_preprocessed.h5ad"), emit: adata
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python3 << 'EOF'
import scanpy as sc
import pandas as pd
import numpy as np
import scvi
import decoupler as dc
import spatialdata as sd
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

try:
    # Check if input is zarr store or h5ad
    input_path = "${adata}"
    logger.info(f"Loading data from: {input_path}")

    if input_path.endswith('.zarr'):
        # Load spatialdata zarr store
        sdata = sd.read_zarr(input_path)
        logger.info(f"Loaded spatialdata with tables: {list(sdata.tables.keys())}")

        # Extract AnnData from first table
        table_key = list(sdata.tables.keys())[0]
        adata = sdata.tables[table_key]
        logger.info(f"Extracted AnnData from table '{table_key}' with shape: {adata.shape}")

    else:
        # Load as H5AD
        adata = sc.read_h5ad(input_path)
        logger.info(f"Loaded H5AD with shape: {adata.shape}")

    # Setup spatial coordinates for ResolVI
    logger.info("Setting up spatial coordinates...")
    if "X_spatial" not in adata.obsm:
        logger.info("'X_spatial' not found in obsm. Attempting to find coordinates...")
        if "spatial" in adata.obsm:
            logger.info("Found coordinates in 'spatial', copying to 'X_spatial'")
            adata.obsm["X_spatial"] = adata.obsm["spatial"].copy()
        else:
            logger.info("Looking for coordinates in observation data columns...")
            coord_pairs = [
                ["x_centroid", "y_centroid"],
                ["x_position", "y_position"],
                ["x", "y"],
                ["X", "Y"],
            ]
            coords_found = False
            for col_x, col_y in coord_pairs:
                if col_x in adata.obs.columns and col_y in adata.obs.columns:
                    logger.info(f"Found coordinates in obs columns: {col_x}, {col_y}")
                    adata.obsm["X_spatial"] = adata.obs[[col_x, col_y]].values
                    coords_found = True
                    break
            if not coords_found:
                logger.error("No spatial coordinates found")
                raise ValueError("No spatial coordinates found in data")
    else:
        logger.info("'X_spatial' already present in obsm")

    logger.info(f"Spatial coordinates shape: {adata.obsm['X_spatial'].shape}")

    # Setup annotation column for ResolVI
    logger.info("Setting up annotation column...")
    if "annotation" not in adata.obs:
        if "cell_type" in adata.obs:
            logger.info("Creating 'annotation' column from 'cell_type'")
            adata.obs["annotation"] = adata.obs["cell_type"]
        else:
            logger.warning("Neither 'annotation' nor 'cell_type' column found")
    else:
        logger.info("'annotation' column already present")

    # Setup counts layer
    logger.info("Setting up counts layer...")
    if "counts" not in adata.layers:
        if "raw_counts" in adata.layers:
            logger.info("Using 'raw_counts' as 'counts'")
            adata.layers["counts"] = adata.layers["raw_counts"]
        elif "count" in adata.layers:
            logger.info("Using 'count' as 'counts'")
            adata.layers["counts"] = adata.layers["count"]
        elif adata.X is not None:
            logger.info("Using X as 'counts'")
            adata.layers["counts"] = adata.X.copy()
        else:
            logger.error("No suitable layer for 'counts' and X is None")
            raise ValueError("No suitable counts data found")
    else:
        logger.info("'counts' layer already present")

    # Basic preprocessing
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    # Save preprocessed data
    adata.write("${prefix}_preprocessed.h5ad")
    logger.info(f"Saved preprocessed data to ${prefix}_preprocessed.h5ad")

except Exception as e:
    logger.error(f"Error processing data: {e}")
    raise

# Write versions
versions = {
    "scanpy": sc.__version__,
    "scvi-tools": scvi.__version__,
    "decoupler": dc.__version__,
    "spatialdata": sd.__version__
}

with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    for tool, version in versions.items():
        f.write(f'    {tool}: {version}\\n')
EOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_preprocessed.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: 1.9.6
        scvi-tools: 1.0.4
        decoupler: 1.4.0
        spatialdata: 0.0.15
    END_VERSIONS
    """
}

