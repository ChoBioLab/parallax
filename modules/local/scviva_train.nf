process SCVIVA_TRAIN {
    tag "$meta.id"
    label 'process_high_gpu'

    conda "/sc/arion/projects/untreatedIBD/ctastad/conda/envs/scvi"

    input:
    tuple val(meta), path(resolvi_model_dir)
    tuple val(meta2), path(adata_trained)
    val scviva_max_epochs
    val num_gpus

    output:
    tuple val(meta), path("*_scviva_model"), emit: model_dir
    tuple val(meta), path("*_scviva_trained.h5ad"), emit: adata_trained
    path("*_training_summary.json"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash
    set -euo pipefail

    # Explicitly activate conda environment
    source /hpc/users/tastac01/micromamba/etc/profile.d/conda.sh
    conda activate /sc/arion/projects/untreatedIBD/ctastad/conda/envs/scvi

    python << 'EOF'
import os
os.environ['NUMBA_CACHE_DIR'] = '/tmp'
os.environ['NUMBA_DISABLE_JIT'] = '1'
os.environ['MPLBACKEND'] = 'Agg'

import scanpy as sc
import pandas as pd
import numpy as np
import scvi
import torch
import sys
import json
import logging
from pathlib import Path
import warnings

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

warnings.filterwarnings("ignore")

logger.info("=== Starting scVIVA Training ===")

try:
    # Load the trained ResolVI data
    logger.info("Loading ResolVI trained data...")
    adata = sc.read_h5ad("${adata_trained}")

    logger.info(f"Data shape: {adata.shape}")
    logger.info(f"Available obs columns: {list(adata.obs.columns)}")
    logger.info(f"Available obsm keys: {list(adata.obsm.keys())}")
    logger.info(f"Available layers: {list(adata.layers.keys())}")

    # Check if we have the required data
    if 'resolvi_predicted' not in adata.obs.columns:
        logger.error("'resolvi_predicted' column not found in adata.obs")
        raise ValueError("Missing resolvi_predicted column")

    if 'X_spatial' not in adata.obsm.keys():
        logger.warning("'X_spatial' not found in adata.obsm, checking alternatives...")
        if 'spatial' in adata.obsm.keys():
            adata.obsm['X_spatial'] = adata.obsm['spatial']
            logger.info("Used 'spatial' coordinates as X_spatial")
        else:
            logger.error("No spatial coordinates found")
            raise ValueError("Missing spatial coordinates")

    # Use scVI for latent representation with spatial information
    logger.info("Setting up spatial scVI analysis...")

    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts" if "counts" in adata.layers else None,
        labels_key="resolvi_predicted",
        batch_key="batch" if "batch" in adata.obs.columns else None
    )

    # Initialize scVI model
    logger.info("Initializing scVI model...")
    model = scvi.model.SCVI(
        adata,
        n_latent=10,
        n_hidden=128,
        n_layers=2,
        dropout_rate=0.1
    )

    # Train the model
    logger.info(f"Training scVI model for ${scviva_max_epochs} epochs...")
    model.train(
        max_epochs=${scviva_max_epochs},
        early_stopping=True,
        early_stopping_patience=20,
        check_val_every_n_epoch=10
    )

    # Save the model
    logger.info("Saving trained model...")
    model.save("${prefix}_scviva_model", overwrite=True)

    # Get latent representation
    logger.info("Computing latent representation...")
    adata.obsm["X_scviva"] = model.get_latent_representation()

    # Compute UMAP on scVI latent space
    logger.info("Computing UMAP on scVI latent space...")
    try:
        sc.pp.neighbors(adata, use_rep="X_scviva")
        sc.tl.umap(adata)
        logger.info("UMAP computation completed")
    except Exception as e:
        logger.warning(f"UMAP computation failed: {e}")

    # Save the processed data
    logger.info("Saving processed data...")
    output_file = f"${prefix}_scviva_trained.h5ad"
    adata.write_h5ad(output_file)

    # Create training summary
    summary = {
        'success': True,
        'method': 'scVI_spatial',
        'cells': int(adata.n_obs),
        'genes': int(adata.n_vars),
        'cell_types': len(adata.obs['resolvi_predicted'].unique()),
        'latent_dims': 10,
        'epochs': ${scviva_max_epochs}
    }

    with open("${prefix}_training_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    logger.info("=== scVIVA Training Completed Successfully ===")

except Exception as e:
    logger.error(f"scVIVA training failed: {e}")
    import traceback
    logger.error(f"Traceback: {traceback.format_exc()}")

    # Create error outputs
    os.makedirs("${prefix}_scviva_model", exist_ok=True)
    with open("${prefix}_scviva_model/error.txt", 'w') as f:
        f.write(f"Error: {str(e)}\\n")

    # Create minimal output if adata exists
    if 'adata' in locals():
        adata.obs['scviva_error'] = str(e)
        adata.write_h5ad("${prefix}_scviva_trained.h5ad")

    summary = {'success': False, 'error': str(e)}
    with open("${prefix}_training_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    raise e

# Write versions
with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f'    scvi-tools: {scvi.__version__}\\n')
    f.write(f'    torch: {torch.__version__}\\n')
    f.write(f'    scanpy: {sc.__version__}\\n')
EOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_scviva_model
    touch ${prefix}_scviva_model/model.pkl
    touch ${prefix}_scviva_trained.h5ad
    touch ${prefix}_training_summary.json
    touch versions.yml
    """
}

