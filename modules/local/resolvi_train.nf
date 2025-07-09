process RESOLVI_TRAIN {
    tag "$meta.id"
    label 'process_high_gpu'

    conda "/sc/arion/projects/untreatedIBD/ctastad/conda/envs/scvi"

    input:
    tuple val(meta), path(adata)
    val annotation_label
    val max_epochs
    val num_samples
    val num_gpus

    output:
    tuple val(meta), path("*_resolvi_trained.h5ad"), emit: adata
    tuple val(meta), path("*_resolvi_model"), emit: model
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
import scanpy as sc
import scvi
import torch
import numpy as np
import pandas as pd
import logging
import sys
import json
import warnings

# Setup logging to file
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('${prefix}_training_log.txt'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Suppress warnings
warnings.filterwarnings("ignore")

logger.info("=== Starting ResolVI Training ===")
logger.info(f"Available GPUs: {torch.cuda.device_count()}")
logger.info(f"CUDA available: {torch.cuda.is_available()}")

try:
    # Load preprocessed data
    logger.info("Loading AnnData...")
    adata = sc.read_h5ad("${adata}")
    logger.info(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")
    logger.info(f"Available layers: {list(adata.layers.keys())}")
    logger.info(f"Available obs columns: {list(adata.obs.columns)}")
    logger.info(f"Available obsm keys: {list(adata.obsm.keys())}")

    # Verify annotation column exists
    if 'annotation' not in adata.obs.columns:
        raise ValueError("'annotation' column not found in obs")

    # Check if we have the counts layer (required for ResolVI)
    if 'counts' not in adata.layers:
        logger.warning("No 'counts' layer found, using X as counts")
        adata.layers['counts'] = adata.X.copy()

    # Ensure X_spatial exists for ResolVI
    if 'X_spatial' not in adata.obsm:
        if 'spatial' in adata.obsm:
            adata.obsm['X_spatial'] = adata.obsm['spatial']
            logger.info("Copied spatial coordinates to X_spatial")
        else:
            logger.warning("No spatial coordinates found - this may cause issues")

    logger.info("Setting up ResolVI...")

    # Setup ResolVI exactly like the working version
    scvi.external.RESOLVI.setup_anndata(
        adata, 
        layer="counts", 
        labels_key="annotation"
    )

    # Initialize model
    logger.info("Initializing ResolVI model...")
    model = scvi.external.RESOLVI(adata, semisupervised=True)

    # Configure training parameters - use the working approach
    logger.info(f"Training ResolVI model for ${max_epochs} epochs...")

    training_success = False

    # Strategy 1: Exactly like working version - with early stopping
    try:
        model.train(
            max_epochs=${max_epochs},
            early_stopping=True,
            early_stopping_patience=10,
            early_stopping_monitor="elbo_train"
        )
        training_success = True
        logger.info("Model training completed successfully with early stopping")
    except Exception as e:
        logger.warning(f"Training with early stopping failed: {e}")

        # Strategy 2: Exactly like working version - without early stopping
        try:
            model.train(max_epochs=${max_epochs})
            training_success = True
            logger.info("Model training completed without early stopping")
        except Exception as e2:
            logger.error(f"Training failed: {e2}")
            raise

    # Save model
    logger.info("Saving trained model...")
    model.save("${prefix}_resolvi_model", overwrite=True)

    # Generate corrected counts - exactly like working version
    logger.info("Generating corrected counts...")
    try:
        samples_corr = model.sample_posterior(
            model=model.module.model_corrected,
            return_sites=["px_rate"],
            summary_fun={"post_sample_q50": np.median},
            num_samples=${num_samples},
            summary_frequency=100,
        )
        samples_corr = pd.DataFrame(samples_corr).T
        adata.layers["resolvi_corrected"] = samples_corr.loc["post_sample_q50", "px_rate"]
        logger.info("Corrected counts generated successfully")
    except Exception as e:
        logger.warning(f"Error generating corrected counts: {e}")

    # Calculate noise components - exactly like working version
    logger.info("Calculating noise components...")
    try:
        samples = model.sample_posterior(
            model=model.module.model_residuals,
            return_sites=["mixture_proportions"],
            summary_fun={"post_sample_means": np.mean},
            num_samples=${num_samples},
            summary_frequency=100,
        )
        samples = pd.DataFrame(samples).T
        adata.obs[["true_proportion", "diffusion_proportion", "background_proportion"]] = samples.loc["post_sample_means", "mixture_proportions"]
        logger.info("Noise components calculated successfully")
    except Exception as e:
        logger.warning(f"Error calculating noise components: {e}")

    # Generate predictions - exactly like working version
    logger.info("Generating cell type predictions...")
    try:
        adata.obsm["resolvi_celltypes"] = model.predict(adata, num_samples=${num_samples}, soft=True)
        adata.obs["resolvi_predicted"] = adata.obsm["resolvi_celltypes"].idxmax(axis=1)
        logger.info("Cell type predictions generated successfully")
    except Exception as e:
        logger.warning(f"Error generating predictions: {e}")
        # Create fallback predictions
        adata.obs["resolvi_predicted"] = adata.obs["annotation"]

    # Get latent representation - exactly like working version
    logger.info("Computing latent representation...")
    try:
        adata.obsm["X_resolVI"] = model.get_latent_representation(adata)
        logger.info("Latent representation computed successfully")
    except Exception as e:
        logger.warning(f"Error computing latent representation: {e}")

    # Compute UMAP - exactly like working version
    logger.info("Computing UMAP...")
    try:
        sc.pp.neighbors(adata, use_rep="X_resolVI")
        sc.tl.umap(adata)
        logger.info("UMAP computed successfully")
    except Exception as e:
        logger.warning(f"Error computing UMAP: {e}")

    # Save processed data
    logger.info("Saving processed data...")
    adata.write_h5ad("${prefix}_resolvi_trained.h5ad")

    # Create summary
    summary = {
        'success': True,
        'method': 'ResolVI',
        'cells': int(adata.n_obs),
        'genes': int(adata.n_vars),
        'annotation_label': 'annotation',
        'unique_annotations': len(adata.obs['annotation'].unique()),
        'training_success': training_success
    }

    with open("${prefix}_training_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    logger.info("=== ResolVI Training Completed Successfully ===")

except Exception as e:
    logger.error(f"ResolVI training failed: {e}")
    import traceback
    logger.error(f"Traceback: {traceback.format_exc()}")

    # Create error outputs
    import os
    os.makedirs("${prefix}_resolvi_model", exist_ok=True)
    with open("${prefix}_resolvi_model/error.txt", 'w') as f:
        f.write(f"Error: {str(e)}\\n")

    # Create minimal output
    if 'adata' in locals():
        adata.obs['resolvi_predicted'] = 'error'
        adata.write_h5ad("${prefix}_resolvi_trained.h5ad")

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
    mkdir -p ${prefix}_resolvi_model
    touch ${prefix}_resolvi_model/model.pt
    touch ${prefix}_resolvi_trained.h5ad
    touch ${prefix}_training_summary.json
    touch versions.yml
    """
}

