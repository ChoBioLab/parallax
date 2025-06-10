process RESOLVI_TRAIN {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'

    conda "bioconda::scvi-tools"
    container 'oras://community.wave.seqera.io/library/pip_decoupler_scanpy_scvi-tools:8124a7e473830fad'

    input:
    tuple val(meta), path(adata)
    val marker_genes
    val max_epochs
    val num_samples
    val num_gpus

    output:
    tuple val(meta), path("resolvi_model/"), emit: model
    tuple val(meta), path("*_trained.h5ad"), emit: adata_processed
    tuple val(meta), path("training_log.txt"), emit: logs
    path "versions.yml", emit: versions

    script:
    def gpu_config = ""
    if (num_gpus != null) {
        if (num_gpus == -1) {
            gpu_config = "accelerator = 'gpu' if torch.cuda.is_available() else 'cpu'; devices = 'auto'"
        } else if (num_gpus == 0) {
            gpu_config = "accelerator = 'cpu'; devices = 'auto'"
        } else {
            gpu_config = "accelerator = 'gpu' if torch.cuda.is_available() else 'cpu'; devices = min(${num_gpus}, torch.cuda.device_count()) if torch.cuda.is_available() else 'auto'"
        }
    } else {
        // Default to 4 GPUs as in LSF script
        gpu_config = "accelerator = 'gpu' if torch.cuda.is_available() else 'cpu'; devices = min(4, torch.cuda.device_count()) if torch.cuda.is_available() else 'auto'"
    }

    """
    #!/usr/bin/env python3

    import scanpy as sc
    import scvi
    import torch
    import numpy as np
    import pandas as pd
    import logging
    import sys

    # Setup logging to file
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('training_log.txt'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger = logging.getLogger(__name__)

    logger.info("=== Starting ResolVI Training ===")
    logger.info(f"Available GPUs: {torch.cuda.device_count()}")
    logger.info(f"CUDA available: {torch.cuda.is_available()}")

    # Load data
    logger.info("Loading AnnData...")
    adata = sc.read_h5ad("${adata}")
    logger.info(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")

    # Setup ResolVI
    logger.info("Setting up ResolVI...")
    scvi.external.RESOLVI.setup_anndata(
        adata, 
        layer="counts", 
        labels_key="annotation"
    )

    # Initialize model
    logger.info("Initializing ResolVI model...")
    model = scvi.external.RESOLVI(adata, semisupervised=True)

    # Configure GPU usage based on LSF script parameters
    ${gpu_config}

    logger.info(f"Training configuration - Accelerator: {accelerator}, Devices: {devices}")
    logger.info(f"Max epochs: ${max_epochs}")

    # Train model with configuration matching LSF script
    logger.info("Starting model training...")
    try:
        model.train(
            max_epochs=${max_epochs},
            accelerator=accelerator,
            devices=devices,
            early_stopping=True,
            early_stopping_patience=10,
            early_stopping_monitor="elbo_train"
        )
        logger.info("Model training completed successfully with early stopping.")
    except Exception as e:
        logger.warning(f"Training with early stopping failed: {e}")
        logger.info("Retrying without early stopping...")
        model.train(
            max_epochs=${max_epochs},
            accelerator=accelerator,
            devices=devices
        )
        logger.info("Model training completed without early stopping.")

    # Save model
    logger.info("Saving trained model...")
    model.save("resolvi_model", overwrite=True)

    # Generate corrected counts
    logger.info("Generating corrected counts...")
    samples_corr = model.sample_posterior(
        model=model.module.model_corrected,
        return_sites=["px_rate"],
        summary_fun={"post_sample_q50": np.median},
        num_samples=${num_samples},
        summary_frequency=100,
    )
    samples_corr = pd.DataFrame(samples_corr).T
    adata.layers["resolvi_corrected"] = samples_corr.loc["post_sample_q50", "px_rate"]

    # Calculate noise components
    logger.info("Calculating noise components...")
    samples = model.sample_posterior(
        model=model.module.model_residuals,
        return_sites=["mixture_proportions"],
        summary_fun={"post_sample_means": np.mean},
        num_samples=${num_samples},
        summary_frequency=100,
    )
    samples = pd.DataFrame(samples).T
    adata.obs[["true_proportion", "diffusion_proportion", "background_proportion"]] = samples.loc["post_sample_means", "mixture_proportions"]

    # Generate predictions
    logger.info("Generating cell type predictions...")
    adata.obsm["resolvi_celltypes"] = model.predict(adata, num_samples=${num_samples}, soft=True)
    adata.obs["resolvi_predicted"] = adata.obsm["resolvi_celltypes"].idxmax(axis=1)

    # Get latent representation
    logger.info("Computing latent representation...")
    adata.obsm["X_resolVI"] = model.get_latent_representation(adata)

    # Compute UMAP
    logger.info("Computing UMAP...")
    sc.pp.neighbors(adata, use_rep="X_resolVI")
    sc.tl.umap(adata)

    # Save processed data
    logger.info("Saving processed data...")
    adata.write_h5ad("${meta.id}_trained.h5ad")

    logger.info("=== ResolVI Training Completed Successfully ===")

    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    scvi-tools: {scvi.__version__}\\n')
        f.write(f'    torch: {torch.__version__}\\n')
        f.write(f'    scanpy: {sc.__version__}\\n')
    """
}
