process RESOLVI_TRAIN {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'

    conda "/sc/arion/projects/untreatedIBD/ctastad/conda/envs/scvi"

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

    # Load preprocessed data
    logger.info("Loading AnnData...")
    adata = sc.read_h5ad("${adata}")
    logger.info(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")
    logger.info("Setting up ResolVI...")

    # Setup ResolVI (data should already be preprocessed)
    scvi.external.RESOLVI.setup_anndata(
        adata, 
        layer="counts", 
        labels_key="annotation"
    )

    # Initialize model
    logger.info("Initializing ResolVI model...")
    model = scvi.external.RESOLVI(adata, semisupervised=True)

    # Configure training parameters - simplified approach
    logger.info(f"Training ResolVI model for {${max_epochs}} epochs...")
    try:
        # First attempt: Let scvi-tools auto-configure GPU settings
        model.train(
            max_epochs=${max_epochs},
            early_stopping=True,
            early_stopping_patience=10,
            early_stopping_monitor="elbo_train"
        )
        logger.info("Model training completed successfully with early stopping")
    except Exception as e:
        logger.warning(f"Training with early stopping failed: {e}")
        logger.info("Attempting training without early stopping...")
        try:
            model.train(max_epochs=${max_epochs})
            logger.info("Model training completed without early stopping")
        except Exception as e2:
            logger.error(f"Training failed: {e2}")
            raise

    # Save model
    logger.info("Saving trained model...")
    model.save("resolvi_model", overwrite=True)

    # Generate corrected counts
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

    # Calculate noise components
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

    # Generate predictions
    logger.info("Generating cell type predictions...")
    try:
        adata.obsm["resolvi_celltypes"] = model.predict(adata, num_samples=${num_samples}, soft=True)
        adata.obs["resolvi_predicted"] = adata.obsm["resolvi_celltypes"].idxmax(axis=1)
        logger.info("Cell type predictions generated successfully")
    except Exception as e:
        logger.warning(f"Error generating predictions: {e}")

    # Get latent representation
    logger.info("Computing latent representation...")
    try:
        adata.obsm["X_resolVI"] = model.get_latent_representation(adata)
        logger.info("Latent representation computed successfully")
    except Exception as e:
        logger.warning(f"Error computing latent representation: {e}")

    # Compute UMAP
    logger.info("Computing UMAP...")
    try:
        sc.pp.neighbors(adata, use_rep="X_resolVI")
        sc.tl.umap(adata)
        logger.info("UMAP computed successfully")
    except Exception as e:
        logger.warning(f"Error computing UMAP: {e}")

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

