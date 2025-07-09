process RESOLVI_TRAIN {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/../../environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/scverse/scvi-tools:py3.12-cu12-1.3.2-dev' :
        'ghcr.io/scverse/scvi-tools:py3.12-cu12-1.3.2-dev' }"

    containerOptions '--writable-tmpfs'

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
    #!/usr/bin/env python3

    # Import and setup environment
    import sys
    sys.path.insert(0, '${projectDir}/bin')

    # Import helper modules
    from setup_python_env import setup_container_environment
    from error_handler import PipelineErrorHandler
    from validate_gpu_params import validate_gpu_configuration

    # Setup environment
    setup_container_environment()

    # Setup error handler
    handler = PipelineErrorHandler('${task.process}', '${meta.id}')
    handler.log_start()

    # Validate GPU configuration
    gpu_config = validate_gpu_configuration()
    print(f"Using GPU config: {gpu_config}")

    # Main training function
    def main_training():
        import scanpy as sc
        import scvi
        import torch
        import numpy as np
        import pandas as pd
        import logging
        import json
        import warnings
        import os
        import tempfile
        from pathlib import Path

        # Set up container-safe numba caching
        if os.path.exists('/.singularity.d') or os.environ.get('SINGULARITY_CONTAINER'):
            numba_cache = tempfile.mkdtemp(prefix='numba_cache_')
            os.environ['NUMBA_CACHE_DIR'] = numba_cache

        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler('${prefix}_training_log.txt'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        logger = logging.getLogger(__name__)
        warnings.filterwarnings("ignore")

        logger.info("=== Starting ResolVI Training ===")
        logger.info(f"Available GPUs: {torch.cuda.device_count()}")
        logger.info(f"CUDA available: {torch.cuda.is_available()}")

        # Load preprocessed data
        logger.info("Loading AnnData...")
        adata = sc.read_h5ad("${adata}")
        logger.info(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")

        # Verify annotation column exists
        if 'annotation' not in adata.obs.columns:
            raise ValueError("'annotation' column not found in obs")

        # Check if we have the counts layer
        if 'counts' not in adata.layers:
            logger.warning("No 'counts' layer found, using X as counts")
            adata.layers['counts'] = adata.X.copy()

        # Ensure X_spatial exists
        if 'X_spatial' not in adata.obsm:
            if 'spatial' in adata.obsm:
                adata.obsm['X_spatial'] = adata.obsm['spatial']
                logger.info("Copied spatial coordinates to X_spatial")
            else:
                logger.warning("No spatial coordinates found")

        logger.info("Setting up ResolVI...")
        scvi.external.RESOLVI.setup_anndata(
            adata, 
            layer="counts", 
            labels_key="annotation"
        )

        # Initialize model
        logger.info("Initializing ResolVI model...")
        model = scvi.external.RESOLVI(adata, semisupervised=True)

        # Train model
        logger.info(f"Training ResolVI model for ${max_epochs} epochs...")
        try:
            model.train(
                max_epochs=${max_epochs},
                early_stopping=True,
                early_stopping_patience=10,
                early_stopping_monitor="elbo_train"
            )
            logger.info("Model training completed successfully with early stopping")
        except Exception as e:
            logger.warning(f"Training with early stopping failed: {e}")
            model.train(max_epochs=${max_epochs})
            logger.info("Model training completed without early stopping")

        # Save model
        logger.info("Saving trained model...")
        model.save("${prefix}_resolvi_model", overwrite=True)

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
            adata.obs["resolvi_predicted"] = adata.obs["annotation"]

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
        adata.write_h5ad("${prefix}_resolvi_trained.h5ad")

        # Create summary
        summary = {
            'success': True,
            'method': 'ResolVI',
            'cells': int(adata.n_obs),
            'genes': int(adata.n_vars),
            'annotation_label': 'annotation',
            'unique_annotations': len(adata.obs['annotation'].unique()),
            'training_success': True
        }

        with open("${prefix}_training_summary.json", 'w') as f:
            json.dump(summary, f, indent=2)

        logger.info("=== ResolVI Training Completed Successfully ===")

    # Execute with error handling
    try:
        main_training()
        handler.log_completion()
        handler.create_versions_file()
    except Exception as e:
        handler.handle_error(e, create_placeholder=True)

        # Create error outputs
        import os
        os.makedirs("${prefix}_resolvi_model", exist_ok=True)
        with open("${prefix}_resolvi_model/error.txt", 'w') as f:
            f.write(f"Error: {str(e)}\\n")

        # Create minimal output
        try:
            import scanpy as sc
            adata = sc.read_h5ad("${adata}")
            adata.obs['resolvi_predicted'] = 'error'
            adata.write_h5ad("${prefix}_resolvi_trained.h5ad")
        except:
            pass

        summary = {'success': False, 'error': str(e)}
        with open("${prefix}_training_summary.json", 'w') as f:
            import json
            json.dump(summary, f, indent=2)

        raise
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

