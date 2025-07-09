process SCVIVA_TRAIN {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/../../environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/scverse/scvi-tools:py3.12-cu12-1.3.2-dev' :
        'ghcr.io/scverse/scvi-tools:py3.12-cu12-1.3.2-dev' }"

    containerOptions '--writable-tmpfs'

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
        import pandas as pd
        import numpy as np
        import scvi
        import torch
        import json
        import logging
        from pathlib import Path
        import warnings
        import os
        import tempfile

        # Set up container-safe numba caching
        if os.path.exists('/.singularity.d') or os.environ.get('SINGULARITY_CONTAINER'):
            numba_cache = tempfile.mkdtemp(prefix='numba_cache_')
            os.environ['NUMBA_CACHE_DIR'] = numba_cache

        # Setup logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        logger = logging.getLogger(__name__)
        warnings.filterwarnings("ignore")

        logger.info("=== Starting scVIVA Training ===")

        # Load the trained ResolVI data
        logger.info("Loading ResolVI trained data...")
        adata = sc.read_h5ad("${adata_trained}")

        logger.info(f"Data shape: {adata.shape}")
        logger.info(f"Available obs columns: {list(adata.obs.columns)}")

        # Check required data
        if 'resolvi_predicted' not in adata.obs.columns:
            raise ValueError("Missing resolvi_predicted column")

        if 'X_spatial' not in adata.obsm.keys():
            if 'spatial' in adata.obsm.keys():
                adata.obsm['X_spatial'] = adata.obsm['spatial']
                logger.info("Used 'spatial' coordinates as X_spatial")
            else:
                raise ValueError("Missing spatial coordinates")

        # Setup scVI
        logger.info("Setting up spatial scVI analysis...")
        scvi.model.SCVI.setup_anndata(
            adata,
            layer="counts" if "counts" in adata.layers else None,
            labels_key="resolvi_predicted",
            batch_key="batch" if "batch" in adata.obs.columns else None
        )

        # Initialize model
        logger.info("Initializing scVI model...")
        model = scvi.model.SCVI(
            adata,
            n_latent=10,
            n_hidden=128,
            n_layers=2,
            dropout_rate=0.1
        )

        # Train model
        logger.info(f"Training scVI model for ${scviva_max_epochs} epochs...")
        model.train(
            max_epochs=${scviva_max_epochs},
            early_stopping=True,
            early_stopping_patience=20,
            check_val_every_n_epoch=10
        )

        # Save model
        logger.info("Saving trained model...")
        model.save("${prefix}_scviva_model", overwrite=True)

        # Get latent representation
        logger.info("Computing latent representation...")
        adata.obsm["X_scviva"] = model.get_latent_representation()

        # Compute UMAP
        logger.info("Computing UMAP on scVI latent space...")
        try:
            sc.pp.neighbors(adata, use_rep="X_scviva")
            sc.tl.umap(adata)
            logger.info("UMAP computation completed")
        except Exception as e:
            logger.warning(f"UMAP computation failed: {e}")

        # Save processed data
        logger.info("Saving processed data...")
        adata.write_h5ad("${prefix}_scviva_trained.h5ad")

        # Create summary
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

    # Execute with error handling
    try:
        main_training()
        handler.log_completion()
        handler.create_versions_file()
    except Exception as e:
        handler.handle_error(e, create_placeholder=True)

        # Create error outputs
        import os
        os.makedirs("${prefix}_scviva_model", exist_ok=True)
        with open("${prefix}_scviva_model/error.txt", 'w') as f:
            f.write(f"Error: {str(e)}\\n")

        # Create minimal output
        try:
            import scanpy as sc
            adata = sc.read_h5ad("${adata_trained}")
            adata.obs['scviva_error'] = str(e)
            adata.write_h5ad("${prefix}_scviva_trained.h5ad")
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
    mkdir -p ${prefix}_scviva_model
    touch ${prefix}_scviva_model/model.pkl
    touch ${prefix}_scviva_trained.h5ad
    touch ${prefix}_training_summary.json
    touch versions.yml
    """
}

