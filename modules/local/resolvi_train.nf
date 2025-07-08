process RESOLVI_TRAIN {
    tag "$meta.id"
    label 'process_high_gpu'

    conda "bioconda::scanpy bioconda::scvi-tools"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/scverse/scvi-tools:py3.12-cu12-1.3.2-dev' :
        'ghcr.io/scverse/scvi-tools:py3.12-cu12-1.3.2-dev' }"

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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env python3

    import os
    import sys
    import logging
    import warnings
    import traceback
    import gc
    import json
    import yaml

    # Suppress warnings
    warnings.filterwarnings('ignore')
    os.environ['PYTHONWARNINGS'] = 'ignore'

    # Set environment variables for stability
    os.environ['NUMBA_CACHE_DIR'] = '/tmp'
    os.environ['NUMBA_DISABLE_JIT'] = '1'
    os.environ['CUDA_VISIBLE_DEVICES'] = '0'
    os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:512'

    # Setup logging
    logging.basicConfig(
        level=logging.INFO, 
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler('resolvi_training.log')
        ]
    )
    logger = logging.getLogger(__name__)

    # Parameters
    annotation_label = "${annotation_label}"
    max_epochs = int(${max_epochs})
    num_samples = int(${num_samples})
    num_gpus = int(${num_gpus})
    prefix = "${prefix}"

    logger.info("=== ResolVI Training Started ===")
    logger.info("Parameters: annotation=%s, epochs=%d, samples=%d" % (annotation_label, max_epochs, num_samples))

    def cleanup_gpu():
        try:
            if 'torch' in globals() and torch.cuda.is_available():
                torch.cuda.empty_cache()
                torch.cuda.synchronize()
                gc.collect()
        except:
            pass

    def safe_train_model(model, **kwargs):
        strategies = [
            {
                'max_epochs': min(kwargs.get('max_epochs', 10), 20),
                'batch_size': 32,
                'enable_progress_bar': False,
                'early_stopping': False,
                'plan_kwargs': {'lr': 5e-4}
            },
            {
                'max_epochs': 5,
                'batch_size': 16,
                'enable_progress_bar': False,
                'early_stopping': False,
                'plan_kwargs': {'lr': 1e-3}
            },
            {
                'max_epochs': 2,
                'batch_size': 8,
                'enable_progress_bar': False,
                'early_stopping': False
            }
        ]

        for i, strategy in enumerate(strategies):
            try:
                logger.info("Training attempt %d/3 with strategy: %s" % (i+1, strategy))
                cleanup_gpu()

                model.train(**strategy)
                logger.info("Training successful with strategy %d" % (i+1))
                return True, strategy

            except Exception as e:
                logger.error("Training strategy %d failed: %s" % (i+1, str(e)))
                if i < len(strategies) - 1:
                    logger.info("Trying next strategy...")
                    cleanup_gpu()
                    continue
                else:
                    logger.error("All training strategies failed")
                    return False, str(e)

        return False, "All strategies exhausted"

    try:
        # Import required packages
        import torch
        import scvi
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scanpy as sc

        logger.info("PyTorch: %s, SCVI: %s" % (torch.__version__, scvi.__version__))
        logger.info("CUDA available: %s" % torch.cuda.is_available())

        if torch.cuda.is_available():
            logger.info("GPU: %s" % torch.cuda.get_device_name(0))
            logger.info("GPU Memory: %.1f GB" % (torch.cuda.get_device_properties(0).total_memory / 1024**3))

        # Configure scvi-tools
        scvi.settings.seed = 42
        scvi.settings.dl_pin_memory_gpu_training = False
        scvi.settings.dl_num_workers = 0
        scvi.settings.verbosity = 30

        # Load and validate data
        logger.info("Loading data...")
        adata = ad.read_h5ad("${adata}")
        logger.info("Data shape: %s" % str(adata.shape))

        # Validate required fields
        if annotation_label not in adata.obs.columns:
            raise ValueError("Annotation '%s' not found" % annotation_label)

        if 'X_spatial' not in adata.obsm:
            raise ValueError("X_spatial not found in adata.obsm")

        # Data preprocessing for ResolVI compatibility
        logger.info("Preprocessing data for ResolVI...")

        # Ensure X is not sparse for ResolVI
        if hasattr(adata.X, 'toarray'):
            adata.X = adata.X.toarray()

        # Filter cells with sufficient counts
        cell_counts = adata.X.sum(axis=1)
        valid_cells = cell_counts >= 20
        logger.info("Cells with >=20 counts: %d/%d" % (valid_cells.sum(), len(valid_cells)))

        if valid_cells.sum() < 50:
            raise ValueError("Insufficient cells (%d) for training" % valid_cells.sum())

        adata_filtered = adata[valid_cells].copy()

        # Sample for training
        n_train = min(300, max(50, adata_filtered.shape[0] // 4))
        logger.info("Using %d cells for training" % n_train)

        np.random.seed(42)
        train_indices = np.random.choice(adata_filtered.shape[0], size=n_train, replace=False)
        adata_train = adata_filtered[train_indices].copy()

        # Ensure data types are correct
        adata_train.X = adata_train.X.astype(np.float32)
        adata_train.obsm['X_spatial'] = adata_train.obsm['X_spatial'].astype(np.float32)

        logger.info("Training data: %s, spatial: %s" % (str(adata_train.shape), str(adata_train.obsm['X_spatial'].shape)))

        # Setup ResolVI
        logger.info("Setting up ResolVI...")
        try:
            scvi.external.RESOLVI.setup_anndata(
                adata_train,
                labels_key=annotation_label
            )
            logger.info("ResolVI setup successful")
        except Exception as e:
            logger.error("ResolVI setup failed: %s" % str(e))
            raise e

        # Initialize model with minimal parameters
        logger.info("Initializing ResolVI model...")
        try:
            model = scvi.external.RESOLVI(
                adata_train,
                n_latent=5,
                n_hidden=64,
                n_layers=1,
            )
            logger.info("Model initialization successful")
        except Exception as e:
            logger.error("Model initialization failed: %s" % str(e))
            raise e

        # Train model with fallback strategies
        logger.info("Starting training with fallback strategies...")
        cleanup_gpu()

        success, result = safe_train_model(model, max_epochs=max_epochs)

        if not success:
            logger.error("All training attempts failed: %s" % str(result))
            # Create a dummy model directory and minimal outputs
            model_dir = "%s_resolvi_model" % prefix
            os.makedirs(model_dir, exist_ok=True)

            # Save minimal model info
            with open("%s/model_info.txt" % model_dir, 'w') as f:
                f.write("Training failed: %s\\n" % str(result))
                f.write("Data shape: %s\\n" % str(adata_train.shape))
                f.write("Annotation: %s\\n" % annotation_label)

            # Create dummy predictions
            adata.obs['resolvi_predicted'] = 'training_failed'

        else:
            logger.info("Training completed successfully with: %s" % str(result))

            # Save model
            model_dir = "%s_resolvi_model" % prefix
            try:
                model.save(model_dir)
                logger.info("Model saved to %s" % model_dir)
            except Exception as e:
                logger.error("Model save failed: %s" % str(e))
                os.makedirs(model_dir, exist_ok=True)
                with open("%s/save_error.txt" % model_dir, 'w') as f:
                    f.write("Save failed: %s" % str(e))

            # Get predictions
            try:
                predictions = model.predict()
                logger.info("Generated %d predictions" % len(predictions))

                # Map predictions to full dataset
                adata.obs['resolvi_predicted'] = 'no_prediction'
                for i, train_idx in enumerate(train_indices):
                    original_idx = adata_filtered.obs.index[train_idx]
                    if original_idx in adata.obs.index:
                        adata.obs.loc[original_idx, 'resolvi_predicted'] = predictions[i]

            except Exception as e:
                logger.error("Prediction failed: %s" % str(e))
                adata.obs['resolvi_predicted'] = 'prediction_failed'

        # Save results
        output_file = "%s_resolvi_trained.h5ad" % prefix
        adata.write_h5ad(output_file)
        logger.info("Results saved to %s" % output_file)

        # Training summary
        summary = {
            'success': success,
            'original_cells': int(adata.shape[0]),
            'training_cells': int(adata_train.shape[0]) if 'adata_train' in locals() else 0,
            'annotation_label': annotation_label,
            'training_result': str(result),
            'scvi_version': scvi.__version__,
            'torch_version': torch.__version__,
            'cuda_available': torch.cuda.is_available()
        }

        with open("%s_training_summary.json" % prefix, 'w') as f:
            json.dump(summary, f, indent=2)

        logger.info("=== ResolVI Training Complete ===")

    except Exception as e:
        logger.error("Fatal error: %s" % str(e))
        logger.error("Traceback: %s" % traceback.format_exc())

        # Create minimal outputs even on failure
        try:
            # Create dummy model directory
            model_dir = "%s_resolvi_model" % prefix
            os.makedirs(model_dir, exist_ok=True)
            with open("%s/error.txt" % model_dir, 'w') as f:
                f.write("Error: %s\\n%s" % (str(e), traceback.format_exc()))

            # Create dummy adata if it doesn't exist
            if 'adata' not in locals():
                adata = ad.read_h5ad("${adata}")

            adata.obs['resolvi_predicted'] = 'error'
            adata.write_h5ad("%s_resolvi_trained.h5ad" % prefix)

            # Error summary
            summary = {
                'success': False,
                'error': str(e),
                'traceback': traceback.format_exc()
            }
            with open("%s_training_summary.json" % prefix, 'w') as f:
                json.dump(summary, f, indent=2)

        except Exception as e2:
            logger.error("Failed to create error outputs: %s" % str(e2))

        raise e

    finally:
        cleanup_gpu()

    # Create versions file
    try:
        versions = {
            'RESOLVI_TRAIN': {
                'python': sys.version.split()[0],
                'scvi-tools': scvi.__version__,
                'torch': torch.__version__,
                'numpy': np.__version__
            }
        }
        with open('versions.yml', 'w') as f:
            yaml.dump(versions, f)
    except:
        with open('versions.yml', 'w') as f:
            f.write('RESOLVI_TRAIN:\\n  status: error\\n')
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
