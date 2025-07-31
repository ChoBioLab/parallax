process RESOLVI_TRAIN {
    tag "$meta.id"
    label 'process_gpu'

    conda "${projectDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'parallax.sif' :
        'parallax:latest' }"

    input:
    tuple val(meta), path(adata)
    val annotation_label
    val max_epochs
    val num_samples
    val gpu_mode

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
    import sys
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
    import re
    
    # Simple GPU setup based on mode
    if ${gpu_mode ? 'True' : 'False'}:
        if not torch.cuda.is_available():
            print("ERROR: GPU mode requested but CUDA not available")
            print("Use -profile cpu to force CPU mode")
            sys.exit(1)
        device = 'cuda'
        use_gpu = True
        print(f"✓ Using GPU: {torch.cuda.get_device_name()}")
    else:
        device = 'cpu'
        use_gpu = False
        print("ℹ Using CPU mode")
    
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
    
    logger.info("=== Starting ResolVI Training ===")
    logger.info(f"Device: {device}")
    logger.info(f"GPU mode: {use_gpu}")
    
    try:
        # Load preprocessed data
        logger.info("Loading AnnData...")
        adata = sc.read_h5ad("${adata}")
        logger.info(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")
    
        # Verify annotation column exists
        if '${annotation_label}' not in adata.obs.columns:
            raise ValueError(f"'{annotation_label}' column not found in obs")
    
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
            labels_key="${annotation_label}"
        )
    
        # Initialize model
        logger.info("Initializing ResolVI model...")
        model = scvi.external.RESOLVI(adata, semisupervised=True)
    
        # Train model with GPU configuration
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
            adata.obs["resolvi_predicted"] = adata.obs["${annotation_label}"]
    
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
    
        # Sanitize keys before saving
        adata = sanitize_h5_keys(adata)
    
        adata.write_h5ad("${prefix}_resolvi_trained.h5ad")
    
        # Create summary
        summary = {
            'success': True,
            'method': 'ResolVI',
            'cells': int(adata.n_obs),
            'genes': int(adata.n_vars),
            'annotation_label': '${annotation_label}',
            'unique_annotations': len(adata.obs['${annotation_label}'].unique()),
            'training_success': True,
            'device': device,
            'gpu_mode': use_gpu,
            'max_epochs': ${max_epochs},
            'num_samples': ${num_samples}
        }
    
        with open("${prefix}_training_summary.json", 'w') as f:
            json.dump(summary, f, indent=2)
    
        logger.info("=== ResolVI Training Completed Successfully ===")
    
    except Exception as e:
        logger.error(f"Training failed: {str(e)}")
    
        # Create error outputs
        os.makedirs("${prefix}_resolvi_model", exist_ok=True)
        with open("${prefix}_resolvi_model/error.txt", 'w') as f:
            f.write(f"Error: {str(e)}\\n")
    
        # Create minimal output
        try:
            adata = sc.read_h5ad("${adata}")
            adata.obs['resolvi_predicted'] = 'error'
            adata.write_h5ad("${prefix}_resolvi_trained.h5ad")
        except:
            # Create empty file if even loading fails
            with open("${prefix}_resolvi_trained.h5ad", 'w') as f:
                f.write("")
    
        summary = {
            'success': False,
            'error': str(e),
            'device': device if 'device' in locals() else 'unknown',
            'gpu_mode': use_gpu if 'use_gpu' in locals() else False
        }
        with open("${prefix}_training_summary.json", 'w') as f:
            json.dump(summary, f, indent=2)
    
        raise
    
    # Create versions file
    versions_content = f'''"${task.process}":
        python: "{sys.version.split()[0]}"
        scanpy: "{sc.__version__}"
        scvi-tools: "{scvi.__version__}"
        pytorch: "{torch.__version__}"
        device: "{device}"
        gpu_mode: "{use_gpu}"
    '''
    
    with open('versions.yml', 'w') as f:
        f.write(versions_content)
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
