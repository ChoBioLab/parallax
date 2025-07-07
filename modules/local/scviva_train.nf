process SCVIVA_TRAIN {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'

    conda "/sc/arion/projects/untreatedIBD/ctastad/conda/envs/scvi"

    input:
    tuple val(meta), path(adata)  // From RESOLVI_TRAIN with predictions
    val max_epochs

    output:
    tuple val(meta), path("scviva_model/"), emit: model
    tuple val(meta), path("*_scviva_trained.h5ad"), emit: adata_processed
    tuple val(meta), path("scviva_training_log.txt"), emit: logs
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import scvi
    import pandas as pd
    import numpy as np
    import logging
    import sys

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('scviva_training_log.txt'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger = logging.getLogger(__name__)

    logger.info("=== Starting scVIVA Training ===")

    # Load data with ResolVI predictions
    adata = sc.read_h5ad("${adata}")
    logger.info(f"Loaded data with shape: {adata.shape}")

    # Use ResolVI predictions as cell type labels for scVIVA
    if "resolvi_predicted" not in adata.obs:
        raise ValueError("ResolVI predictions not found. Run ResolVI training first.")

    # Setup scVIVA with ResolVI predictions
    logger.info("Setting up scVIVA with ResolVI predictions...")

    # Compute environment features using preprocessing_anndata
    from scvi.external import SCVIVA

    # Setup scVIVA - it will compute environment features automatically
    SCVIVA.setup_anndata(
        adata,
        layer="counts",
        labels_key="resolvi_predicted",  # Use ResolVI predictions
        spatial_key="X_spatial"
    )

    # Initialize scVIVA model
    logger.info("Initializing scVIVA model...")
    model = SCVIVA(adata)

    # Train scVIVA model
    logger.info(f"Training scVIVA model for {${max_epochs}} epochs...")
    model.train(
        max_epochs=${max_epochs},
        early_stopping=True,
        early_stopping_patience=15,
        early_stopping_monitor="elbo_train"
    )

    # Save model
    logger.info("Saving scVIVA model...")
    model.save("scviva_model", overwrite=True)

    # Get latent representation
    logger.info("Computing scVIVA latent representation...")
    adata.obsm["X_scVIVA"] = model.get_latent_representation()

    # Compute neighborhood-aware UMAP
    logger.info("Computing neighborhood-aware UMAP...")
    sc.pp.neighbors(adata, use_rep="X_scVIVA", key_added="scviva")
    sc.tl.umap(adata, neighbors_key="scviva")
    adata.obsm["X_umap_scviva"] = adata.obsm["X_umap"].copy()

    # Save processed data
    adata.write_h5ad("${meta.id}_scviva_trained.h5ad")

    logger.info("=== scVIVA Training Completed ===")

    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    scvi-tools: {scvi.__version__}\\n')
    """
}
