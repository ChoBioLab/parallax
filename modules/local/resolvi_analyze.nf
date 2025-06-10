process RESOLVI_ANALYZE {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'

    conda "bioconda::scvi-tools bioconda::decoupler"
    container 'oras://community.wave.seqera.io/library/pip_decoupler_scanpy_scvi-tools:8124a7e473830fad'

    input:
    tuple val(meta), path(model_dir), path(adata)

    output:
    tuple val(meta), path("*_analyzed.h5ad"), emit: adata_final
    tuple val(meta), path("de_results/"), emit: de_results, optional: true
    tuple val(meta), path("da_results/"), emit: da_results, optional: true
    tuple val(meta), path("*_analysis_report.txt"), emit: reports
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3

    import scanpy as sc
    import scvi
    import pandas as pd
    import numpy as np
    import decoupler as dc
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os
    import logging
    import warnings

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    # Suppress warnings
    warnings.filterwarnings("ignore")

    # Load data and model
    logger.info("Loading AnnData and ResolVI model...")
    adata = sc.read_h5ad("${adata}")
    model = scvi.external.RESOLVI.load("${model_dir}", adata)

    # Generate annotation match column
    adata.obs["annotation_match"] = (
        adata.obs["annotation"] == adata.obs["resolvi_predicted"]
    )
    agreement_pct = (adata.obs["annotation_match"].sum() / len(adata)) * 100
    logger.info(f"Annotation agreement: {agreement_pct:.2f}%")

    # Setup spatial neighbors for differential niche abundance
    logger.info("Computing spatial neighbors for DNA analysis...")
    try:
        sc.pp.neighbors(
            adata,
            use_rep="X_spatial",
            n_neighbors=15,
            key_added="spatial",
        )
        logger.info("Spatial neighbors computed successfully.")
    except Exception as e:
        logger.error(f"Error computing spatial neighbors: {e}")

    # Differential expression analysis
    cell_types = adata.obs["resolvi_predicted"].unique().tolist()
    logger.info(f"Predicted cell types for DE analysis: {cell_types}")

    if len(cell_types) >= 2:
        os.makedirs("de_results", exist_ok=True)
        cell_types_to_compare_de = cell_types[:2]
        logger.info(f"Comparing cell types for DE: {cell_types_to_compare_de[0]} vs {cell_types_to_compare_de[1]}")

        try:
            de_results = model.differential_expression(
                adata,
                groupby="resolvi_predicted",
                group1=cell_types_to_compare_de[0],
                group2=cell_types_to_compare_de[1],
                weights="importance",
                pseudocounts=1e-2,
                delta=0.05,
                filter_outlier_cells=True,
                mode="change",
                test_mode="three",
            )
            de_csv_path = f"de_results/de_{cell_types_to_compare_de[0]}_vs_{cell_types_to_compare_de[1]}.csv"
            de_results.to_csv(de_csv_path)
            logger.info(f"DE results saved to {de_csv_path}")

            # Generate volcano plot
            plt.figure(figsize=(10, 8))
            dc.plot_volcano_df(
                de_results,
                x="lfc_mean",
                y="proba_not_de",
                sign_thr=0.1,
                lFCs_thr=0.4,
                top=30,
                figsize=(10, 8),
            )
            plt.title(f"DE Genes: {cell_types_to_compare_de[0]} vs {cell_types_to_compare_de[1]}")
            de_plot_path = f"de_results/de_volcano_{cell_types_to_compare_de[0]}_vs_{cell_types_to_compare_de[1]}.png"
            plt.savefig(de_plot_path, dpi=300, bbox_inches="tight")
            plt.close()
            logger.info(f"DE volcano plot saved to {de_plot_path}")

        except Exception as e:
            logger.error(f"Error during DE analysis: {e}")
    else:
        logger.warning("Not enough unique cell types (need >= 2) for DE. Skipping.")

    # Differential niche abundance analysis
    logger.info("Performing differential niche abundance (DNA) analysis...")

    if len(cell_types) >= 2 and "spatial_connectivities" in adata.obsp:
        os.makedirs("da_results", exist_ok=True)
        cell_types_to_compare_da = cell_types[:2]
        logger.info(f"Comparing cell types for DNA: {cell_types_to_compare_da[0]} vs {cell_types_to_compare_da[1]}")

        try:
            da_results = model.differential_niche_abundance(
                adata,
                groupby="resolvi_predicted",
                group1=cell_types_to_compare_da[0],
                group2=cell_types_to_compare_da[1],
                neighbor_key="spatial_connectivities",
                test_mode="three",
                delta=0.05,
                pseudocounts=3e-2,
            )
            da_csv_path = f"da_results/da_{cell_types_to_compare_da[0]}_vs_{cell_types_to_compare_da[1]}.csv"
            da_results.to_csv(da_csv_path)
            logger.info(f"DNA results saved to {da_csv_path}")

            # Generate volcano plot
            plt.figure(figsize=(10, 8))
            dc.plot_volcano_df(
                da_results,
                x="lfc_mean",
                y="proba_not_de",
                sign_thr=0.1,
                lFCs_thr=0.5,
                top=10,
                figsize=(10, 8),
            )
            plt.title(f"DNA: {cell_types_to_compare_da[0]} vs {cell_types_to_compare_da[1]}")
            da_plot_path = f"da_results/da_volcano_{cell_types_to_compare_da[0]}_vs_{cell_types_to_compare_da[1]}.png"
            plt.savefig(da_plot_path, dpi=300, bbox_inches="tight")
            plt.close()
            logger.info(f"DNA volcano plot saved to {da_plot_path}")

        except Exception as e:
            logger.error(f"Error in DNA analysis: {e}")
    elif "spatial_connectivities" not in adata.obsp:
        logger.warning("Spatial connectivities not found. Skipping DNA analysis.")
    else:
        logger.warning("Not enough unique cell types (need >= 2) for DNA. Skipping.")

    # Generate analysis summary report
    summary_message = "\\n=== ResolVI Analysis Summary ===\\n"
    summary_message += f"Total cells analyzed: {len(adata)}\\n"
    if "annotation_match" in adata.obs:
        summary_message += f"Cell type annotation agreement: {agreement_pct:.2f}%\\n"
    if "true_proportion" in adata.obs:
        summary_message += f"Avg true signal proportion: {adata.obs['true_proportion'].mean():.2f}\\n"
        summary_message += f"Avg diffusion proportion: {adata.obs['diffusion_proportion'].mean():.2f}\\n"
        summary_message += f"Avg background proportion: {adata.obs['background_proportion'].mean():.2f}\\n"
    else:
        summary_message += "Noise proportion data not available.\\n"
    summary_message += "=============================="

    with open("${meta.id}_analysis_report.txt", "w") as f:
        f.write(summary_message)

    logger.info(summary_message)

    # Save the analyzed data
    adata.write_h5ad("${meta.id}_analyzed.h5ad")
    logger.info("Analysis completed successfully!")

    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    scvi-tools: {scvi.__version__}\\n')
        f.write(f'    decoupler: {dc.__version__}\\n')
        f.write(f'    scanpy: {sc.__version__}\\n')
    """
}
