process SCVIVA_ANALYZE {
    tag "$meta.id"
    label 'process_high'

    conda "/sc/arion/projects/untreatedIBD/ctastad/conda/envs/scvi"

    input:
    tuple val(meta), path(model_dir), path(adata)
    val comparisons  // Flexible comparison inputs

    output:
    tuple val(meta), path("*_scviva_analyzed.h5ad"), emit: adata_final
    tuple val(meta), path("scviva_de_results/"), emit: de_results
    tuple val(meta), path("*_scviva_analysis_report.txt"), emit: reports
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import scvi
    import pandas as pd
    import numpy as np
    import json
    import os
    import logging
    import matplotlib.pyplot as plt

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    logger.info("=== Starting scVIVA Niche-Aware DE Analysis ===")

    # Load data and model
    adata = sc.read_h5ad("${adata}")
    model = scvi.external.SCVIVA.load("${model_dir}", adata)

    # Create output directory
    os.makedirs("scviva_de_results", exist_ok=True)

    # Parse comparisons
    comparisons_input = "${comparisons}"
    comparisons = []

    if comparisons_input and comparisons_input != "null":
        try:
            if comparisons_input.endswith('.json'):
                with open(comparisons_input, 'r') as f:
                    comparisons = json.load(f)
            elif comparisons_input.endswith('.csv'):
                df = pd.read_csv(comparisons_input)
                comparisons = df.to_dict('records')
            else:
                comparisons = json.loads(comparisons_input)
        except Exception as e:
            logger.warning(f"Error parsing comparisons: {e}")
            comparisons = []

    # Default comparisons if none provided
    if not comparisons and "resolvi_predicted" in adata.obs:
        cell_types = adata.obs["resolvi_predicted"].unique().tolist()
        if len(cell_types) >= 2:
            comparisons = [
                {"group1": cell_types[0], "group2": cell_types[1], "name": f"{cell_types[0]}_vs_{cell_types[1]}"}
            ]

    # Perform niche-aware differential expression
    logger.info(f"Running {len(comparisons)} niche-aware DE comparisons...")

    for comp in comparisons:
        try:
            group1 = comp.get("group1")
            group2 = comp.get("group2")
            comp_name = comp.get("name", f"{group1}_vs_{group2}")

            logger.info(f"Niche-aware DE comparison: {comp_name}")

            # Perform niche-aware DE analysis using scVIVA
            de_results = model.differential_expression(
                adata,
                groupby="resolvi_predicted",
                group1=group1,
                group2=group2,
                delta=0.05,
                batch_size=None,
                all_stats=True,
                batch_correction=True
            )

            # Filter and identify marker genes
            de_results = de_results[
                (de_results['proba_not_de'] <= 0.05) &
                (de_results['bayes_factor'] >= 3) &
                (np.abs(de_results['lfc_mean']) >= 0.25)
            ]

            # Save results
            de_csv_path = f"scviva_de_results/scviva_de_{comp_name}.csv"
            de_results.to_csv(de_csv_path)
            logger.info(f"scVIVA DE results saved to {de_csv_path}")

            # Create volcano plot
            plt.figure(figsize=(10, 8))
            plt.scatter(de_results['lfc_mean'], -np.log10(de_results['proba_not_de']), alpha=0.6)
            plt.xlabel('Log2 Fold Change')
            plt.ylabel('-log10(P-value)')
            plt.title(f'Niche-aware DE: {comp_name}')
            plt.savefig(f"scviva_de_results/volcano_{comp_name}.png", dpi=300, bbox_inches='tight')
            plt.close()

        except Exception as e:
            logger.error(f"Error in scVIVA DE analysis for {comp_name}: {e}")
            with open(f"scviva_de_results/scviva_de_failed_{comp_name}.txt", "w") as f:
                f.write(f"scVIVA DE analysis failed: {e}")

    # Save analyzed data
    adata.write_h5ad("${meta.id}_scviva_analyzed.h5ad")

    # Generate report
    with open("${meta.id}_scviva_analysis_report.txt", "w") as f:
        f.write(f"scVIVA niche-aware DE analyses completed: {len(comparisons)}")

    logger.info("=== scVIVA Analysis Completed ===")

    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    scvi-tools: {scvi.__version__}\\n')
    """
}

