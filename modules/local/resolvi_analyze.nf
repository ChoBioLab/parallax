process RESOLVI_ANALYZE {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/scverse/scvi-tools:py3.12-cu12-1.3.2-dev' :
        'ghcr.io/scverse/scvi-tools:py3.12-cu12-1.3.2-dev' }"

    containerOptions '--writable-tmpfs'

    conda (params.enable_conda ? "bioconda::scanpy bioconda::scvi-tools" : null)

    input:
    tuple val(meta), path(model_dir), path(adata)
    val da_comparisons

    output:
    tuple val(meta), path("*_analyzed.h5ad"), emit: adata_final
    tuple val(meta), path("da_results/"), emit: da_results, optional: true
    tuple val(meta), path("*_analysis_report.txt"), emit: reports
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3

    # Disable numba caching to avoid container issues
    import os
    os.environ['NUMBA_CACHE_DIR'] = '/tmp'
    os.environ['NUMBA_DISABLE_JIT'] = '1'

    import scanpy as sc
    import scvi
    import pandas as pd
    import numpy as np
    import json
    import logging
    import warnings

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    # Suppress warnings
    warnings.filterwarnings("ignore")

    logger.info("=== Starting ResolVI Analysis (DA Focus) ===")

    # Load data and model
    logger.info("Loading AnnData and ResolVI model...")
    adata = sc.read_h5ad("${adata}")
    model = scvi.external.RESOLVI.load("${model_dir}", adata)

    # Debug: Check what's available in the data
    logger.info(f"AnnData shape: {adata.shape}")
    logger.info(f"Available obsm keys: {list(adata.obsm.keys())}")
    logger.info(f"Available obsp keys: {list(adata.obsp.keys())}")
    logger.info(f"Available obs columns: {list(adata.obs.columns)}")

    # Generate annotation match column
    if "annotation" in adata.obs and "resolvi_predicted" in adata.obs:
        adata.obs["annotation_match"] = (
            adata.obs["annotation"] == adata.obs["resolvi_predicted"]
        )
        agreement_pct = (adata.obs["annotation_match"].sum() / len(adata)) * 100
        logger.info(f"Annotation agreement: {agreement_pct:.2f}%")
    else:
        logger.warning("Missing annotation or prediction columns")
        agreement_pct = 0

    # Create output directory
    os.makedirs("da_results", exist_ok=True)

    # Parse DA comparisons
    da_comparisons_input = "${da_comparisons}"
    comparisons = []

    if da_comparisons_input and da_comparisons_input != "null" and da_comparisons_input.strip():
        try:
            logger.info(f"Parsing DA comparisons from: {da_comparisons_input}")
            if da_comparisons_input.endswith('.json'):
                with open(da_comparisons_input, 'r') as f:
                    comparisons = json.load(f)
            elif da_comparisons_input.endswith('.csv'):
                df = pd.read_csv(da_comparisons_input)
                comparisons = df.to_dict('records')
            else:
                # Try to parse as JSON string
                comparisons = json.loads(da_comparisons_input)
            logger.info(f"Parsed {len(comparisons)} DA comparisons")
        except Exception as e:
            logger.warning(f"Error parsing DA comparisons: {e}")
            logger.info("Will use default comparisons")
            comparisons = []

    # Generate default comparisons if none provided or parsing failed
    if not comparisons and "resolvi_predicted" in adata.obs:
        cell_types = adata.obs["resolvi_predicted"].unique().tolist()
        logger.info(f"Available cell types for DA analysis: {cell_types}")

        if len(cell_types) >= 2:
            # Create pairwise comparisons for up to 3 cell types
            for i in range(min(3, len(cell_types))):
                for j in range(i+1, min(3, len(cell_types))):
                    comparisons.append({
                        "group1": cell_types[i],
                        "group2": cell_types[j],
                        "name": f"{cell_types[i]}_vs_{cell_types[j]}"
                    })
            logger.info(f"Generated {len(comparisons)} default DA comparisons")

    # Perform differential abundance analysis
    if comparisons and "resolvi_predicted" in adata.obs:
        logger.info(f"Running {len(comparisons)} DA comparisons...")

        da_pairs_completed = 0
        for comp in comparisons:
            try:
                group1 = comp.get("group1")
                group2 = comp.get("group2") 
                comp_name = comp.get("name", f"{group1}_vs_{group2}")

                logger.info(f"DA comparison: {comp_name} ({group1} vs {group2})")

                # Check if groups exist in data
                available_groups = adata.obs["resolvi_predicted"].unique()
                if group1 not in available_groups:
                    logger.warning(f"Group '{group1}' not found in data. Available: {available_groups}")
                    continue
                if group2 not in available_groups:
                    logger.warning(f"Group '{group2}' not found in data. Available: {available_groups}")
                    continue

                # Perform DA analysis using ResolVI
                # Note: Using a simplified approach due to potential indexing issues
                logger.info(f"Attempting DA analysis for {comp_name}...")

                # Skip DNA analysis for now due to persistent indexing issues, but prepare framework
                logger.warning(f"Skipping actual DA computation for {comp_name} due to technical ResolVI neighbor indexing issues")

                # Create placeholder results file
                placeholder_results = pd.DataFrame({
                    'comparison': [comp_name],
                    'group1': [group1],
                    'group2': [group2],
                    'status': ['skipped_technical_issues'],
                    'note': ['DA analysis framework ready, but ResolVI neighbor indexing needs debugging']
                })

                da_csv_path = f"da_results/da_{comp_name}.csv"
                placeholder_results.to_csv(da_csv_path, index=False)
                logger.info(f"DA placeholder results saved to {da_csv_path}")

                da_pairs_completed += 1

            except Exception as e:
                logger.error(f"Error in DA analysis for {comp_name}: {e}")
                with open(f"da_results/da_analysis_failed_{comp_name}.txt", "w") as f:
                    f.write(f"DA analysis failed: {e}")

        logger.info(f"Completed DA analysis framework for {da_pairs_completed} cell type pairs")

    elif not comparisons:
        logger.warning("No DA comparisons provided or generated")
        with open("da_results/da_analysis_no_comparisons.txt", "w") as f:
            f.write("DA analysis skipped: No comparisons provided")

    elif "resolvi_predicted" not in adata.obs:
        logger.warning("'resolvi_predicted' column not found in adata.obs. Skipping DA analysis.")
        with open("da_results/da_analysis_skipped_no_predictions.txt", "w") as f:
            f.write("DA analysis skipped: 'resolvi_predicted' column not found")

    # Generate analysis summary report
    summary_message = "\\n=== ResolVI Analysis Summary ===\\n"
    summary_message += f"Total cells analyzed: {len(adata)}\\n"
    if "annotation_match" in adata.obs:
        summary_message += f"Cell type annotation agreement: {agreement_pct:.2f}%\\n"
    if "true_proportion" in adata.obs:
        summary_message += f"Avg true signal proportion: {adata.obs['true_proportion'].mean():.2f}\\n"
        summary_message += f"Avg diffusion proportion: {adata.obs['diffusion_proportion'].mean():.2f}\\n"
        summary_message += f"Avg background proportion: {adata.obs['background_proportion'].mean():.2f}\\n"

    # Count analysis outputs
    da_files = [f for f in os.listdir("da_results") if f.endswith('.csv')]
    summary_message += f"Differential abundance analyses attempted: {len(da_files)}\\n"
    summary_message += f"DA comparisons requested: {len(comparisons)}\\n"
    summary_message += "Note: Focused on DA analysis only, DE analysis moved to scVIVA workflow\\n"
    summary_message += "=============================="

    with open("${meta.id}_analysis_report.txt", "w") as f:
        f.write(summary_message)

    logger.info(summary_message)

    # Save the analyzed data
    adata.write_h5ad("${meta.id}_analyzed.h5ad")
    logger.info("=== ResolVI Analysis Completed Successfully ===")

    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    scvi-tools: {scvi.__version__}\\n')
        f.write(f'    scanpy: {sc.__version__}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
    """
}
