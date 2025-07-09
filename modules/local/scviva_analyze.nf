process SCVIVA_ANALYZE {
    tag "$meta.id"
    label 'process_medium'

    conda "/sc/arion/projects/untreatedIBD/ctastad/conda/envs/scvi"

    input:
    tuple val(meta), path(scviva_model_dir)
    tuple val(meta2), path(adata_trained)
    val scviva_comparisons

    output:
    tuple val(meta), path("*_scviva_results.h5ad"), emit: results
    tuple val(meta), path("*_scviva_de_results.csv"), emit: de_results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash
    set -euo pipefail

    # Explicitly activate conda environment
    source /hpc/users/tastac01/micromamba/etc/profile.d/conda.sh
    conda activate /sc/arion/projects/untreatedIBD/ctastad/conda/envs/scvi

    python << 'EOF'
import os
os.environ['NUMBA_CACHE_DIR'] = '/tmp'
os.environ['NUMBA_DISABLE_JIT'] = '1'
os.environ['MPLBACKEND'] = 'Agg'

import scanpy as sc
import pandas as pd
import numpy as np
import scvi
import torch
import sys
import json
import logging
from pathlib import Path
import warnings

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

warnings.filterwarnings("ignore")

logger.info("=== Starting scVIVA Analysis ===")

try:
    # Load the scVIVA trained data
    logger.info("Loading scVIVA trained data...")
    adata = sc.read_h5ad("${adata_trained}")

    logger.info(f"Data shape: {adata.shape}")
    logger.info(f"Available obs columns: {list(adata.obs.columns)}")
    logger.info(f"Available obsm keys: {list(adata.obsm.keys())}")

    # Check if we have the required data
    if 'resolvi_predicted' not in adata.obs.columns:
        logger.error("'resolvi_predicted' column not found in adata.obs")
        raise ValueError("Missing resolvi_predicted column")

    # Load the trained scVI model
    logger.info("Loading scVI model...")
    model = scvi.model.SCVI.load("${scviva_model_dir}", adata)

    # Handle comparisons parameter
    comparisons_input = "${scviva_comparisons}"
    comparisons = []

    if comparisons_input and comparisons_input != "default_pairwise" and os.path.exists(comparisons_input):
        logger.info(f"Loading comparisons from: {comparisons_input}")
        if comparisons_input.endswith('.json'):
            with open(comparisons_input, 'r') as f:
                comparisons_data = json.load(f)
                if isinstance(comparisons_data, list):
                    comparisons = comparisons_data
                elif 'comparisons' in comparisons_data:
                    comparisons = comparisons_data['comparisons']
        elif comparisons_input.endswith('.csv'):
            comp_df = pd.read_csv(comparisons_input)
            comparisons = [
                {
                    'group1': row['group1'],
                    'group2': row['group2'],
                    'covariate': row.get('covariate', None)
                }
                for _, row in comp_df.iterrows()
            ]
    else:
        # Generate default pairwise comparisons
        logger.info("Generating default pairwise comparisons...")
        cell_types = adata.obs['resolvi_predicted'].unique()
        cell_types = [ct for ct in cell_types if ct != 'error' and not pd.isna(ct)]

        for i, ct1 in enumerate(cell_types):
            for ct2 in cell_types[i+1:]:
                comparisons.append({
                    'group1': ct1,
                    'group2': ct2,
                    'covariate': None
                })

    logger.info(f"Running {len(comparisons)} comparisons...")

    # Perform differential expression analysis
    de_results_list = []

    for comp in comparisons:
        logger.info(f"Running comparison: {comp['group1']} vs {comp['group2']}")

        try:
            # Create binary labels for this comparison
            cell_type_mask = adata.obs['resolvi_predicted'].isin([comp['group1'], comp['group2']])

            if cell_type_mask.sum() == 0:
                logger.warning(f"No cells found for comparison {comp['group1']} vs {comp['group2']}")
                continue

            # Temporarily add comparison labels
            adata.obs['temp_de_labels'] = 'other'
            adata.obs.loc[adata.obs['resolvi_predicted'] == comp['group1'], 'temp_de_labels'] = comp['group1']
            adata.obs.loc[adata.obs['resolvi_predicted'] == comp['group2'], 'temp_de_labels'] = comp['group2']

            # Filter to only cells in this comparison
            adata_sub = adata[cell_type_mask].copy()

            if adata_sub.n_obs < 10:
                logger.warning(f"Too few cells for comparison {comp['group1']} vs {comp['group2']} ({adata_sub.n_obs} cells)")
                continue

            # Perform DE analysis
            de_df = model.differential_expression(
                adata_sub,
                groupby='temp_de_labels',
                group1=comp['group1'],
                group2=comp['group2'],
                delta=0.25,
                batch_size=None
            )

            # Add comparison information
            de_df['comparison'] = f"{comp['group1']}_vs_{comp['group2']}"
            de_df['group1'] = comp['group1']
            de_df['group2'] = comp['group2']
            de_df['gene'] = de_df.index

            de_results_list.append(de_df)
            logger.info(f"Found {len(de_df)} DE genes for {comp['group1']} vs {comp['group2']}")

        except Exception as e:
            logger.warning(f"DE analysis failed for {comp['group1']} vs {comp['group2']}: {e}")
            continue

    # Clean up temporary column
    if 'temp_de_labels' in adata.obs.columns:
        del adata.obs['temp_de_labels']

    # Combine all DE results
    if de_results_list:
        all_de_results = pd.concat(de_results_list, ignore_index=True)

        # Sort by significance and effect size
        all_de_results = all_de_results.sort_values(['comparison', 'bayes_factor'], ascending=[True, False])

        # Save DE results
        de_output_file = f"${prefix}_scviva_de_results.csv"
        all_de_results.to_csv(de_output_file, index=False)
        logger.info(f"Saved {len(all_de_results)} DE results to {de_output_file}")
    else:
        logger.warning("No DE results generated - creating empty file")
        empty_de = pd.DataFrame(columns=['gene', 'comparison', 'group1', 'group2', 'lfc_mean', 'lfc_std', 'bayes_factor', 'prob_de'])
        empty_de.to_csv(f"${prefix}_scviva_de_results.csv", index=False)

    # Save the results data
    results_output_file = f"${prefix}_scviva_results.h5ad"
    adata.write_h5ad(results_output_file)
    logger.info(f"Results data saved to: {results_output_file}")

    logger.info("=== scVIVA Analysis Completed Successfully ===")

except Exception as e:
    logger.error(f"scVIVA analysis failed: {e}")
    import traceback
    logger.error(f"Traceback: {traceback.format_exc()}")

    # Create empty DE results on error
    empty_de = pd.DataFrame(columns=['gene', 'comparison', 'group1', 'group2', 'lfc_mean', 'lfc_std', 'bayes_factor', 'prob_de'])
    empty_de.to_csv(f"${prefix}_scviva_de_results.csv", index=False)

    # Create minimal results output
    if 'adata' in locals():
        adata.obs['analysis_error'] = str(e)
        adata.write_h5ad(f"${prefix}_scviva_results.h5ad")

    raise e

# Write versions
with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f'    scvi-tools: {scvi.__version__}\\n')
    f.write(f'    torch: {torch.__version__}\\n')
    f.write(f'    scanpy: {sc.__version__}\\n')
EOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_scviva_results.h5ad
    touch ${prefix}_scviva_de_results.csv
    touch versions.yml
    """
}

