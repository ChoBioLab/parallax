process RESOLVI_VISUALIZE {
    tag "$meta.id"
    label 'process_medium'

    conda "/sc/arion/projects/untreatedIBD/ctastad/conda/envs/scvi"

    input:
    tuple val(meta), path(adata)
    val marker_genes

    output:
    tuple val(meta), path("plots/*"), emit: plots
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genes_arg = marker_genes ? marker_genes.join(',') : ""

    """
    python3 << 'EOF'
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import decoupler as dc
import os
import logging
import warnings

# Setup
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore")
plt.ioff()
os.makedirs("plots", exist_ok=True)

# Load data
adata = sc.read_h5ad("${meta.id}_analyzed.h5ad")

# Generate basic plots
try:
    if "X_umap" in adata.obsm:
        sc.pl.umap(adata, color="annotation", show=False)
        plt.savefig("plots/umap_annotation.png", dpi=300, bbox_inches="tight")
        plt.close()
        logger.info("Generated UMAP plot")
except Exception as e:
    logger.error(f"Error generating UMAP: {e}")

try:
    if "X_spatial" in adata.obsm:
        coords = adata.obsm["X_spatial"]
        plt.figure(figsize=(10, 8))
        plt.scatter(coords[:, 0], coords[:, 1], s=1, alpha=0.7)
        plt.title("Spatial Distribution")
        plt.savefig("plots/spatial_overview.png", dpi=300, bbox_inches="tight")
        plt.close()
        logger.info("Generated spatial plot")
except Exception as e:
    logger.error(f"Error generating spatial plot: {e}")

logger.info("Visualization completed!")
EOF

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    seaborn: \$(python -c "import seaborn; print(seaborn.__version__)")
    pandas: \$(python -c "import pandas; print(pandas.__version__)")
    numpy: \$(python -c "import numpy; print(numpy.__version__)")
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p plots
    touch plots/umap_annotation.png
    touch plots/spatial_overview.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: 1.9.6
        matplotlib: 3.7.1
        decoupler: 1.4.0
    END_VERSIONS
    """
}

