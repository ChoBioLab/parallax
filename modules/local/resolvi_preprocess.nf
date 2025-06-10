process RESOLVI_PREPROCESS {
    tag "$meta.id"
    label 'process_medium'

    container 'oras://community.wave.seqera.io/library/pip_decoupler_scanpy_scvi-tools:8124a7e473830fad'

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("*_preprocessed.h5ad"), emit: adata
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python3 << 'EOF'
import scanpy as sc
import pandas as pd
import numpy as np
import scvi
import decoupler as dc

# Load data
adata = sc.read_h5ad("${adata}")

# Basic preprocessing
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# Save preprocessed data
adata.write("${prefix}_preprocessed.h5ad")

# Write versions
versions = {
    "scanpy": sc.__version__,
    "scvi-tools": scvi.__version__,
    "decoupler": dc.__version__
}

with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    for tool, version in versions.items():
        f.write(f'    {tool}: {version}\\n')
EOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_preprocessed.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: 1.9.6
        scvi-tools: 1.0.4
        decoupler: 1.4.0
    END_VERSIONS
    """
}

