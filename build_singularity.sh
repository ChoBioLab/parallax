#!/bin/bash

# Build script for parallax Singularity container
# This script builds the Singularity container for the parallax pipeline

set -e

echo "Building parallax Singularity container..."

# Check if singularity is available
if ! command -v singularity &>/dev/null; then
	echo "ERROR: Singularity is not installed or not in PATH"
	echo "Please install Singularity first: https://sylabs.io/guides/3.0/user-guide/installation.html"
	exit 1
fi

# Check if we have sudo privileges (needed for building)
if ! sudo -n true 2>/dev/null; then
	echo "ERROR: This script requires sudo privileges to build the container"
	echo "Please run with sudo or ensure you have passwordless sudo configured"
	exit 1
fi

# Build the container
echo "Building container from parallax.def..."
sudo singularity build parallax.sif parallax.def

# Check if build was successful
if [ -f "parallax.sif" ]; then
	echo "✓ Container built successfully: parallax.sif"
	echo "✓ Container size: $(du -h parallax.sif | cut -f1)"

	# Test the container
	echo "Testing container..."
	singularity exec parallax.sif python -c "import torch; print(f'PyTorch: {torch.__version__}')"
	singularity exec parallax.sif python -c "import scvi; print(f'scvi-tools: {scvi.__version__}')"
	singularity exec parallax.sif python -c "import spatialdata; print(f'spatialdata: {spatialdata.__version__}')"

	echo "✓ Container test passed!"
	echo ""
	echo "You can now run the pipeline with:"
	echo "nextflow run main.nf -profile singularity --input samplesheet.csv --outdir results"
else
	echo "✗ Container build failed"
	exit 1
fi
