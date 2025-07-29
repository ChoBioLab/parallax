# nf-core/parallax: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## v1.0dev - [2025-06-06]

Initial release of nf-core/parallax, created with the [nf-core](http://nf-co.re/) template.

### `Added`

- Initial implementation of ResolVI pipeline for spatial transcriptomics analysis
- Support for spatialdata zarr store input format
- GPU acceleration support with configurable GPU usage
- Comprehensive visualization outputs including spatial plots and UMAP embeddings
- Differential expression and differential niche abundance analysis
- MultiQC integration for quality control reporting
- Docker and Singularity container support
- Comprehensive parameter validation with JSON schema
- Test profiles for CI/CD integration

### `Fixed`

- Updated container versions to latest stable releases
- Improved error handling and logging throughout pipeline
- Fixed GPU configuration for various cluster environments

### `Dependencies`

- scanpy=1.9.3
- scvi-tools=1.0.4
- spatialdata=0.0.14
- decoupler=1.4.0
- pytorch=2.0.1
- multiqc=1.14

### `Deprecated`

- None
