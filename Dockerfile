FROM nfcore/base:2.1
LABEL authors="Christopher Tastad" \
      description="Docker image containing all software requirements for the nf-core/parallax pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/parallax-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name parallax-1.0dev > parallax-1.0dev.yml

# Install additional Python packages not available in conda
RUN pip install --no-cache-dir \
    spatialdata-io==0.0.9 \
    spatialdata-plot==0.0.4

WORKDIR /data

