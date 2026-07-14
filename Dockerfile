# Use a lightweight Micromamba base
FROM mambaorg/micromamba:1.5.8

# gcc is required by DOLFINx/FEniCSx for JIT-compiling variational forms at runtime
USER root
RUN apt-get update && apt-get install -y --no-install-recommends build-essential && \
    rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Copy the environment file into the container
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Install dependencies and clean up to keep the image small
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Set the working directory
WORKDIR /app

# Install processforge and API dependencies from local source.
# Copy package metadata first so dependency installation is cached
# across source-code changes.
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md ./
COPY --chown=$MAMBA_USER:$MAMBA_USER src/ src/
RUN micromamba run -n base python -m pip install -e ".[api,coolprop]"

# Copy the rest of the repository (application code, scripts, tests)
COPY --chown=$MAMBA_USER:$MAMBA_USER . .

# Copy and register the OpenMC data-fetch startup script
RUN chmod +x ./scripts/fetch_openmc_data.sh

# Data volume for OpenMC cross sections.
#
# Docker:   docker run -v /host/openmc_data:/data ...
# Railway:  Attach a Railway Volume at /data in the service Volume settings.
#
# On first start, set OPENMC_DATA_URL (cross-section archive) to download
# data automatically. Subsequent starts skip the download (files already
# on the volume). OPENMC_CROSS_SECTIONS is exported automatically by the
# startup script.
VOLUME /data

# Default data root — override with -e OPENMC_DATA_ROOT=<path>
ENV OPENMC_DATA_ROOT=/data

# Root for run outputs (results zarr + provider artifacts). Defaults to the
# /data volume so `docker run -v <host>:/data ...` captures outputs on the host.
# Override with -e PROCESSFORGE_OUTPUT_DIR=<path>.
ENV PROCESSFORGE_OUTPUT_DIR=/data

# Ensure the conda environment is activated for the entrypoint
ARG MAMBA_DOCKERFILE_COMMAND=activate
USER $MAMBA_USER

# Port used by the API server (pf-serve / python -m processforge.api.serve)
EXPOSE 9000

# Execute commands inside the activated conda environment. Default to CLI mode.
# The fetch_openmc_data.sh script runs first and downloads any missing data files.
# Examples:
#   docker run --rm <image>                        # runs processforge
#   docker run --rm -p 9000:9000 <image> pf-serve # runs API server
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/app/scripts/fetch_openmc_data.sh"]
CMD ["processforge"]