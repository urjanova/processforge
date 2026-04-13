# Use a lightweight Micromamba base
FROM mambaorg/micromamba:1.5.8

# Copy the environment file into the container
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Install dependencies and clean up to keep the image small
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Set the working directory
WORKDIR /app

# Copy the rest of your repository code
COPY --chown=$MAMBA_USER:$MAMBA_USER . .

# Install processforge and API dependencies from local source
RUN micromamba run -n base python -m pip install -e ".[api,coolprop]"

# Copy and register the OpenMC data-fetch startup script
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/ ./scripts/
RUN chmod +x ./scripts/fetch_openmc_data.sh

# Data volume for OpenMC cross sections and DAGMC geometry files.
#
# Docker:   docker run -v /host/openmc_data:/data ...
# Railway:  Attach a Railway Volume at /data in the service Volume settings.
#
# On first start, set OPENMC_DATA_URL (cross-section archive) and/or
# OPENMC_DAGMC_URL (geometry .h5m) to download data automatically.
# Subsequent starts skip the download (files already on the volume).
# OPENMC_CROSS_SECTIONS is exported automatically by the startup script.
VOLUME /data

# Default data root — override with -e OPENMC_DATA_ROOT=<path>
ENV OPENMC_DATA_ROOT=/data

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