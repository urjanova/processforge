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
RUN micromamba run -n base python -m pip install -e ".[api]"

# Ensure the conda environment is activated for the entrypoint
ARG MAMBA_DOCKERFILE_COMMAND=activate
USER $MAMBA_USER

# Port used by the API server (pf-serve / python -m processforge.api.serve)
EXPOSE 9000

# Execute commands inside the activated conda environment. Default to CLI mode.
# Examples:
#   docker run --rm <image>                        # runs processforge
#   docker run --rm -p 9000:9000 <image> pf-serve # runs API server
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["processforge"]