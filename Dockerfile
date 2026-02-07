# syntax=docker/dockerfile:1

# Micromamba-based image with a single environment containing both Python and R.
# This keeps Python+R dependency management in one place (conda-forge).

FROM mambaorg/micromamba:1.5.8

# Auto-activate the target env in subsequent RUN commands
ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /app

# Add bioconda channel
RUN micromamba config append channels bioconda

# Install Python + R + core libs into the base env.
# Note: Signac is available on conda-forge for many platforms; if solving fails,
# you can drop r-signac here and install it in R instead.
RUN micromamba install -y -n base -c conda-forge \
    python=3.11 \
    pip \
    r-base \
    r-optparse \
    r-seurat \
    r-signac \
    bioconductor-genomicranges \
    numpy \
    pandas \
    scipy \
    anndata \
    scanpy \
    && micromamba clean -a -y

# Copy the project
COPY py/ /app/py/
COPY r/ /app/r/
COPY pyproject.toml uv.lock README.md /app/

# Ensure R CLIs are executable
RUN chmod +x /app/r/disc2r.R /app/r/cli.R || true

# Expose converters as first-class commands in the container.
# These wrappers ensure the micromamba base env is used even when not interactively activated.
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*
RUN printf '%s\n' \
    '#!/usr/bin/env bash' \
    'cd /app || exit 1' \
    'exec micromamba run -n base python -m py.cli "$@"' \
    > /usr/local/bin/py2disc \
    && chmod +x /usr/local/bin/py2disc \
    && ln -sf /usr/local/bin/py2disc /usr/local/bin/snaptosign \
    && printf '%s\n' \
    '#!/usr/bin/env bash' \
    'exec micromamba run -n base Rscript /app/r/disc2r.R "$@"' \
    > /usr/local/bin/disc2r \
    && chmod +x /usr/local/bin/disc2r
USER mambauser

# Optional: install additional Python deps that are not on conda-forge
# (e.g., snapatac2). This is best-effort; comment out if you prefer conda-only.
RUN pip install --no-cache-dir snapatac2 || true

# Use the micromamba entrypoint so `docker run IMAGE <cmd>` works normally.
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["bash"]
# (kept inside the image at /app/.venv)
