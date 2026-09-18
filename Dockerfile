FROM python:3.12-slim

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends curl \
    && rm -rf /var/lib/apt/lists/*

COPY requirements-serve.txt .
RUN pip install --no-cache-dir -r requirements-serve.txt \
    && pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu

# Reference genome for SBS-96 trinucleotide context lookups at inference
# time (see docs/phase1_features.md). Fetched at build time so `docker run`
# needs no network access and no pre-populated host data/ directory.
RUN mkdir -p data/reference \
    && curl -s -o data/reference/hg19.2bit \
       https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

# Code + the already-trained model artifacts (committed to the repo --
# small, deterministic, see docs/phase3_uncertainty.md). The full data
# pipeline (scripts/build_dataset.py etc.) is NOT run inside this image;
# this image only serves what Phases 0-4 already produced.
COPY scripts/ scripts/
COPY api/ api/
COPY models/ models/
COPY results/figures/ results/figures/

EXPOSE 8000
# PORT is injected by some hosting platforms (e.g. Render) to tell the
# container which port to listen on; default to 8000 for local `docker run`.
CMD ["sh", "-c", "cd api && uvicorn main:app --host 0.0.0.0 --port ${PORT:-8000}"]
