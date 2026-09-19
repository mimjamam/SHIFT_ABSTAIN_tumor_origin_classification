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
# Mirrored as a GitHub release asset on this repo rather than fetched from
# UCSC directly: UCSC's server returned a 257-byte stub instead of the real
# file on Render's build network (likely blocking/rate-limiting some cloud
# IP ranges), while GitHub's release CDN (Azure Blob-backed) works reliably.
# -f -L: fail on HTTP errors, follow the redirect to the actual blob URL.
# The size check is a second line of defense -- a truncated/corrupt
# download must fail the build loudly, never ship (this is exactly how the
# UCSC failure was caught instead of shipping a broken image).
RUN mkdir -p data/reference \
    && curl -f -L --retry 3 --retry-delay 5 -o data/reference/hg19.2bit \
       https://github.com/mimjamam/SHIFT_ABSTAIN_tumor_origin_classification/releases/download/reference-data/hg19.2bit \
    && actual_size=$(stat -c%s data/reference/hg19.2bit) \
    && echo "downloaded hg19.2bit: ${actual_size} bytes" \
    && [ "$actual_size" -gt 800000000 ]

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
