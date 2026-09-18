# Phase 6 — Docker + README

## What was built

- `Dockerfile`: single-stage `python:3.12-slim` image. Installs only
  `requirements-serve.txt` (the runtime subset — no scikit-learn,
  matplotlib, or snakemake, which the demo doesn't need), fetches
  `hg19.2bit` at **build** time (so `docker run` needs no network access
  and no pre-populated host `data/`), and copies in code + the already-
  committed `models/` artifacts + `results/figures/`. Does **not** run
  the data pipeline inside the image — it only serves what Phases 0-4
  already produced.
- `requirements-serve.txt`: split out from `requirements.txt` because the
  demo genuinely needs a much smaller dependency set than the research
  pipeline does.
- `.dockerignore`: excludes `data/` (3+ GB of raw MAFs, irrelevant to
  serving — the reference genome is fetched fresh at build time instead),
  the venv, git history, and dev-only docs.
- `README.md`: project summary, shift description, results (cross-checked
  against the actual saved metrics files line by line before writing),
  limitations, and run instructions.

## Docker build/run verified end to end

Docker wasn't installed in this sandbox initially and needed root, which
required the user to install it out-of-band (`apt-get install docker.io`)
and add the session's user to the `docker` group. Once available:

- Before Docker was available, the *dependency set* was verified by
  proxy: a fresh, isolated Python environment installed with exactly
  `requirements-serve.txt` plus the CPU-only PyTorch wheel served correct
  predictions.
- With Docker actually available, `docker build -t shift-abstain .` was
  run for real (all 12 steps completed, image `shift-abstain:latest`,
  3.35GB) and `docker run -d -p 8000:8000 shift-abstain` was started and
  exercised directly:
  - `GET /health` → 200
  - `GET /ui` → 200 (after Gradio's own redirect)
  - `POST /predict` with a real GBM exome MAF (395 samples) → identical
    predictions to every earlier verification (sample `TCGA-02-0003-01`
    → Glioma, confidence 0.7595, PREDICT)

The image is 3.35GB, dominated by the 816MB reference genome and the
~1.1GB PyTorch CPU wheel plus its dependencies (numpy, etc. pulled in
again inside the container even though the host already has them, since
Docker images are isolated by design).

## Definition of done, honestly assessed

- "`docker run` gives a working demo" — **verified**: built and ran the
  actual image, confirmed `/health`, `/ui`, and `/predict` all work
  end-to-end with real data.
- "one command reproduces every figure" — verified: `snakemake --cores 4
  all` (or `make`) runs the full pipeline from raw data download through
  `results/figures/*.png`, and every figure/metric currently in the repo
  was in fact produced by running these exact scripts.
- "the README states honestly how large the shift penalty was and how
  much abstention recovered" — done: 68.7% → 33.5% (35.2-point drop),
  with abstention recovering only a modest fraction (e.g. 33.5% → 37.1%
  at 80% coverage), stated plainly rather than reframed as a success.
