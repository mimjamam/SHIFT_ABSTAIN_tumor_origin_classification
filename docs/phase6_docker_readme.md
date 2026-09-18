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

## Important honesty note: Docker was not runnable in this environment

There is no Docker (or Podman) installed in this sandbox, and installing
it requires root, which wasn't available. **The `docker build`/`docker
run` path itself could not be executed and verified end-to-end here.**

What *was* verified, as the closest available substitute: a fresh,
isolated Python environment was built with `uv venv`, installed with
*exactly* `requirements-serve.txt` plus the CPU-only PyTorch wheel (the
same two install steps the Dockerfile runs), and the API was started and
exercised against real data from that environment alone --

- `GET /health` → 200
- `GET /ui` → 200 (after Gradio's own redirect)
- `POST /predict` with a real GBM exome MAF (395 samples) → identical
  predictions to the full dev environment (e.g. sample
  `TCGA-02-0003-01` → Glioma, confidence 0.7595, PREDICT)

This confirms the *dependency set* the Dockerfile installs is complete
and sufficient to serve correct predictions. It does **not** confirm the
Dockerfile's build steps themselves (base image, `apt-get`, `COPY` paths,
the `hg19.2bit` curl step, the final `CMD`) are free of typos or
environment-specific issues, since those can only be caught by an actual
`docker build`. **Before relying on this for a real demo, run `docker
build -t shift-abstain .` and `docker run -p 8000:8000 shift-abstain`
yourself and confirm `http://localhost:8000` comes up.**

## Definition of done, honestly assessed

- "`docker run` gives a working demo" — Dockerfile written and its
  dependency set verified by proxy (above), but the actual `docker
  build`/`run` commands are **unverified** in this environment.
- "one command reproduces every figure" — verified: `snakemake --cores 4
  all` (or `make`) runs the full pipeline from raw data download through
  `results/figures/*.png`, and every figure/metric currently in the repo
  was in fact produced by running these exact scripts.
- "the README states honestly how large the shift penalty was and how
  much abstention recovered" — done: 68.7% → 33.5% (35.2-point drop),
  with abstention recovering only a modest fraction (e.g. 33.5% → 37.1%
  at 80% coverage), stated plainly rather than reframed as a success.
