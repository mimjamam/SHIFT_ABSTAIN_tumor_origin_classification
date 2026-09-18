# Phase 5 — FastAPI backend + Gradio UI

## What was built

- `scripts/compute_thresholds.py` — precomputes the `softmax_response`
  threshold (the best-ranking abstention score from Phase 4) at 80/90/95%
  coverage, separately for the exome and panel regimes, from the exact
  same Phase 4 test-set predictions. Writes `models/thresholds.json`.
- `api/inference.py` — loads the 5-member deep ensemble, scaler, class
  list, and thresholds once at process start, then turns an uploaded MAF
  into predictions using the *same* feature extraction and normalization
  code as training (`scripts/extract_features.py`,
  `scripts/build_dataset.py`'s footprint constants) — not a separate,
  simplified inference path that could silently drift from what was
  evaluated.
- `api/main.py` — FastAPI app:
  - `POST /predict` (multipart MAF upload + `assay_type` + `coverage`
    form fields) → predicted type, confidence, PREDICT/ABSTAIN decision,
    and top-3 probabilities for every sample in the file.
  - `GET /figures/{name}` — serves the Phase 4 evaluation figures.
  - `GET /health`.
  - A Gradio Blocks UI mounted at `/ui` (root `/` redirects there) so the
    whole demo runs as one process on one port:
    `uvicorn api.main:app`.

## Assay-type -> footprint/regime mapping (reuses Phase 1-4 constants exactly)

| Assay type | Footprint | Regime | Threshold source |
|---|---|---|---|
| Whole exome | 38 Mb | exome | Phase 4's `test_exome_id` score distribution |
| MSK-IMPACT IMPACT341 | 0.901587 Mb | panel | Phase 4's `test_panel_shifted` score distribution |
| MSK-IMPACT IMPACT410 | 1.021743 Mb | panel | Phase 4's `test_panel_shifted` score distribution |

Picking the wrong assay type for an upload would silently use the wrong
footprint (and the wrong calibration regime) — this mirrors a real
failure mode a deployed version of this tool would have to guard against,
not just a demo detail.

## Smoke-tested against real data

- `exome` / GBM TCGA MAF (395 samples): predictions returned for all 395,
  correctly favoring Glioma for most samples.
- `IMPACT341` / a slice of the real MSK-IMPACT MAF (6,687 samples):
  predictions returned for all, split 5,861 PREDICT / 826 ABSTAIN
  (~12.4% abstained) at the 90%-coverage operating point — in the
  expected ballpark for a threshold calibrated at 90% coverage on a
  similar but not identical sample set.
- Verified `/health`, `/ui` (200 after Gradio's own redirect), and
  `/predict` all work with the app running as a single `uvicorn` process.

## What's next (Phase 6)

Dockerfile so `docker run` gives a working demo in one command, and the
final README. The reference genome (`data/reference/hg19.2bit`, 816MB) is
needed at inference time for SBS-96 context lookups and isn't committed
to git — Phase 6 needs to decide how it gets into the container (baked
into the image vs. downloaded at build time vs. a mounted volume).
