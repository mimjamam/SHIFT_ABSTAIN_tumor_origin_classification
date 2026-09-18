# SHIFT-ABSTAIN

Tumor-origin classification from somatic mutation spectra, trained on
whole-exome data and evaluated under a real covariate shift to
targeted-panel sequencing — with a calibrated abstention policy that
knows when to defer instead of guessing.

This is a portfolio project, not a leaderboard entry. The goal was an
honest answer to two questions: how much does a real-world exome-to-panel
shift actually hurt a tumor-type classifier, and how much of that loss
can a selective-prediction policy recover. **The honest answer: shift
hurts a lot (68.7% → 33.5% accuracy), and abstention recovers only a
modest fraction of that, not most of it** — see [Results](#results).

## The shift

- **In-distribution**: [TCGA PanCancer Atlas](https://www.cbioportal.org/)
  somatic mutations, whole-exome sequencing (WES). 19 cancer types used
  here, ~38 Mb captured footprint.
- **Shifted**: [MSK-IMPACT](https://www.mskcc.org/msk-impact) targeted
  panel sequencing (341- or 410-gene versions), ~0.90–1.02 Mb captured
  footprint — about 1/40th the genomic territory.

Both cohorts were pulled from the public
[cBioPortal datahub](https://github.com/cBioPortal/datahub). Both are
GRCh37. The shift is a real deployment scenario: a model trained on
research-grade exome data being asked to work on the targeted panel data
that's actually used in clinics.

Features (per sample): 96-channel SBS trinucleotide mutational-context
rates, tumor mutational burden (TMB), and 30 pan-cancer driver-gene
mutation flags (restricted to genes present on the smaller MSK panel
version, so every flag means the same thing for every sample). All
count-based features are normalized by captured footprint (mutations per
Mb) — **this is the scientific crux of the project**: without it, a
model can "solve" panel-vs-exome by counting mutations instead of
learning tumor biology. Full derivation and caveats, including how the
panel footprint was reverse-engineered directly from MSK's own published
TMB values rather than approximated: [`docs/normalization.md`](docs/normalization.md).

14 cancer types are shared across both cohorts after taxonomy
reconciliation (NSCLC, Breast, Colorectal, Glioma, Prostate, Pancreatic,
Renal, Esophagogastric, Bladder, Melanoma, Thyroid, Ovarian, Endometrial,
Head & Neck) — **16,018 samples total: 8,215 exome / 7,803 panel**.

## Models

- Gradient-boosted trees (`sklearn.HistGradientBoostingClassifier`) and a
  2-layer MLP — small models, deliberately, per project scope.
- Uncertainty: a **5-seed deep ensemble** and **MC-dropout** (30
  stochastic passes), both on the MLP.
- Abstention scores: softmax response, predictive entropy, and ensemble
  disagreement (mutual information / BALD), each evaluated head-to-head.

Trained on TCGA (exome) only — 70/15/15 train/val/test, stratified.
MSK-IMPACT (panel) is held out **entirely** from training and threshold
selection, used only as the final shifted-evaluation set, matching a
realistic deployment where you don't get labeled shifted data to tune on.

## Results

Deep-ensemble MLP, full test sets (n=1,233 exome / n=7,803 panel):

| Coverage | In-distribution accuracy | Shifted accuracy |
|---|---|---|
| 100% (no abstention) | 68.7% | 33.5% |
| 90% | 72.5% | 35.0% |
| 80% | 77.6% | 37.1% |
| 50% | 91.1% | 43.3% |
| 10% | 99.2% | 56.4% |

**Abstention helps but does not close the shift gap.** Even discarding
90% of panel samples and keeping only the most-confident 10%, shifted
accuracy (56.4%) still falls short of the *unfiltered* in-distribution
accuracy (68.7%). This is a negative result relative to the project's
original hypothesis, reported as-is rather than reframed.

The mechanism is visible directly in calibration:

| Confidence variant | ID ECE | Shifted ECE |
|---|---|---|
| Single model | 0.033 | 0.404 |
| Deep ensemble | 0.035 | 0.299 |
| MC-dropout | 0.080 | 0.236 |

The model is close to well-calibrated in-distribution and becomes badly
**overconfident** under shift (predictions made at ~95% confidence on
panel data are right only ~51% of the time — see
`results/figures/calibration.png`). That's the root cause of the limited
abstention recovery: a confidence signal that has decoupled from actual
correctness can't rank samples for deferral as well as one that hasn't.

Among the three abstention scores, `softmax_response` and
`predictive_entropy` (closely tied) beat `ensemble_disagreement` on both
splits (AURC 0.114–0.117 vs 0.163 in-distribution; 0.547–0.554 vs 0.560
shifted) — a specific, useful negative result about which baseline earns
its complexity here.

Full numbers: `results/metrics_baseline.json`,
`results/metrics_uncertainty.json`, `results/metrics_evaluation.json`,
`results/accuracy_vs_coverage.csv`. Figures:
`results/figures/{risk_coverage_curves,calibration,aurc_comparison}.png`.
Per-phase writeups with more detail and caveats: `docs/phase*.md`.

## Limitations

- **Panel footprint precision**: the MSK-IMPACT 341/410-gene footprints
  (0.901587 / 1.021743 Mb) were reverse-engineered from MSK's own
  published per-sample TMB values, not independently verified against a
  BED file of the actual capture intervals. See `docs/normalization.md`.
- **Two shifts are conflated**: TCGA MAFs report silent/intronic/UTR
  variants; the MSK-IMPACT MAF reports almost none (clinical panels
  generally don't report synonymous variants). SBS-96 features are built
  from *all* reported SNVs in each file, so the reported shift penalty
  reflects both the real capture-footprint difference *and* this
  pipeline-reporting difference — they cannot be fully disentangled from
  public MAFs alone. See `docs/normalization.md`.
- **Cohort composition differs for real, non-technical reasons**: MSK-IMPACT
  testing is disproportionately ordered for metastatic/refractory
  patients, who carry higher mutation burden; TCGA is a broader surgical
  cohort. Median panel TMB is *higher* than median exome TMB despite the
  panel covering far fewer genes — a real selection effect, not a
  pipeline artifact, that's folded into the reported shift.
- **14-class taxonomy is a simplification** of what TCGA's 32 studies and
  MSK's 58 cancer-type labels actually distinguish (e.g. TCGA keeps
  LIHC/CHOL separate where MSK's "Hepatobiliary Cancer" combines them).
- **Small models by design** (project scope), on 127 features. A larger
  feature set (e.g. copy-number, structural variants, or the full
  driver-gene panel rather than the 30-gene subset used here) would
  likely improve both accuracy and calibration, in and out of
  distribution.
- The Docker image is 3.35GB, dominated by the 816MB reference genome and
  the CPU-only PyTorch install — fine for a local demo, not optimized for
  distribution.

## Running it

### Reproduce the pipeline

```
pip install -r requirements.txt
snakemake --cores 4 all        # or: make
```

This downloads both cohorts (~3.3 GB of MAFs + 816 MB reference genome),
extracts features, trains all models, and regenerates every figure and
metrics file in `results/`.

### Run the demo

```
pip install -r requirements-serve.txt
pip install torch --index-url https://download.pytorch.org/whl/cpu
cd api && uvicorn main:app --host 0.0.0.0 --port 8000
```

or, with Docker:

```
docker build -t shift-abstain .
docker run -p 8000:8000 shift-abstain
```

Then open `http://localhost:8000` (redirects to the Gradio UI at `/ui`),
or `POST` a MAF file to `http://localhost:8000/predict` directly
(`multipart/form-data`: `maf_file`, `assay_type` ∈
`{exome, IMPACT341, IMPACT410}`, `coverage` ∈ `{80, 90, 95}`).

## Repository layout

```
scripts/          feature extraction, training, evaluation (the pipeline)
api/               FastAPI backend + Gradio UI
models/            trained model weights, scaler, thresholds (committed -- small, deterministic)
results/           metrics JSON, accuracy-vs-coverage table, figures (committed)
data/              raw + processed data (gitignored -- regenerate via the pipeline)
docs/              one write-up per project phase, with the reasoning and caveats behind each decision
Snakefile, Makefile  pipeline entry points
Dockerfile         one-command demo
```
