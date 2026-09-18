# SHIFT-ABSTAIN — tumor-origin classification that knows when to defer

## Context
Portfolio project for a computational oncology lab application (cancer genomics,
tumor evolution, glioma). Research area: selective prediction / abstention under
distribution shift. 2-week project, ~15 hours of the user's attention. Must end
in a deployed demo and an honest, reproducible result — not a leaderboard score.

## Objective
Train a tumor-type classifier on somatic mutation spectra from whole-exome data
(TCGA), then evaluate it on targeted-panel data (MSK-IMPACT). The panel-vs-exome
mismatch is a real covariate shift. Show that (a) accuracy degrades under that
shift, and (b) a calibrated abstention policy recovers most of the loss by
deferring a minority of cases.

## Data
- In-distribution: TCGA PanCancer Atlas somatic mutations (exome), from
  cBioPortal datahub / cbioportal.org/datasets.
- Shifted: an MSK-IMPACT targeted-panel cohort, same source.
- Phase 0 is a data spike: download the smallest viable files first, confirm
  required MAF columns exist, print per-cancer-type sample counts for both
  cohorts, and report back before writing any model code.
- If either cohort is not actually obtainable, STOP and say so. Never
  substitute simulated data.

## Features
Per sample: 96-channel SBS trinucleotide context counts, TMB, and a small set
of driver-gene mutation flags restricted to genes present in BOTH the exome
gene set and the panel gene list. Normalize counts by captured footprint
(per-Mb) so panel samples aren't trivially separable by depth alone. Document
the normalization choice in the README — it's the scientific crux.

## Models
Small: gradient-boosted trees + a 2-layer MLP. Uncertainty via deep ensemble
(5 seeds) plus MC-dropout for the MLP. Abstention baselines: softmax response,
predictive entropy, ensemble disagreement.

## Evaluation
Risk–coverage curves, AURC, selective accuracy at 80/90/95% coverage, ECE —
all reported separately for in-distribution and shifted test sets. Plain
accuracy-vs-coverage table too. Report negative results as-is.

## Deliverables
1. Reproducible pipeline (Snakemake preferred; Makefile acceptable fallback).
2. `results/` with figures and a metrics JSON.
3. FastAPI backend + Gradio UI: upload a MAF, return predicted type, calibrated
   confidence, ABSTAIN/PREDICT decision, and the risk–coverage plot.
4. Dockerfile that builds and runs the demo locally in one command.
5. README with the shift description, results table, and limitations section.

## Phases (stop for review after 0, 2, 4)
0. Data spike
1. Features
2. Baseline model — STOP for review
3. Uncertainty + abstention
4. Evaluation / figures — STOP for review
5. API / UI
6. Docker / README

Commit at the end of each phase.

## Hard rules
- Never fabricate numbers. Never fill a failed download with synthetic data.
  Never write a result into the README that hasn't actually been produced.
- Prefer boring, small dependencies. Ask before installing anything heavy.
- If a phase is ballooning past budget, propose a cut rather than silently
  expanding scope.

## Definition of done
`docker run` gives a working demo, one command reproduces every figure, and
the README states honestly how large the shift penalty was and how much
abstention recovered.
