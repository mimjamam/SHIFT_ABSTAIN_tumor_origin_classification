# Phase 1 — feature pipeline

## What was built

- `scripts/sbs96.py` — canonical COSMIC SBS-96 trinucleotide category list
  and the pyrimidine-normalization logic.
- `scripts/extract_features.py` — per-MAF feature extraction: raw SBS-96
  counts (via `hg19.2bit` context lookups with `py2bit`), nonsynonymous
  mutation count, and driver-gene flags.
- `scripts/build_dataset.py` — orchestrates extraction across all 19 TCGA
  studies + MSK-IMPACT, applies per-Mb footprint normalization (see
  `docs/normalization.md`), maps both cohorts onto the shared 14-class
  taxonomy, and writes `data/processed/features.parquet`.
- `scripts/fetch_gene_panels.py`, `scripts/fetch_cbioportal_lfs.py`,
  `scripts/download_tcga_mafs.sh` — reproducible data acquisition.
- `Snakefile` / `Makefile` — pipeline entry points (`snakemake --cores 4
  all` or `make`).

## Driver genes

30 well-established pan-cancer driver genes (TP53, KRAS, PIK3CA, PTEN,
APC, EGFR, BRAF, IDH1, IDH2, VHL, RB1, NF1, ARID1A, CDKN2A, SMAD4, STK11,
ATM, BRCA1, BRCA2, CTNNB1, KMT2D, NOTCH1, FBXW7, MYC, ERBB2, MET, KIT,
CDH1, GATA3, RUNX1), restricted to genes present in **IMPACT341** (the
smaller of the two MSK panel versions in this cohort, and a strict subset
of IMPACT410 — confirmed by direct set comparison) so every flag is valid
for every MSK-IMPACT sample regardless of which panel version it used. All
30 candidates happened to already be in IMPACT341, so none were dropped.

## Correctness checks performed

- SBS-96 context lookups: 0 reference-mismatch errors across all
  20 files (19 TCGA + MSK-IMPACT), ~3.1M SNVs total — every context3 lookup
  from `hg19.2bit` matched the MAF's own `Reference_Allele`, confirming
  build/coordinate consistency.
- TMB cross-validation: our nonsynonymous-count-based TMB correlates with
  MSK's own precomputed `TMB_NONSYNONYMOUS` at Pearson r=0.9987 (see
  `docs/normalization.md` for how this also let us recover MSK's exact
  per-panel-version footprint).
- Driver-gene flag rates are in the range expected from the cancer
  genomics literature (TP53 ~42-49%, PIK3CA ~15%, IDH1 higher in the
  exome cohort where glioma is proportionally better represented).

## Final dataset

16,018 samples total (8,215 exome / 7,803 panel) across 14 shared cancer
types. Per-class counts:

| Class | Exome (TCGA) | Panel (MSK-IMPACT) |
|---|---|---|
| NSCLC | 1,031 | 1,620 |
| Breast | 1,009 | 1,291 |
| Colorectal | 528 | 997 |
| Glioma | 905 | 545 |
| Prostate | 493 | 623 |
| Pancreatic | 176 | 491 |
| Renal | 695 | 325 |
| Esophagogastric | 618 | 319 |
| Bladder | 409 | 413 |
| Melanoma | 440 | 358 |
| Thyroid | 485 | 222 |
| Ovarian | 409 | 216 |
| Endometrial | 515 | 215 |
| Head and Neck | 502 | 168 |

135 columns: `sample_id`, 96 `sbs_*` per-Mb rate columns, `n_nonsyn`,
30 `driver_*` binary flags, `n_snv_total`, `cancer_type`, `source_study`,
`cohort` (`exome`/`panel`), `tmb`, `footprint_mb`, `panel_version`
(MSK-IMPACT samples only).

## One more real-world observation worth flagging

Median TMB is *higher* in the panel cohort (4.9/Mb) than the exome cohort
(1.7/Mb), despite the panel covering far fewer genes. This is not a bug —
it's a real, expected selection effect: MSK-IMPACT testing is
disproportionately ordered for metastatic/refractory patients (who tend to
carry a higher mutation burden, partly from prior therapy exposure and
partly from enrichment for MSI/hypermutator cases relevant to
immunotherapy eligibility), while TCGA is a broader surgical cohort. This
is part of the real distribution shift the project is studying, not
something to correct away.
