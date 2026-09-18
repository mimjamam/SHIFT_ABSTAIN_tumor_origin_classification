# Phase 0 — data spike findings

Date: 2026-09-18

## Source

Both cohorts come from the public [cBioPortal datahub](https://github.com/cBioPortal/datahub)
(`github.com/cBioPortal/datahub`). Data files are stored as Git LFS objects
behind a custom LFS server declared in the repo's `.lfsconfig`
(`https://nsssw8k94d.execute-api.us-east-1.amazonaws.com/`). `git-lfs` isn't
installed in this environment and installing it needs root (no passwordless
sudo available), so `scripts/fetch_cbioportal_lfs.py` talks to the LFS batch
API directly: it reads the pointer file's oid/size, requests a presigned S3
URL, and downloads that. Verified working for both cohorts below.

## In-distribution: TCGA PanCancer Atlas (exome)

32 per-cancer-type study folders (`public/<code>_tcga_pan_can_atlas_2018/`).
Downloaded `data_clinical_sample.txt` for all 32 to get exact sample counts,
and pulled one full `data_mutations.txt` (GBM) to confirm MAF columns.

- Genome build: **GRCh37** (confirmed from `NCBI_Build` column).
- MAF columns present: `Hugo_Symbol`, `Chromosome`, `Start_Position`,
  `Reference_Allele`, `Tumor_Seq_Allele2`, `Variant_Classification`,
  `Variant_Type`, `Tumor_Sample_Barcode`, plus (bonus) a VEP `CONTEXT` column
  giving an 11bp sequence context around each variant — could save us from
  needing a reference FASTA for this side, though for consistency with the
  panel side (which has no `CONTEXT` column) we may still want to derive
  trinucleotide context the same way for both cohorts.
- Total samples across all 32 studies: **10,967**.

Per-study sample counts (rows in `data_clinical_sample.txt`):

| Study | N | Study | N | Study | N | Study | N |
|---|---|---|---|---|---|---|---|
| acc | 92 | esca | 182 | lusc | 487 | sarc | 255 |
| blca | 411 | gbm | 592 | meso | 87 | skcm | 448 |
| brca | 1084 | hnsc | 523 | ov | 585 | stad | 440 |
| cesc | 297 | kich | 65 | paad | 184 | tgct | 149 |
| chol | 36 | kirc | 512 | pcpg | 178 | thca | 500 |
| coadread | 594 | kirp | 283 | prad | 494 | thym | 123 |
| dlbc | 48 | laml | 200 | | | ucec | 529 |
| | | lgg | 514 | | | ucs | 57 |
| | | luad | 566 | | | uvm | 80 |

## Shifted: MSK-IMPACT targeted panel

Study `public/msk_impact_2017/` (Zehir et al. 2017 cohort).

- Genome build: **GRCh37** (same as TCGA — no liftover needed).
- MAF columns present: same core set as TCGA (`Hugo_Symbol`, `Chromosome`,
  `Start_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`,
  `Variant_Classification`, `Variant_Type`, `Tumor_Sample_Barcode`) — **no**
  `CONTEXT` column, so trinucleotide context for this side needs a reference
  genome lookup.
- `data_clinical_sample.txt` has a `CANCER_TYPE` column (58 distinct values)
  and a `TMB_NONSYNONYMOUS` column (pre-computed TMB, useful as a sanity
  check against whatever we compute ourselves).
- `data_gene_panel_matrix.txt` records which panel version each sample used:
  **IMPACT341** (n=2,809) or **IMPACT410** (n=8,136). Panel gene lists (not
  footprint sizes) are retrievable from the live cBioPortal REST API
  (`https://www.cbioportal.org/api/gene-panels/IMPACT341`, etc.) — confirmed
  341 and 410 genes respectively.
- Total samples: **10,945**.

Top cancer types by sample count:

| Cancer type | N | Cancer type | N |
|---|---|---|---|
| Non-Small Cell Lung Cancer | 1,668 | Esophagogastric Cancer | 341 |
| Breast Cancer | 1,337 | Germ Cell Tumor | 288 |
| Colorectal Cancer | 1,007 | Thyroid Cancer | 231 |
| Prostate Cancer | 717 | Ovarian Cancer | 224 |
| Glioma | 553 | Endometrial Cancer | 218 |
| Pancreatic Cancer | 502 | Head and Neck Cancer | 186 |
| Soft Tissue Sarcoma | 443 | Cancer of Unknown Primary | 186 |
| Bladder Cancer | 423 | ...(43 more, smaller) | |
| Melanoma | 365 | | |
| Renal Cell Carcinoma | 361 | | |
| Hepatobiliary Cancer | 355 | | |

## Proposed shared class set

TCGA study codes and MSK `CANCER_TYPE` don't share a vocabulary and don't
split cancers identically (e.g. MSK's "Hepatobiliary Cancer" spans what TCGA
keeps separate as LIHC + CHOL). Mapping to the coarsest common categories,
these 14 classes have reasonable sample counts (>150) on both sides:

| Class | TCGA source | TCGA N | MSK N |
|---|---|---|---|
| NSCLC | luad + lusc | 1,053 | 1,668 |
| Breast | brca | 1,084 | 1,337 |
| Colorectal | coadread | 594 | 1,007 |
| Glioma | gbm + lgg | 1,106 | 553 |
| Prostate | prad | 494 | 717 |
| Pancreatic | paad | 184 | 502 |
| Renal | kirc + kirp + kich | 860 | 361 |
| Esophagogastric | stad + esca | 622 | 341 |
| Bladder | blca | 411 | 423 |
| Melanoma | skcm | 448 | 365 |
| Thyroid | thca | 500 | 231 |
| Ovarian | ov | 585 | 224 |
| Endometrial | ucec | 529 | 218 |
| Head and Neck | hnsc | 523 | 186 |

This is a proposal, not yet implemented — needs sign-off before Phase 1.

## Open items for Phase 1 (not yet resolved)

1. **Panel footprint size (Mb)**, needed to normalize mutation counts per-Mb
   so panel samples aren't trivially separable by depth. Not present in any
   file pulled so far — MSK papers report approximate captured footprint
   sizes for IMPACT341/410/468, but I have not yet located an authoritative,
   citable source (e.g. an actual BED/interval file) for the exact figure.
   Need to find one, or compute it from panel gene exon coordinates.
2. **Trinucleotide context for MSK-IMPACT** needs a GRCh37 reference genome
   (FASTA or 2bit) lookup — this is a heavier dependency (hundreds of MB to a
   few GB depending on format), flagging per project rule to ask before
   installing/downloading anything heavy.
3. `git-lfs` is not installed and there's no passwordless sudo in this
   environment — `scripts/fetch_cbioportal_lfs.py` works around this, but if
   the user wants proper `git lfs` support later it needs `sudo apt install
   git-lfs`.
