# Normalization choice (the scientific crux)

## The problem

TCGA exome samples typically carry hundreds to thousands of called
mutations; MSK-IMPACT panel samples typically carry a handful, because the
panel only sequences ~410 genes instead of the whole exome. If a classifier
is handed raw mutation counts (including raw SBS-96 counts), it can learn
to separate "exome" from "panel" — and therefore, indirectly, learn cohort
identity instead of tumor biology — just from the sheer count of mutations,
without learning anything about mutational signatures or driver biology.
That would make the shift-robustness result meaningless: the model would
"work" on panel data only by recognizing it as a different data source and
adjusting, not by generalizing.

## The fix

Every count-based feature (the 96 SBS channels, and the mutation count
behind TMB) is divided by the captured footprint of the assay, in
megabases, before it reaches the model. This turns absolute mutation counts
into mutation *rates per Mb of captured sequence* — the same units TMB is
already reported in — so a panel sample and an exome sample with the same
underlying mutation rate produce comparable feature magnitudes.

This does **not** eliminate the covariate shift (that's the point of the
project — the shift is still there and the model still needs to cope with
it), but it removes the *trivial, uninteresting* version of the shift
(depth) so what's left is the *scientifically real* version: which genes a
341-to-410-gene panel can see at all, how its variant-calling pipeline
differs from an exome pipeline (see caveat below), and how a coarser gene
list changes what mutational-signature information is even present in the
data.

## Footprint values used

| Cohort | Footprint | Source |
|---|---|---|
| TCGA (exome) | 38 Mb | Size of the Agilent SureSelect Human All Exon 38 Mb kit, the capture kit most commonly cited as the denominator in pan-cancer TCGA TMB analyses. |
| MSK-IMPACT, IMPACT341 samples | 0.901587 Mb | Reverse-engineered from MSK's own data (see below). |
| MSK-IMPACT, IMPACT410 samples | 1.021743 Mb | Reverse-engineered from MSK's own data (see below). |

**How the panel figures were derived.** We first tried to find a published
per-version footprint. The primary Zehir et al. 2017 *Nature Medicine*
paper explicitly declines to give one ("...the total genomic area where
mutations were reported, according to the version of the assay that was
used"), and the only concrete number we could find anywhere — "the probes
target approximately 1.5Mb" in the FDA De Novo classification summary for
MSK-IMPACT (DEN170058) — describes the newer 468-gene panel, not the 341-
or 410-gene versions actually used in `msk_impact_2017`.

But `data_clinical_sample.txt` in that same cohort ships a precomputed
`TMB_NONSYNONYMOUS` column — MSK's own reported clinical TMB for every
sample. Since TMB = nonsynonymous mutation count / footprint, dividing our
own nonsynonymous count for each sample by MSK's reported TMB recovers
the exact footprint MSK used to compute it:

```
implied_footprint = n_nonsyn / TMB_NONSYNONYMOUS
```

Grouped by which panel version each sample used (`data_gene_panel_matrix.txt`),
this value is a near-perfect constant:

| Panel version | n samples | Implied footprint (mode) | Samples exactly matching mode |
|---|---|---|---|
| IMPACT341 | 1,914 | 0.901587 Mb | 1,859 (97.1%) |
| IMPACT410 | 5,868 | 1.021743 Mb | 5,679 (96.8%) |

The small remainder that doesn't exactly match is most likely explained by
`TMB_NONSYNONYMOUS` being reported at limited decimal precision, or minor
differences between our nonsynonymous-classification rules and MSK's
internal ones on a handful of edge-case variants — not by a different
underlying footprint. We use these two values as per-sample footprints
(looked up by each sample's actual panel version), which is more precise
than a single blanket number and is fully reproducible from the public
data (see `scripts/build_dataset.py::build_msk`).

## A second, non-normalizable source of shift

Independent of footprint, the two pipelines report different things: TCGA
MAFs include silent, intronic, UTR, and flanking-region mutations (tens of
thousands per study); the MSK-IMPACT MAF reports almost none (1 `Silent`
row across the entire 10,945-sample cohort — clinical panels generally
don't report synonymous variants because they're not actionable). The
SBS-96 features are built from *all* reported single-base substitutions in
each file, so this pipeline-reporting difference is folded into the
mutational-context features along with the real capture-footprint
difference. We did not attempt to filter TCGA down to "what a clinical
panel pipeline would have reported," both because there's no principled
way to do that from public MAFs alone, and because this reporting gap is
itself a genuine, real-world part of what makes panel data harder to
classify from — not an artifact to be hidden. It's flagged here, and in
the README limitations section, so the reported shift penalty is
understood to include both effects rather than being misread as capturing
only the "clean" biological one.
