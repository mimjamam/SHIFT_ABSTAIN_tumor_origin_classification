"""
Extract per-sample features from one MAF file: SBS-96 trinucleotide context
raw counts, nonsynonymous mutation count (for TMB), and driver-gene flags.

This produces RAW counts only (no per-Mb normalization) — normalization by
captured footprint happens later in build_dataset.py, once, in one place,
so the choice is easy to find and change. See docs/normalization.md.

Usage (library): extract_maf_features(maf_path, tb, driver_genes) -> DataFrame
indexed by Tumor_Sample_Barcode.
"""
import sys
import os
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from sbs96 import CATEGORIES, sbs96_category

MAF_COLS = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Variant_Classification",
    "Variant_Type",
    "Tumor_Sample_Barcode",
]

# Protein-coding-altering classes -- the standard "nonsynonymous" set used
# for TMB (matches MSK-IMPACT's own TMB_NONSYNONYMOUS convention, and the
# Chalmers et al. 2017 definition commonly used in the TMB literature).
NONSYN_CLASSES = {
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Nonstop_Mutation",
    "Splice_Site",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Translation_Start_Site",
}

VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}


def _load_maf(maf_path: str) -> pd.DataFrame:
    df = pd.read_csv(
        maf_path, sep="\t", comment="#", usecols=MAF_COLS,
        dtype={"Chromosome": str, "Start_Position": "int64"},
        low_memory=False,
    )
    return df


def extract_maf_features(maf_path: str, tb, driver_genes: list[str]) -> pd.DataFrame:
    df = _load_maf(maf_path)
    df = df[df["Chromosome"].isin(VALID_CHROMS)]

    # --- SBS-96 raw counts, from all clean single-base substitutions ---
    snv = df[
        (df["Variant_Type"] == "SNP")
        & (df["Reference_Allele"].str.len() == 1)
        & (df["Tumor_Seq_Allele2"].str.len() == 1)
        & (df["Reference_Allele"].isin(list("ACGT")))
        & (df["Tumor_Seq_Allele2"].isin(list("ACGT")))
    ].copy()

    n_context_mismatch = 0
    contexts = [
        tb.sequence(f"chr{chrom}", int(pos) - 2, int(pos) + 1).upper()
        for chrom, pos in zip(snv["Chromosome"], snv["Start_Position"])
    ]
    snv["_context3"] = contexts

    categories = []
    for context3, ref, alt in zip(snv["_context3"], snv["Reference_Allele"], snv["Tumor_Seq_Allele2"]):
        cat = sbs96_category(context3, ref, alt)
        if cat is None:
            n_context_mismatch += 1
        categories.append(cat)
    snv["_sbs_category"] = categories
    snv_clean = snv.dropna(subset=["_sbs_category"])

    sbs_counts = (
        snv_clean.groupby(["Tumor_Sample_Barcode", "_sbs_category"]).size().unstack(fill_value=0)
    )
    for cat in CATEGORIES:
        if cat not in sbs_counts.columns:
            sbs_counts[cat] = 0
    sbs_counts = sbs_counts[CATEGORIES]
    sbs_counts.columns = [f"sbs_{c}" for c in CATEGORIES]

    # --- nonsynonymous count per sample (for TMB) ---
    nonsyn = df[df["Variant_Classification"].isin(NONSYN_CLASSES)]
    n_nonsyn = nonsyn.groupby("Tumor_Sample_Barcode").size().rename("n_nonsyn")

    # --- driver-gene flags: any nonsynonymous mutation in that gene ---
    driver_hits = nonsyn[nonsyn["Hugo_Symbol"].isin(driver_genes)]
    driver_flags = (
        driver_hits.groupby(["Tumor_Sample_Barcode", "Hugo_Symbol"]).size().unstack(fill_value=0)
    )
    for gene in driver_genes:
        if gene not in driver_flags.columns:
            driver_flags[gene] = 0
    driver_flags = (driver_flags[driver_genes] > 0).astype(int)
    driver_flags.columns = [f"driver_{g}" for g in driver_genes]

    # --- assemble: every sample that has at least one row in the MAF ---
    all_samples = pd.Index(df["Tumor_Sample_Barcode"].unique(), name="Tumor_Sample_Barcode")
    out = pd.DataFrame(index=all_samples)
    out = out.join(sbs_counts).fillna(0)
    out = out.join(n_nonsyn).fillna(0)
    out = out.join(driver_flags).fillna(0)

    sbs_cols = [f"sbs_{c}" for c in CATEGORIES]
    out[sbs_cols] = out[sbs_cols].astype(int)
    out["n_nonsyn"] = out["n_nonsyn"].astype(int)
    driver_cols = [f"driver_{g}" for g in driver_genes]
    out[driver_cols] = out[driver_cols].astype(int)

    out["n_snv_total"] = snv_clean.groupby("Tumor_Sample_Barcode").size().reindex(out.index).fillna(0).astype(int)
    out = out.copy()  # defragment after repeated join/assign above
    out.attrs["n_context_mismatch"] = n_context_mismatch
    out.attrs["n_snv_input"] = len(snv)
    return out
