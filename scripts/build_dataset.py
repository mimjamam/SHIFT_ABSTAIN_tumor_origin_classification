"""
Build the final modeling table: one row per sample, with SBS-96 per-Mb
mutation rates, TMB, and driver-gene flags, labeled with cohort
(exome/panel) and a shared 14-class cancer-type taxonomy.

Normalization: SBS-96 counts and TMB are both divided by the captured
footprint in Mb. TCGA (exome) uses 38 Mb -- the size of the Agilent
SureSelect Human All Exon kit used across most TCGA studies, and the
standard denominator in pan-cancer TCGA TMB literature. MSK-IMPACT (panel)
uses per-panel-version footprints (0.901587 Mb for IMPACT341, 1.021743 Mb
for IMPACT410) reverse-engineered from MSK's own precomputed
TMB_NONSYNONYMOUS column in data_clinical_sample.txt: dividing our
nonsynonymous mutation count by their TMB recovers the exact denominator
they used, and it's a near-constant within each panel version (matches on
>96% of samples). See docs/normalization.md for the full derivation and
caveats.

Usage: python3 build_dataset.py
Writes: data/processed/features.parquet
"""
import os
import sys
import pandas as pd
import py2bit

sys.path.insert(0, os.path.dirname(__file__))
from extract_features import extract_maf_features

ROOT = os.path.join(os.path.dirname(__file__), "..")
REFERENCE_2BIT = os.path.join(ROOT, "data", "reference", "hg19.2bit")

EXOME_FOOTPRINT_MB = 38.0

# Reverse-engineered from MSK's own TMB_NONSYNONYMOUS column: see
# docs/normalization.md for the derivation (n_nonsyn / TMB_NONSYNONYMOUS,
# which is a near-constant within each panel version).
PANEL_FOOTPRINT_MB = {
    "IMPACT341": 0.901587,
    "IMPACT410": 1.021743,
}

DRIVER_GENES = [
    "TP53", "KRAS", "PIK3CA", "PTEN", "APC", "EGFR", "BRAF", "IDH1", "IDH2",
    "VHL", "RB1", "NF1", "ARID1A", "CDKN2A", "SMAD4", "STK11", "ATM",
    "BRCA1", "BRCA2", "CTNNB1", "KMT2D", "NOTCH1", "FBXW7", "MYC", "ERBB2",
    "MET", "KIT", "CDH1", "GATA3", "RUNX1",
]

# TCGA PanCancer Atlas study code -> shared class
TCGA_STUDY_TO_CLASS = {
    "luad": "NSCLC", "lusc": "NSCLC",
    "brca": "Breast",
    "coadread": "Colorectal",
    "gbm": "Glioma", "lgg": "Glioma",
    "prad": "Prostate",
    "paad": "Pancreatic",
    "kirc": "Renal", "kirp": "Renal", "kich": "Renal",
    "stad": "Esophagogastric", "esca": "Esophagogastric",
    "blca": "Bladder",
    "skcm": "Melanoma",
    "thca": "Thyroid",
    "ov": "Ovarian",
    "ucec": "Endometrial",
    "hnsc": "Head and Neck",
}

# MSK-IMPACT CANCER_TYPE -> shared class
MSK_CANCERTYPE_TO_CLASS = {
    "Non-Small Cell Lung Cancer": "NSCLC",
    "Breast Cancer": "Breast",
    "Colorectal Cancer": "Colorectal",
    "Glioma": "Glioma",
    "Prostate Cancer": "Prostate",
    "Pancreatic Cancer": "Pancreatic",
    "Renal Cell Carcinoma": "Renal",
    "Esophagogastric Cancer": "Esophagogastric",
    "Bladder Cancer": "Bladder",
    "Melanoma": "Melanoma",
    "Thyroid Cancer": "Thyroid",
    "Ovarian Cancer": "Ovarian",
    "Endometrial Cancer": "Endometrial",
    "Head and Neck Cancer": "Head and Neck",
}


def normalize_block(df: pd.DataFrame, footprint_mb) -> pd.DataFrame:
    """Divide SBS-96 counts and TMB numerator by footprint. footprint_mb may
    be a scalar (applied to every row) or a per-row Series aligned to df's
    index (e.g. panel-version-specific footprints)."""
    sbs_cols = [c for c in df.columns if c.startswith("sbs_")]
    df[sbs_cols] = df[sbs_cols].div(footprint_mb, axis=0)
    df["tmb"] = df["n_nonsyn"] / footprint_mb
    df["footprint_mb"] = footprint_mb
    return df


def build_tcga() -> pd.DataFrame:
    tb = py2bit.open(REFERENCE_2BIT)
    frames = []
    for study_code, cls in TCGA_STUDY_TO_CLASS.items():
        study_dir = f"{study_code}_tcga_pan_can_atlas_2018"
        maf_path = os.path.join(ROOT, "data", "tcga", study_dir, "data_mutations.txt")
        if not os.path.exists(maf_path):
            print(f"MISSING: {maf_path} -- skipping {study_code}")
            continue
        feats = extract_maf_features(maf_path, tb, DRIVER_GENES)
        feats["cancer_type"] = cls
        feats["source_study"] = study_dir
        feats["cohort"] = "exome"
        feats = normalize_block(feats, EXOME_FOOTPRINT_MB)
        frames.append(feats)
        print(f"{study_dir}: {len(feats)} samples -> {cls} "
              f"(context mismatches: {feats.attrs.get('n_context_mismatch', 0)}/{feats.attrs.get('n_snv_input', 0)})")
    return pd.concat(frames)


def build_msk() -> pd.DataFrame:
    tb = py2bit.open(REFERENCE_2BIT)
    maf_path = os.path.join(ROOT, "data", "msk_impact_2017", "data_mutations.txt")
    clinical_path = os.path.join(ROOT, "data", "msk_impact_2017", "data_clinical_sample.txt")
    panel_path = os.path.join(ROOT, "data", "msk_impact_2017", "data_gene_panel_matrix.txt")

    feats = extract_maf_features(maf_path, tb, DRIVER_GENES)
    print(f"msk_impact_2017: {len(feats)} samples with >=1 mutation "
          f"(context mismatches: {feats.attrs.get('n_context_mismatch', 0)}/{feats.attrs.get('n_snv_input', 0)})")

    clinical = pd.read_csv(clinical_path, sep="\t", comment="#")
    clinical = clinical.set_index("SAMPLE_ID")
    feats["cancer_type"] = clinical["CANCER_TYPE"].reindex(feats.index).map(MSK_CANCERTYPE_TO_CLASS)
    feats["source_study"] = "msk_impact_2017"
    feats["cohort"] = "panel"
    feats = feats.dropna(subset=["cancer_type"])

    panel_matrix = pd.read_csv(panel_path, sep="\t").set_index("SAMPLE_ID")
    panel_version = panel_matrix["mutations"].reindex(feats.index)
    feats["panel_version"] = panel_version
    feats = feats.dropna(subset=["panel_version"])
    footprint = feats["panel_version"].map(PANEL_FOOTPRINT_MB)

    feats = normalize_block(feats, footprint)
    return feats


if __name__ == "__main__":
    tcga = build_tcga()
    msk = build_msk()
    combined = pd.concat([tcga, msk])
    combined.index.name = "sample_id"

    out_dir = os.path.join(ROOT, "data", "processed")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "features.parquet")
    combined.reset_index().to_parquet(out_path, index=False)

    print("\n=== final dataset ===")
    print(f"total samples: {len(combined)}")
    print(combined.groupby(["cohort", "cancer_type"]).size().unstack(fill_value=0).T)
    print(f"\nwrote {out_path}")
