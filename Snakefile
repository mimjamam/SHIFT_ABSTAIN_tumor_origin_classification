"""
SHIFT-ABSTAIN pipeline. Run with: snakemake --cores 4 --use-conda=False all

Phases map to rule groups: data acquisition -> features -> (later phases add
train/eval/figure rules as they land).
"""

TCGA_STUDIES = [
    "luad", "lusc", "brca", "coadread", "gbm", "lgg", "prad", "paad",
    "kirc", "kirp", "kich", "stad", "esca", "blca", "skcm", "thca", "ov",
    "ucec", "hnsc",
]

rule all:
    input:
        "data/processed/features.parquet",


rule fetch_reference:
    output:
        "data/reference/hg19.2bit",
    shell:
        "mkdir -p data/reference && "
        "curl -s -o {output} https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit"


rule fetch_gene_panels:
    output:
        "data/gene_panels/IMPACT341.txt",
        "data/gene_panels/IMPACT410.txt",
    script:
        "scripts/fetch_gene_panels.py"


rule fetch_tcga_maf:
    output:
        "data/tcga/{study}_tcga_pan_can_atlas_2018/data_mutations.txt",
    params:
        repo_path=lambda wc: f"public/{wc.study}_tcga_pan_can_atlas_2018/data_mutations.txt",
    shell:
        "python3 scripts/fetch_cbioportal_lfs.py {params.repo_path} {output}"


rule fetch_tcga_clinical:
    output:
        "data/tcga/{study}_tcga_pan_can_atlas_2018/data_clinical_sample.txt",
    params:
        repo_path=lambda wc: f"public/{wc.study}_tcga_pan_can_atlas_2018/data_clinical_sample.txt",
    shell:
        "python3 scripts/fetch_cbioportal_lfs.py {params.repo_path} {output}"


rule fetch_msk_impact:
    output:
        mutations="data/msk_impact_2017/data_mutations.txt",
        clinical_sample="data/msk_impact_2017/data_clinical_sample.txt",
        clinical_patient="data/msk_impact_2017/data_clinical_patient.txt",
        gene_panel_matrix="data/msk_impact_2017/data_gene_panel_matrix.txt",
    shell:
        """
        python3 scripts/fetch_cbioportal_lfs.py public/msk_impact_2017/data_mutations.txt {output.mutations}
        python3 scripts/fetch_cbioportal_lfs.py public/msk_impact_2017/data_clinical_sample.txt {output.clinical_sample}
        python3 scripts/fetch_cbioportal_lfs.py public/msk_impact_2017/data_clinical_patient.txt {output.clinical_patient}
        python3 scripts/fetch_cbioportal_lfs.py public/msk_impact_2017/data_gene_panel_matrix.txt {output.gene_panel_matrix}
        """


rule build_features:
    input:
        reference="data/reference/hg19.2bit",
        gene_panels=["data/gene_panels/IMPACT341.txt", "data/gene_panels/IMPACT410.txt"],
        tcga_mafs=expand(
            "data/tcga/{study}_tcga_pan_can_atlas_2018/data_mutations.txt", study=TCGA_STUDIES
        ),
        tcga_clinical=expand(
            "data/tcga/{study}_tcga_pan_can_atlas_2018/data_clinical_sample.txt", study=TCGA_STUDIES
        ),
        msk_mutations="data/msk_impact_2017/data_mutations.txt",
        msk_clinical="data/msk_impact_2017/data_clinical_sample.txt",
        msk_panel_matrix="data/msk_impact_2017/data_gene_panel_matrix.txt",
    output:
        "data/processed/features.parquet",
    shell:
        "python3 scripts/build_dataset.py"
