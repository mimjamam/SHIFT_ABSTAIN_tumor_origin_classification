# Fallback for anyone without Snakemake. `snakemake --cores 4 all` is
# preferred (see Snakefile) -- this just calls the same scripts in order.
.PHONY: all features clean

all: features

features: data/processed/features.parquet

data/processed/features.parquet:
	python3 scripts/fetch_gene_panels.py
	bash scripts/download_tcga_mafs.sh
	python3 scripts/build_dataset.py

clean:
	rm -rf data/processed
