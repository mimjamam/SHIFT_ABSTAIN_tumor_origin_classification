# Fallback for anyone without Snakemake. `snakemake --cores 4 all` is
# preferred (see Snakefile) -- this just calls the same scripts in order.
.PHONY: all features models evaluate serve clean

all: evaluate

features: data/processed/features.parquet

data/processed/features.parquet:
	python3 scripts/fetch_gene_panels.py
	bash scripts/download_tcga_mafs.sh
	python3 scripts/build_dataset.py

models: features
	python3 scripts/train_baseline.py
	python3 scripts/train_ensemble.py
	python3 scripts/compute_thresholds.py

evaluate: models
	python3 scripts/evaluate.py

serve:
	cd api && uvicorn main:app --host 0.0.0.0 --port 8000

clean:
	rm -rf data/processed results/predictions
