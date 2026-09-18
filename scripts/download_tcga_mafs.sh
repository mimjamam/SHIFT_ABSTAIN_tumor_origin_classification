#!/usr/bin/env bash
# Download data_mutations.txt for the TCGA PanCancer Atlas studies used in
# the 14-class shared taxonomy (see docs/phase0_data_spike.md).
set -euo pipefail
cd "$(dirname "$0")/.."

STUDIES="luad lusc brca coadread gbm lgg prad paad kirc kirp kich stad esca blca skcm thca ov ucec hnsc"

for s in $STUDIES; do
  dir="${s}_tcga_pan_can_atlas_2018"
  dest="data/tcga/${dir}/data_mutations.txt"
  if [ -f "$dest" ]; then
    echo "already have $dest, skipping"
    continue
  fi
  python3 scripts/fetch_cbioportal_lfs.py "public/${dir}/data_mutations.txt" "$dest"
done
