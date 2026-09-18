"""
Fetch MSK-IMPACT gene panel definitions from the public cBioPortal REST API
and save each as a plain list of Hugo gene symbols, one per line.

Usage: python3 fetch_gene_panels.py
"""
import json
import subprocess
import os

PANELS = ["IMPACT341", "IMPACT410"]
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "gene_panels")

if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    for panel in PANELS:
        url = f"https://www.cbioportal.org/api/gene-panels/{panel}"
        out = subprocess.run(["curl", "-s", "--max-time", "30", url], capture_output=True, text=True)
        data = json.loads(out.stdout)
        genes = sorted(g["hugoGeneSymbol"] for g in data["genes"])
        dest = os.path.join(OUT_DIR, f"{panel}.txt")
        with open(dest, "w") as f:
            f.write("\n".join(genes) + "\n")
        print(f"{panel}: {len(genes)} genes -> {dest}")
