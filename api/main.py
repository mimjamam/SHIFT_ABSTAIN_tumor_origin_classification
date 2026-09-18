"""
FastAPI backend + Gradio UI for the SHIFT-ABSTAIN demo, served as one app
on one port: `uvicorn api.main:app`.

- POST /predict: upload a MAF, get back predicted tumor type, calibrated
  confidence, and a PREDICT/ABSTAIN decision for every sample in the file.
- GET /figures/{name}: serves the Phase 4 evaluation figures.
- GET /: the Gradio UI, mounted into this same app.
"""
import os
import shutil
import tempfile

import gradio as gr
import pandas as pd
from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, RedirectResponse

from inference import ASSAY_TYPES, get_predictor

ROOT = os.path.join(os.path.dirname(__file__), "..")
FIGURES_DIR = os.path.join(ROOT, "results", "figures")

app = FastAPI(title="SHIFT-ABSTAIN", description="Tumor-origin classification that knows when to defer")


@app.get("/health")
def health():
    return {"status": "ok"}


@app.get("/")
def root():
    return RedirectResponse(url="/ui")


@app.post("/predict")
async def predict(
    maf_file: UploadFile = File(...),
    assay_type: str = Form(...),
    coverage: int = Form(90),
):
    if assay_type not in ASSAY_TYPES:
        raise HTTPException(400, f"assay_type must be one of {list(ASSAY_TYPES)}")
    if coverage not in (80, 90, 95):
        raise HTTPException(400, "coverage must be 80, 90, or 95")

    with tempfile.NamedTemporaryFile(suffix=".maf", delete=False) as tmp:
        shutil.copyfileobj(maf_file.file, tmp)
        tmp_path = tmp.name
    try:
        predictor = get_predictor()
        results = predictor.predict_maf(tmp_path, assay_type, coverage)
    finally:
        os.unlink(tmp_path)

    if not results:
        raise HTTPException(422, "no usable single-nucleotide variants found in the uploaded MAF")
    return {"n_samples": len(results), "results": results}


@app.get("/figures/{name}")
def figure(name: str):
    path = os.path.join(FIGURES_DIR, name)
    if not os.path.exists(path) or not name.endswith(".png"):
        raise HTTPException(404, "figure not found")
    return FileResponse(path)


def _gradio_predict(maf_file, assay_type_label, coverage_label):
    if maf_file is None:
        return pd.DataFrame(), "Upload a MAF file first."

    label_to_key = {v["label"]: k for k, v in ASSAY_TYPES.items()}
    assay_type = label_to_key[assay_type_label]
    coverage = int(coverage_label.rstrip("%"))

    predictor = get_predictor()
    results = predictor.predict_maf(maf_file.name, assay_type, coverage)
    if not results:
        return pd.DataFrame(), "No usable single-nucleotide variants found in this file."

    df = pd.DataFrame([
        {
            "sample_id": r["sample_id"],
            "predicted_type": r["predicted_type"],
            "confidence": round(r["confidence"], 3),
            "decision": r["decision"],
            "n_mutations": r["n_mutations_used"],
            "tmb": round(r["tmb"], 2),
        }
        for r in results
    ])

    first = results[0]
    top3_str = ", ".join(f"{t['cancer_type']} ({t['probability']:.2f})" for t in first["top3"])
    summary = (
        f"Regime: {first['regime']} (footprint {first['footprint_mb']:.3f} Mb) | "
        f"threshold @ {coverage}% coverage: {first['threshold_used']:.3f}\n"
        f"First sample ({first['sample_id']}) top-3: {top3_str}"
    )
    return df, summary


with gr.Blocks(title="SHIFT-ABSTAIN") as ui:
    gr.Markdown(
        "# SHIFT-ABSTAIN\n"
        "Tumor-origin classification trained on TCGA exome data, evaluated under "
        "real covariate shift to MSK-IMPACT targeted-panel data. Upload a MAF file, "
        "pick the assay it came from, and get a prediction with a calibrated "
        "PREDICT/ABSTAIN decision instead of a bare guess.\n\n"
        "**Read this first:** shift substantially degrades both accuracy and "
        "calibration here (68.7% -> 33.5% full-coverage accuracy; see "
        "`docs/phase4_evaluation.md`). Abstention helps but does not close that "
        "gap -- treat panel-regime predictions accordingly."
    )
    with gr.Row():
        with gr.Column():
            maf_input = gr.File(label="MAF file (data_mutations.txt format)")
            assay_dropdown = gr.Dropdown(
                choices=[v["label"] for v in ASSAY_TYPES.values()],
                value=ASSAY_TYPES["exome"]["label"],
                label="Assay type",
            )
            coverage_dropdown = gr.Dropdown(
                choices=["80%", "90%", "95%"], value="90%",
                label="Coverage operating point (higher coverage = more predictions, lower selective accuracy)",
            )
            submit_btn = gr.Button("Predict", variant="primary")
        with gr.Column():
            output_table = gr.Dataframe(label="Predictions")
            output_summary = gr.Textbox(label="Detail", lines=3)

    submit_btn.click(
        _gradio_predict,
        inputs=[maf_input, assay_dropdown, coverage_dropdown],
        outputs=[output_table, output_summary],
    )

    gr.Markdown("## Risk-coverage curves (from Phase 4 evaluation)")
    gr.Image(os.path.join(FIGURES_DIR, "risk_coverage_curves.png"), label=None, show_label=False)


app = gr.mount_gradio_app(app, ui, path="/ui")
