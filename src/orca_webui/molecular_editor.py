"""Repository-owned molecular sketcher using Gradio HTML and local Ketcher assets."""
from pathlib import Path
from uuid import uuid4

import gradio as gr
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from rdkit import Chem


APP_ROOT = Path(__file__).resolve().parent
KETCHER_DIR = APP_ROOT / "assets" / "ketcher"

# The template does not contain the molecule value. Updating SMILES must not replace
# the iframe, lose undo history, or interpret chemical/user input as HTML.
EDITOR_TEMPLATE = """
<div class="molecule-editor">
  <iframe title="Molecule drawing editor" class="molecule-canvas"></iframe>
  <div class="molecule-actions">
    <button type="button" class="use-drawing" disabled>Use drawing</button>
    <span class="editor-status" role="status" aria-live="polite">Loading molecular editor…</span>
  </div>
</div>
"""

EDITOR_CSS = """
.molecule-canvas { width: 100%; height: 580px; border: 1px solid var(--border-color-primary); border-radius: 6px; background: white; }
.molecule-actions { display: flex; align-items: center; flex-wrap: wrap; gap: 12px; padding-top: 10px; }
.use-drawing { padding: 8px 16px; border: 1px solid var(--button-primary-border-color); border-radius: 6px; background: var(--button-primary-background-fill); color: var(--button-primary-text-color); cursor: pointer; font: inherit; }
.use-drawing:disabled { opacity: 0.5; cursor: default; }
.editor-status { font-size: 0.9em; color: var(--body-text-color-subdued); }
.editor-status[data-error="true"] { color: var(--error-text-color, #b42318); }
"""


def mount_molecular_editor(app):
    """Mount permanent assets separately from the generated /static viewers."""
    if (KETCHER_DIR / "index.html").is_file():
        app.mount("/molecule-editor", StaticFiles(directory=KETCHER_DIR, html=True), name="molecule-editor")
    else:
        @app.get("/molecule-editor/{asset_path:path}", include_in_schema=False)
        def editor_not_installed(asset_path: str):
            return HTMLResponse(
                "<p>The drawing editor is not installed. You can still enter SMILES directly.</p>",
                status_code=503,
            )
    return app


def molecular_editor():
    """Build the drawing surface using Gradio's native HTML component API."""
    return gr.HTML(
        value={"smiles": "", "load_id": None},
        label="Molecule",
        html_template=EDITOR_TEMPLATE,
        css_template=EDITOR_CSS,
        js_on_load=(APP_ROOT / "frontend" / "molecular_editor.js").read_text(encoding="utf-8"),
        editor_available=(KETCHER_DIR / "index.html").is_file(),
        elem_id="molecular-editor",
    )


def _canonical_smiles(smiles):
    """Validate text from either side of the editor without interpreting HTML."""
    if not isinstance(smiles, str):
        raise ValueError("SMILES must be text.")
    smiles = smiles.strip()
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or not mol.GetNumAtoms():
            raise ValueError("Invalid SMILES.")
        smiles = Chem.MolToSmiles(mol, canonical=True)
    return smiles


def load_smiles_into_editor(smiles):
    """Validate before changing the drawing; each explicit load gets a new ID."""
    try:
        smiles = _canonical_smiles("" if smiles is None else smiles)
    except (TypeError, ValueError, RuntimeError):
        gr.Warning("The SMILES string is invalid. The drawing has not been changed.")
        return gr.skip()
    return {"smiles": smiles, "load_id": uuid4().hex}


def smiles_from_editor(value):
    """Accept only molecule text from the browser, then validate it with RDKit."""
    if not isinstance(value, dict) or not isinstance(value.get("smiles"), str):
        gr.Warning("The drawing could not be read. Please try Use drawing again.")
        return gr.skip()
    try:
        return _canonical_smiles(value["smiles"])
    except (TypeError, ValueError, RuntimeError):
        gr.Warning("The drawing does not describe a valid molecule. Please check its atoms and bonds.")
        return gr.skip()
