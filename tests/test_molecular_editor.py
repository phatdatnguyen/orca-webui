"""Contracts between the local sketcher, the SMILES field, and conformer generation."""
from pathlib import Path

import pytest
from rdkit import Chem

pytest.importorskip("gradio")

from fastapi import FastAPI
from orca_webui import conformer_generation
from orca_webui import molecular_editor as editor


@pytest.fixture(autouse=True)
def collect_warnings(monkeypatch):
    messages = []
    monkeypatch.setattr(editor.gr, "Warning", messages.append)
    return messages


def test_explicit_smiles_load_canonicalizes_and_can_repeat():
    first = editor.load_smiles_into_editor(" OCC ")
    second = editor.load_smiles_into_editor(" OCC ")
    assert first["smiles"] == second["smiles"] == "CCO"
    assert first["load_id"] and first["load_id"] != second["load_id"]


@pytest.mark.parametrize("smiles", [None, "", "   "])
def test_blank_smiles_load_is_an_explicit_clear(smiles):
    command = editor.load_smiles_into_editor(smiles)
    assert command["smiles"] == ""
    assert command["load_id"]


@pytest.mark.parametrize("smiles", ["not a molecule", "<script>alert(1)</script>", 123, [], {}])
def test_invalid_smiles_load_keeps_current_drawing(smiles, collect_warnings):
    assert editor.load_smiles_into_editor(smiles) == editor.gr.skip()
    assert collect_warnings and "has not been changed" in collect_warnings[0]


@pytest.mark.parametrize("smiles", ["[NH4+]", "N[C@@H](C)C(=O)O", "c1ccccc1"])
def test_editor_roundtrip_preserves_charge_stereochemistry_and_aromaticity(smiles):
    command = editor.load_smiles_into_editor(smiles)
    result = editor.smiles_from_editor({"smiles": command["smiles"], "load_id": None})
    assert result == Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True)


@pytest.mark.parametrize("payload", [None, "CCO", [], {}, {"smiles": None}, {"smiles": 123},
                                      {"smiles": "not a molecule"}])
def test_invalid_editor_payload_keeps_current_smiles(payload, collect_warnings):
    assert editor.smiles_from_editor(payload) == editor.gr.skip()
    assert collect_warnings


def test_cleared_drawing_clears_smiles():
    assert editor.smiles_from_editor({"smiles": "  ", "load_id": None}) == ""


def test_editor_events_use_matching_component_values():
    with editor.gr.Blocks() as blocks:
        directory = editor.gr.State()
        files = editor.gr.State()
        status = editor.gr.Markdown()
        conformer_generation.conformer_generation_tab_content(directory, files, status)
    handlers = {event.fn: event for event in blocks.fns.values()}
    from_drawing = handlers[editor.smiles_from_editor]
    to_drawing = handlers[editor.load_smiles_into_editor]
    generate = handlers[conformer_generation.on_generate_conformers]
    surface = from_drawing.inputs[0]
    smiles_field = from_drawing.outputs[0]

    assert to_drawing.inputs == [smiles_field]
    assert to_drawing.outputs == [surface]
    assert generate.inputs[1] is smiles_field
    assert surface.preprocess({"smiles": "O", "load_id": None}) == {"smiles": "O", "load_id": None}
    assert surface.postprocess(editor.load_smiles_into_editor("O"))["smiles"] == "O"
    # SMILES travels as structured data; including it in the HTML template would
    # rebuild the iframe on every load and could interpret molecule text as markup.
    assert "${value}" not in surface.html_template
    assert "{{value}}" not in surface.html_template
    assert surface.props["editor_available"] is True


def test_editor_smiles_can_generate_a_real_conformer(tmp_path):
    class Progress:
        def __call__(self, fraction, desc=None):
            pass

        def tqdm(self, iterable, total=None, desc=None):
            return iterable

    smiles = editor.smiles_from_editor({"smiles": "O", "load_id": None})
    status, files, table = conformer_generation.on_generate_conformers(
        str(tmp_path), smiles, 0, 1, 1, 0.1, 0.5, "water", "xyz", progress=Progress())
    assert "color:green" in status, status
    assert files == ["water_1.xyz"]
    assert len(table) == 1
    assert "O " in (tmp_path / "water_1.xyz").read_text()


def test_editor_assets_are_served_from_the_package(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    app = editor.mount_molecular_editor(FastAPI())
    route = next(route for route in app.routes if route.name == "molecule-editor")
    filename, details = route.app.lookup_path("index.html")
    assert Path(filename) == editor.KETCHER_DIR / "index.html"
    assert details.st_size > 0
    assert not (tmp_path / "assets").exists()
    # StaticFiles must keep requests inside the package's editor directory.
    assert route.app.lookup_path("../../molecular_editor.py")[1] is None


def test_incomplete_source_checkout_keeps_smiles_entry_usable(tmp_path, monkeypatch):
    monkeypatch.setattr(editor, "KETCHER_DIR", tmp_path / "missing")
    app = editor.mount_molecular_editor(FastAPI())
    fallback = next(route for route in app.routes if route.path == "/molecule-editor/{asset_path:path}")
    response = fallback.endpoint("index.html")
    assert response.status_code == 503
    assert b"SMILES" in response.body
    assert editor.molecular_editor().props["editor_available"] is False
