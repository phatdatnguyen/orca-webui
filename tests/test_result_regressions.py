"""Result handling with incomplete cclib records and user-provided export names."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("gradio")
pytest.importorskip("cclib")

from cclib.parser.utils import convertor
from orca_webui import result


def _load(monkeypatch, tmp_path, **attributes):
    data = SimpleNamespace(**attributes)
    (tmp_path / "sample.log").write_text("Synthetic calculation output\n")
    monkeypatch.setattr(result.cclib.io, "ccopen", lambda path: SimpleNamespace(parse=lambda: data))
    outputs = result.on_load_result_file(str(tmp_path), "sample.log")
    assert len(outputs) == 25
    assert "color:green" in outputs[0], outputs[0]
    return outputs


def test_scf_energy_does_not_require_molecular_orbitals(monkeypatch, tmp_path):
    outputs = _load(monkeypatch, tmp_path, scfenergies=[convertor(-75.0, "hartree", "eV")])
    assert outputs[2]["visible"] is True
    assert outputs[3] == "-75.0000 (hartree)"
    assert outputs[5] is None


def test_orbital_units_and_unrestricted_spin_channels(monkeypatch, tmp_path):
    orbital_energies = convertor(np.array([[-1.0, -0.5], [-0.9, -0.4]]), "hartree", "eV")
    outputs = _load(monkeypatch, tmp_path, moenergies=orbital_energies, homos=[1, 0])
    table = outputs[5]
    assert list(table["Energy (hartree)"]) == ["-1.0000", "-0.5000", "-0.9000", "-0.4000"]
    assert list(table["Molecular orbital"]) == [
        "MO 1 (alpha)", "MO 2 (alpha) (HOMO)", "MO 1 (beta) (HOMO)", "MO 2 (beta)"]


def test_orbital_table_without_homo_metadata(monkeypatch, tmp_path):
    outputs = _load(monkeypatch, tmp_path, moenergies=[[-27.2114]])
    assert list(outputs[5]["Molecular orbital"]) == ["MO 1"]


@pytest.mark.parametrize("with_ir", [False, True])
def test_partial_frequency_output_preserves_available_results(monkeypatch, tmp_path, with_ir):
    attributes = {"scfenergies": [-100.0], "vibfreqs": [1000.0, 1500.0], "zpve": 0.0123}
    if with_ir:
        attributes["vibirs"] = [5.0, 10.0]
    outputs = _load(monkeypatch, tmp_path, **attributes)
    assert outputs[2]["visible"] is True
    assert outputs[8]["visible"] is True
    assert list(outputs[9]["Frequencies"]) == [1000.0, 1500.0]
    assert (outputs[10] is not None) == with_ir
    assert outputs[11].iloc[0].tolist() == ["0.0123", "", "", ""]


def test_empty_optional_results_do_not_crash_loader(monkeypatch, tmp_path):
    outputs = _load(monkeypatch, tmp_path, scfenergies=[], moenergies=[], vibfreqs=[])
    assert outputs[2]["visible"] is False
    assert outputs[8]["visible"] is False


def test_optimization_plot_does_not_retain_pyplot_figures(monkeypatch, tmp_path):
    import matplotlib.pyplot as plt

    previous_figures = plt.get_fignums()
    outputs = _load(monkeypatch, tmp_path, scfenergies=[-100.0, -102.0])
    assert outputs[6]["visible"] is True
    assert plt.get_fignums() == previous_figures
    assert outputs[7].axes[0].get_xlabel() == "Optimization step"


def test_unknown_log_reports_recognition_error_and_clears_outputs(monkeypatch, tmp_path):
    (tmp_path / "unknown.log").write_text("Not an ORCA output")
    monkeypatch.setattr(result.cclib.io, "ccopen", lambda path: None)
    outputs = result.on_load_result_file(str(tmp_path), "unknown.log")
    assert len(outputs) == 25
    assert "not a recognized" in outputs[0]
    assert "Upgrade cclib" not in outputs[0]
    assert outputs[1] is None
    for i in [2, 6, 8, 13, 18]:
        assert outputs[i]["visible"] is False


@pytest.mark.parametrize("handler", [result.on_show_absorption_spectrum,
                                     result.on_show_emission_spectrum,
                                     result.on_show_ecd_spectrum])
def test_spectra_do_not_require_optional_configuration_or_symmetry(handler):
    data = SimpleNamespace(etenergies=[20000.0, 25000.0], etoscs=[0.3, 0.2], etrotats=[-1.0, 2.0])
    table, figure = handler(data)
    assert figure is not None
    assert set(table["Transitions"]) == {""}
    assert set(table["Symmetry"]) == {""}


def test_ecd_does_not_require_oscillator_strengths():
    table, figure = result.on_show_ecd_spectrum(SimpleNamespace(etenergies=[20000.0], etrotats=[-1.0]))
    assert table["Oscillator Strength"].isna().all()
    assert figure is not None


def test_transition_metadata_remains_aligned_when_partial():
    data = SimpleNamespace(etenergies=[20000.0, 25000.0], etoscs=[0.3, 0.2],
                           etsecs=[[((4, 0), (5, 0), 0.7)]], etsyms=["Singlet"])
    table, _ = result.on_show_absorption_spectrum(data)
    assert list(table["Symmetry"]) == ["", "Singlet"]
    assert table.iloc[1]["Transitions"] == "5 -> 6 (coefficient: 0.7000)"


@pytest.mark.parametrize("handler", [result.on_show_absorption_spectrum,
                                     result.on_show_emission_spectrum,
                                     result.on_show_ecd_spectrum])
@pytest.mark.parametrize("energies,strengths", [([], []), ([0.0], [1.0]),
                                               ([20000.0, 25000.0], [1.0])])
def test_invalid_or_empty_transition_data_does_not_raise(monkeypatch, handler, energies, strengths):
    messages = []
    monkeypatch.setattr(result.gr, "Warning", messages.append)
    data = SimpleNamespace(etenergies=energies, etoscs=strengths, etrotats=strengths)
    assert handler(data) == (None, None)
    if energies:
        assert messages


@pytest.mark.parametrize("value", ["nan", "inf", "-inf", "", None, "invalid"])
def test_nmr_reference_falls_back_for_nonfinite_values(value):
    assert result._parse_reference(value, 32.9) == 32.9


def test_export_preserves_working_directory_and_refreshes_files(tmp_path):
    status, files = result.on_export_data(str(tmp_path), "peaks", pd.DataFrame({"Intensity": [1.0]}))
    assert "color:green" in status
    assert files == ["peaks.csv"]
    assert pd.read_csv(tmp_path / "peaks.csv")["Intensity"].tolist() == [1.0]


@pytest.mark.parametrize("name", ["../outside", "/tmp/outside", "", "nested/file"])
def test_export_rejects_invalid_names(tmp_path, name):
    status, _ = result.on_export_data(str(tmp_path), name, pd.DataFrame({"Intensity": [1.0]}))
    assert "color:red" in status
    assert list(tmp_path.iterdir()) == []


def test_export_rejects_missing_data(tmp_path):
    status, files = result.on_export_data(str(tmp_path), "peaks", None)
    assert "no result data" in status
    assert files == []


def test_file_refresh_preserves_selected_result():
    update = result.on_working_directory_file_list_change(["a.log", "z.log", "peaks.csv"], "z.log")
    assert update["value"] == "z.log"


def test_changing_directory_clears_loaded_result_state_and_exports():
    outputs = result.on_working_directory_change()
    assert len(outputs) == 25
    assert outputs[0] == ""
    assert outputs[1] is None
    for i in [2, 6, 8, 13, 18]:
        assert outputs[i]["visible"] is False
    for i in [12, 16, 24]:
        assert outputs[i]["interactive"] is False
    for i in [5, 9, 14, 19]:
        assert outputs[i] is None


def test_directory_clear_and_loader_share_gradio_output_wiring():
    with result.gr.Blocks() as blocks:
        directory = result.gr.State()
        files = result.gr.State()
        status = result.gr.Markdown()
        result.result_tab_content(directory, files, status)
    loader_events = [event for event in blocks.fns.values() if event.fn is result.on_load_result_file]
    clear_events = [event for event in blocks.fns.values() if event.fn is result.on_working_directory_change]
    assert len(clear_events) == 1
    assert len(loader_events) == 2
    assert clear_events[0].outputs == loader_events[0].outputs == loader_events[1].outputs
