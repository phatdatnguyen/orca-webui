"""Tests for the spectrum-generation helpers (Plotly figures + math kernels)."""
import numpy as np
import plotly.graph_objects as go

from orca_webui import utils


# --- broadening kernels ------------------------------------------------------

def test_lorentzian_ir_peaks_at_position():
    # Value at the peak equals the intensity; falls off away from it.
    assert utils.lorentzian_ir(1000.0, 1000.0, 1.0, width=10) == 1.0
    assert utils.lorentzian_ir(1050.0, 1000.0, 1.0, width=10) < 0.1


def test_gaussian_normalized_at_center():
    assert utils.gaussian(5.0, 5.0, fwhm=2.0) == 1.0
    assert utils.gaussian(100.0, 5.0, fwhm=2.0) < 1e-6


def test_lorentzian_nmr_peaks_at_center():
    assert utils.lorentzian_nmr(2.0, 2.0, gamma=0.5) == 1.0


# --- IR spectrum -------------------------------------------------------------

def test_ir_spectrum_returns_figure():
    fig = utils.generate_ir_spectrum_interactive([1000.0, 1600.0], [50.0, 80.0])
    assert isinstance(fig, go.Figure)
    y = np.asarray(fig.data[0].y)
    assert not np.isnan(y).any()


def test_ir_spectrum_none_when_empty():
    assert utils.generate_ir_spectrum_interactive([], []) is None


# --- absorption / emission spectrum -----------------------------------------

def test_absorption_spectrum_returns_figure():
    fig = utils.generate_absorption_emission_spectrum_interactive(
        wavelengths=[250.0, 400.0], oscs=[0.5, 0.2])
    assert isinstance(fig, go.Figure)
    y = np.asarray(fig.data[0].y)
    assert not np.isnan(y).any()


def test_absorption_spectrum_none_when_empty():
    assert utils.generate_absorption_emission_spectrum_interactive(
        wavelengths=[], oscs=[]) is None


def test_absorption_spectrum_zero_oscillators_no_nan():
    """Regression: all-zero oscillator strengths must not produce a NaN spectrum."""
    fig = utils.generate_absorption_emission_spectrum_interactive(
        wavelengths=[250.0, 400.0], oscs=[0.0, 0.0])
    assert isinstance(fig, go.Figure)
    y = np.asarray(fig.data[0].y)
    assert not np.isnan(y).any()
    assert np.allclose(y, 0.0)


# --- NMR spectrum ------------------------------------------------------------

def test_nmr_spectrum_returns_figure():
    shifts_data = [(1, "H", 1.2), (2, "H", 1.2), (3, "H", 3.6)]
    fig = utils.generate_nmr_spectrum_interactive(shifts_data, element_symbol="H")
    assert isinstance(fig, go.Figure)
    y = np.asarray(fig.data[0].y)
    assert not np.isnan(y).any()
    assert y.max() > 0


def test_nmr_spectrum_none_for_missing_element():
    shifts_data = [(1, "H", 1.2)]
    assert utils.generate_nmr_spectrum_interactive(shifts_data, element_symbol="C") is None
