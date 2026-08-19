"""Unit checks, schema, sanity constraints, and rigidity-path KGO."""
from pathlib import Path

import numpy as np
import pytest

import ParticleRigidityCalculationTools as PRCT

from atmosphericRadiationDoseAndFlux import particle as particle_mod
from atmosphericRadiationDoseAndFlux.doseAndFluxCalculator import (
    calculate_energy_integrated_flux,
    calculate_from_energy_spec,
    calculate_from_energy_spec_array,
    calculate_from_rigidity_spec_array,
)
from atmosphericRadiationDoseAndFlux.particleResponse import fullListOfDoseResponseTypes

FIXTURES = Path(__file__).resolve().parent / "fixtures"
KGO_RTOL = 1e-9
KGO_ATOL = 1e-6
EXPECTED_EXTRA_COLS = ("SEU", "SEL")


def _default_energy_bins():
    return 10 ** (0.1 * (np.array(range(1, 52)) - 1) + 1)


def _proton_rigidity_bins_from_default_energy():
    bins = _default_energy_bins()
    p = particle_mod.Particle("proton")
    return np.array(
        PRCT.convertParticleEnergyToRigidity(
            bins, particleMassAU=p.atomicMass, particleChargeAU=p.atomicCharge
        )
    )


def test_calculate_energy_integrated_flux():
    fluxes = np.array([2.0, 3.0], dtype=np.float64)
    d_e = np.array([1.0, 2.0], dtype=np.float64)
    out = calculate_energy_integrated_flux(fluxes, d_e)
    expected = np.array([2.0 * np.pi, 3.0 * 2.0 * np.pi])
    np.testing.assert_allclose(out, expected, rtol=0.0, atol=0.0)


def test_output_schema_flat_proton():
    bins = _default_energy_bins()
    alts = np.linspace(0, 15, 4)
    df = calculate_from_energy_spec_array(
        bins, np.full(50, 1.0), alts, particleName="proton"
    )
    assert len(df) == len(alts)
    assert "altitude (km)" in df.columns
    for c in fullListOfDoseResponseTypes:
        assert c in df.columns
    for c in EXPECTED_EXTRA_COLS:
        assert c in df.columns
    assert df["SEU"].iloc[0] == pytest.approx(df["tn2"].iloc[0] * 1e-13)
    assert df["SEL"].iloc[0] == pytest.approx(df["tn2"].iloc[0] * 1e-8)


def test_non_negative_outputs_physical_spectrum():
    bins = _default_energy_bins()
    alts = np.linspace(0, 60, 12)
    df = calculate_from_energy_spec_array(
        bins, np.full(50, 1.0), alts, particleName="proton"
    )
    numeric = df.drop(columns=["altitude (km)"])
    assert (numeric.to_numpy() >= 0).all()


def test_vertical_cutoff_rigidity_rejected():
    """High-level API documents unsupported cutoff; array path may still apply cutoff internally."""
    with pytest.raises(Exception, match="Vertical cutoff rigidity"):
        calculate_from_energy_spec(
            lambda x: x**-7,
            [0.0],
            particleName="proton",
            verticalCutOffRigidity=1.0,
        )


def test_invalid_particle_name():
    bins = _default_energy_bins()
    with pytest.raises(Exception, match="only protons and alpha"):
        calculate_from_energy_spec_array(
            bins, np.full(50, 1.0), [0.0], particleName="muon"
        )


def test_kgo_rigidity_flat_proton_60km():
    path = FIXTURES / "kgo_rigidity_flat_proton_60km.npz"
    if not path.is_file():
        pytest.skip(f"Missing {path}; run scripts/generate_kgo_reference.py")
    kgo = np.load(path)
    rigidity_bins = _proton_rigidity_bins_from_default_energy()
    df = calculate_from_rigidity_spec_array(
        rigidity_bins, np.full(50, 1.0), [60.0], particleName="proton"
    )
    np.testing.assert_allclose(
        df["altitude (km)"].to_numpy(), kgo["altitude_km"], rtol=0, atol=0
    )
    for key in list(fullListOfDoseResponseTypes) + list(EXPECTED_EXTRA_COLS):
        np.testing.assert_allclose(
            df[key].to_numpy(), kgo[key], rtol=KGO_RTOL, atol=KGO_ATOL
        )


def test_finite_outputs_power_law():
    bins = _default_energy_bins()
    df = calculate_from_energy_spec_array(
        bins,
        np.array([x**-7 for x in (bins[1:] + bins[:-1]) / 2]),
        [10.0, 20.0],
        particleName="proton",
    )
    num = df.select_dtypes(include=[np.floating])
    assert np.isfinite(num.to_numpy()).all()
