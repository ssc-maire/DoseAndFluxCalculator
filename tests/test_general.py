"""Regression (KGO) and additivity tests for dose/flux calculator."""
from pathlib import Path

import numpy as np
import pytest

import ParticleRigidityCalculationTools as PRCT

import atmosphericRadiationDoseAndFlux.doseAndFluxCalculator as DAFcalc
from atmosphericRadiationDoseAndFlux import particle as particle_mod
from atmosphericRadiationDoseAndFlux.doseAndFluxCalculator import (
    calculate_from_energy_spec_array,
    calculate_from_rigidity_spec_array,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures"

# KGO files are generated on a reference environment; allow tiny FP slack for other platforms.
KGO_RTOL = 1e-9
KGO_ATOL = 1e-6

DOSE_TYPES = ["adose", "edose", "dosee", "tn1", "tn2", "tn3", "SEU", "SEL"]
ADDIVITY_RTOL = 1e-8
ADDIVITY_ATOL = 1e-5


def _default_bins_and_altitudes():
    input_energy_bins = 10 ** (0.1 * (np.array(range(1, 52)) - 1) + 1)
    input_altitudes = np.linspace(0, 60, 60)
    return input_energy_bins, input_altitudes


def _load_kgo(stem: str) -> np.lib.npyio.NpzFile:
    path = FIXTURES / f"{stem}.npz"
    if not path.is_file():
        pytest.skip(f"Missing KGO fixture: {path} (run scripts/generate_kgo_reference.py)")
    return np.load(path)


def _assert_df_matches_kgo(df, kgo: np.lib.npyio.NpzFile) -> None:
    np.testing.assert_allclose(
        df["altitude (km)"].to_numpy(),
        kgo["altitude_km"],
        rtol=0.0,
        atol=0.0,
        err_msg="altitude grid mismatch",
    )
    for key in DOSE_TYPES:
        np.testing.assert_allclose(
            df[key].to_numpy(),
            kgo[key],
            rtol=KGO_RTOL,
            atol=KGO_ATOL,
            err_msg=f"KGO mismatch for column {key}",
        )


def run_over_single_step_spectrum_only(particle_species: str):
    input_energy_bins, input_altitudes = _default_bins_and_altitudes()
    input_fluxes = np.append(1.0, np.full(49, 0.0))
    return calculate_from_energy_spec_array(
        input_energy_bins, input_fluxes, input_altitudes, particleName=particle_species
    )


def run_over_flat_spectrum(particle_species: str):
    input_energy_bins, input_altitudes = _default_bins_and_altitudes()
    input_fluxes = np.full(50, 1.0)
    return calculate_from_energy_spec_array(
        input_energy_bins, input_fluxes, input_altitudes, particleName=particle_species
    )


def _energy_bin_midpoints(input_energy_bins: np.ndarray) -> np.ndarray:
    return (input_energy_bins[1:] + input_energy_bins[:-1]) / 2.0


def run_powerlaw_em7_spectrum(particle_species: str):
    bins, alts = _default_bins_and_altitudes()
    fluxes = _energy_bin_midpoints(bins) ** (-7.0)
    return calculate_from_energy_spec_array(bins, fluxes, alts, particleName=particle_species)


def run_powerlaw_em2_proton():
    bins, alts = _default_bins_and_altitudes()
    fluxes = _energy_bin_midpoints(bins) ** (-2.0)
    return calculate_from_energy_spec_array(bins, fluxes, alts, particleName="proton")


def _rigidity_bins_for_species(species: str, energy_bins: np.ndarray) -> np.ndarray:
    p = particle_mod.Particle(species)
    return np.array(
        PRCT.convertParticleEnergyToRigidity(
            energy_bins, particleMassAU=p.atomicMass, particleChargeAU=p.atomicCharge
        )
    )


def run_rigidity_powerlaw_rm7(particle_species: str):
    bins, alts = _default_bins_and_altitudes()
    rb = _rigidity_bins_for_species(particle_species, bins)
    rmid = (rb[1:] + rb[:-1]) / 2.0
    fluxes = rmid ** (-7.0)
    return calculate_from_rigidity_spec_array(rb, fluxes, alts, particleName=particle_species)


def _assert_additivity(dfp, dfa, dfboth) -> None:
    for dose_type in DOSE_TYPES:
        np.testing.assert_allclose(
            dfp[dose_type].to_numpy() + dfa[dose_type].to_numpy(),
            dfboth[dose_type].to_numpy(),
            rtol=ADDIVITY_RTOL,
            atol=ADDIVITY_ATOL,
            err_msg=f"proton+alpha != both for {dose_type}",
        )


def test_function_input_rigidity():
    output_df_proton = DAFcalc.calculate_from_rigidity_spec(
        lambda x: x**-7, [60.0], particleName="proton"
    )
    output_df_alpha = DAFcalc.calculate_from_rigidity_spec(
        lambda x: x**-7, [60.0], particleName="alpha"
    )
    output_df_both = DAFcalc.calculate_from_rigidity_spec(
        lambda x: x**-7, [60.0], particleName="both"
    )
    _assert_additivity(output_df_proton, output_df_alpha, output_df_both)


def test_function_input_rigidity2():
    input_rigidity_fn = lambda x: 27593.36 * (x**-2.82844)
    output_df_proton = DAFcalc.calculate_from_rigidity_spec(
        input_rigidity_fn, [60.0], particleName="proton"
    )
    output_df_alpha = DAFcalc.calculate_from_rigidity_spec(
        input_rigidity_fn, [60.0], particleName="alpha"
    )
    output_df_both = DAFcalc.calculate_from_rigidity_spec(
        input_rigidity_fn, [60.0], particleName="both"
    )
    _assert_additivity(output_df_proton, output_df_alpha, output_df_both)


def test_function_input_energy():
    output_df_proton = DAFcalc.calculate_from_energy_spec(
        lambda x: x**-7, [60.0], particleName="proton"
    )
    output_df_alpha = DAFcalc.calculate_from_energy_spec(
        lambda x: x**-7, [60.0], particleName="alpha"
    )
    output_df_both = DAFcalc.calculate_from_energy_spec(
        lambda x: x**-7, [60.0], particleName="both"
    )
    _assert_additivity(output_df_proton, output_df_alpha, output_df_both)


@pytest.mark.parametrize(
    "species,fixture",
    [
        ("proton", "kgo_energy_singlestep_proton"),
        ("alpha", "kgo_energy_singlestep_alpha"),
    ],
)
def test_kgo_single_step_spectrum(species, fixture):
    df = run_over_single_step_spectrum_only(species)
    _assert_df_matches_kgo(df, _load_kgo(fixture))


@pytest.mark.parametrize(
    "species,fixture",
    [
        ("proton", "kgo_energy_flat_proton"),
        ("alpha", "kgo_energy_flat_alpha"),
        ("both", "kgo_energy_flat_both"),
    ],
)
def test_kgo_flat_spectrum(species, fixture):
    df = run_over_flat_spectrum(species)
    _assert_df_matches_kgo(df, _load_kgo(fixture))


@pytest.mark.parametrize(
    "species,fixture",
    [
        ("proton", "kgo_energy_powerlaw_em7_proton"),
        ("alpha", "kgo_energy_powerlaw_em7_alpha"),
        ("both", "kgo_energy_powerlaw_em7_both"),
    ],
)
def test_kgo_energy_powerlaw_em7(species, fixture):
    """KGO: j(E) ∝ E^-7 at bin midpoints (steep cosmic-ray-like spectrum)."""
    df = run_powerlaw_em7_spectrum(species)
    _assert_df_matches_kgo(df, _load_kgo(fixture))


def test_kgo_energy_powerlaw_em2_proton():
    """KGO: softer spectrum j(E) ∝ E^-2 (strongly non-flat across bins)."""
    df = run_powerlaw_em2_proton()
    _assert_df_matches_kgo(df, _load_kgo("kgo_energy_powerlaw_em2_proton"))


@pytest.mark.parametrize(
    "species,fixture",
    [
        ("proton", "kgo_rigidity_powerlaw_rm7_proton"),
        ("alpha", "kgo_rigidity_powerlaw_rm7_alpha"),
    ],
)
def test_kgo_rigidity_powerlaw_rm7(species, fixture):
    """KGO: j(R) ∝ R^-7 at rigidity bin midpoints (species-specific rigidity grid)."""
    df = run_rigidity_powerlaw_rm7(species)
    _assert_df_matches_kgo(df, _load_kgo(fixture))


def test_additivity_flat_spectrum_all_altitudes():
    proton_df = run_over_flat_spectrum("proton")
    alpha_df = run_over_flat_spectrum("alpha")
    both_df = run_over_flat_spectrum("both")
    _assert_additivity(proton_df, alpha_df, both_df)


def test_additivity_powerlaw_em7_all_altitudes():
    proton_df = run_powerlaw_em7_spectrum("proton")
    alpha_df = run_powerlaw_em7_spectrum("alpha")
    both_df = run_powerlaw_em7_spectrum("both")
    _assert_additivity(proton_df, alpha_df, both_df)
