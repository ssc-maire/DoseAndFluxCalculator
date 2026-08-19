"""Dose and flux regression tests for atmosphericRadiationDoseAndFlux.

KGO rows are stored in DOSE_COLUMNS order:
adose, edose, dosee, tn1, tn2, tn3, SEU, SEL.
Energy inputs are MeV/n; rigidity inputs are total GV.
"""

from atmosphericRadiationDoseAndFlux.doseAndFluxCalculator import calculate_from_energy_spec_array
import numpy as np
import pytest
import ParticleRigidityCalculationTools as PRCT
import atmosphericRadiationDoseAndFlux.doseAndFluxCalculator as DAFcalc
from atmosphericRadiationDoseAndFlux.particle import Particle
from tests.kgo_data import (
    DOSE_COLUMNS,
    ENERGY_POWER_LAW_KGO,
    FLAT_ENERGY_KGO,
    GLE_LIKE_RIGIDITY_KGO,
    RIGIDITY_POWER_LAW_KGO,
    SINGLE_STEP_ENERGY_KGO,
)


KGO_RTOL = 1e-7
KGO_ATOL = 1e-16
KGO_ALTITUDES_KM = [0.0, 18.5928, 60.0]  # ground, ~FL360, and 60 km
PARTICLE_NAMES = ["proton", "alpha"]
KGO_PARTICLE_NAMES = ["proton", "alpha", "both"]


def _power_law_spectrum(x):
    return x ** -7


def _gle_like_rigidity_spectrum(rigidity_gv):
    return 27593.36 * (rigidity_gv ** -2.82844)


def _maire_energy_bins_MeV_n():
    return 10**(0.1*(np.array(range(1, 52))-1)+1)


def assert_dose_rows_match_kgo(output_df, kgo_by_altitude, rtol=KGO_RTOL, atol=KGO_ATOL):
    """Compare each stored altitude row to the matching output row."""
    for altitude_km, expected_values in kgo_by_altitude.items():
        matches = output_df[np.isclose(output_df["altitude (km)"].to_numpy(dtype=float), altitude_km)]
        assert len(matches) == 1, f"expected one row at {altitude_km} km, got {len(matches)}"
        np.testing.assert_allclose(
            matches.iloc[0][DOSE_COLUMNS].to_numpy(dtype=float),
            expected_values,
            rtol=rtol,
            atol=atol,
        )


def assert_species_sum_to_both(proton_df, alpha_df, both_df, rtol=1e-12, atol=0.0):
    """particleName='both' must equal proton + alpha."""
    for column in DOSE_COLUMNS:
        np.testing.assert_allclose(
            proton_df[column].to_numpy(dtype=float) + alpha_df[column].to_numpy(dtype=float),
            both_df[column].to_numpy(dtype=float),
            rtol=rtol,
            atol=atol,
        )


def assert_seu_sel_derived_from_tn2(output_df):
    """SEU and SEL are scaled >10 MeV neutron flux."""
    np.testing.assert_allclose(output_df["SEU"], output_df["tn2"] * 1e-13, rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(output_df["SEL"], output_df["tn2"] * 1e-8, rtol=1e-12, atol=0.0)


def getDefaultInputParameters():
    inputEnergyBins = _maire_energy_bins_MeV_n()
    inputLatitudes = np.full(60, 0.0)
    inputLongitudes = np.full(60, 0.0)
    inputAltitudes = np.linspace(0, 60, 60)
    return inputEnergyBins, inputLatitudes, inputLongitudes, inputAltitudes


def runOverSingleStepSpectrumOnly(particleSpecies: str):
    """Energy-array path: flux only in the first MeV/n bin, 60 altitude layers."""
    inputEnergyBins, inputLatitudes, inputLongitudes, inputAltitudes = getDefaultInputParameters()
    inputFluxes = np.append(1, np.full(49, 0))
    return calculate_from_energy_spec_array(
        inputEnergyBins, inputFluxes, inputAltitudes, particleName=particleSpecies
    )


def runOverFlatSpectrum(particleSpecies: str):
    """Energy-array path: unit flux in every MeV/n bin, 60 altitude layers."""
    inputEnergyBins, inputLatitudes, inputLongitudes, inputAltitudes = getDefaultInputParameters()
    inputFluxes = np.full(50, 1)
    return calculate_from_energy_spec_array(
        inputEnergyBins, inputFluxes, inputAltitudes, particleName=particleSpecies
    )


def run_energy_power_law(particle_name):
    """Energy-function path: j(E) = E**-7."""
    return DAFcalc.calculate_from_energy_spec(
        _power_law_spectrum, KGO_ALTITUDES_KM, particleName=particle_name
    )


def run_rigidity_power_law(particle_name):
    """Rigidity-function path: j(R) = R**-7."""
    return DAFcalc.calculate_from_rigidity_spec(
        _power_law_spectrum, KGO_ALTITUDES_KM, particleName=particle_name
    )


def run_gle_like_rigidity(particle_name):
    """Rigidity-function path: GLE-like power law."""
    return DAFcalc.calculate_from_rigidity_spec(
        _gle_like_rigidity_spectrum, KGO_ALTITUDES_KM, particleName=particle_name
    )


SPECTRUM_CASES = [
    pytest.param(runOverSingleStepSpectrumOnly, SINGLE_STEP_ENERGY_KGO, id="single_step_energy"),
    pytest.param(runOverFlatSpectrum, FLAT_ENERGY_KGO, id="flat_energy"),
    pytest.param(run_energy_power_law, ENERGY_POWER_LAW_KGO, id="energy_power_law"),
    pytest.param(run_rigidity_power_law, RIGIDITY_POWER_LAW_KGO, id="rigidity_power_law"),
    pytest.param(run_gle_like_rigidity, GLE_LIKE_RIGIDITY_KGO, id="gle_like_rigidity"),
]


@pytest.mark.parametrize("run_spectrum, kgo", SPECTRUM_CASES)
@pytest.mark.parametrize("particle_name", KGO_PARTICLE_NAMES)
def test_spectrum_matches_kgo(particle_name, run_spectrum, kgo):
    """Pin dose/flux outputs for each spectrum and particle type."""
    output_df = run_spectrum(particle_name)
    assert_dose_rows_match_kgo(output_df, kgo[particle_name])
    assert_seu_sel_derived_from_tn2(output_df)


@pytest.mark.parametrize("run_spectrum, _kgo", SPECTRUM_CASES)
def test_spectrum_proton_and_alpha_sum_to_both(run_spectrum, _kgo):
    """particleName='both' is the sum of the proton and alpha runs."""
    proton_df = run_spectrum("proton")
    alpha_df = run_spectrum("alpha")
    both_df = run_spectrum("both")
    assert_species_sum_to_both(proton_df, alpha_df, both_df)


def test_default_alpha_rigidity_bins_use_total_kinetic_energy():
    """At the same MeV/n, alpha total rigidity is twice the proton value (A/Z = 2)."""
    energy_bins = _maire_energy_bins_MeV_n()
    proton = Particle("proton")
    alpha = Particle("alpha")

    proton_rigidity_bins = np.array(PRCT.convertPerNucleonEnergyToTotalRigidity(
        energy_bins,
        particleMassAU=proton.atomicMass,
        particleChargeAU=proton.atomicCharge,
    ))
    alpha_rigidity_bins = np.array(PRCT.convertPerNucleonEnergyToTotalRigidity(
        energy_bins,
        particleMassAU=alpha.atomicMass,
        particleChargeAU=alpha.atomicCharge,
    ))
    naive_alpha_rigidity_bins = np.array(PRCT.convertParticleEnergyToRigidity(
        energy_bins,
        particleMassAU=alpha.atomicMass,
        particleChargeAU=alpha.atomicCharge,
    ))

    np.testing.assert_allclose(alpha_rigidity_bins, 2.0 * proton_rigidity_bins, rtol=1e-12)
    assert np.all(alpha_rigidity_bins > naive_alpha_rigidity_bins)


@pytest.mark.parametrize("particle_name", PARTICLE_NAMES)
def test_energy_and_rigidity_paths_are_consistent(particle_name):
    """The same physical spectrum should give the same dose via E/n or R."""
    energy_bins = _maire_energy_bins_MeV_n()
    energy_midpoints = (energy_bins[1:] + energy_bins[:-1]) / 2
    fluxes_per_mev_n = energy_midpoints ** -7
    particle = Particle(particle_name)

    energy_dose = calculate_from_energy_spec_array(
        energy_bins,
        fluxes_per_mev_n,
        KGO_ALTITUDES_KM,
        particleName=particle_name,
    )
    rigidity_bins = np.array(PRCT.convertPerNucleonEnergyToTotalRigidity(
        energy_bins,
        particleMassAU=particle.atomicMass,
        particleChargeAU=particle.atomicCharge,
    ))
    rigidity_spec = PRCT.convertPerNucleonEnergySpecToTotalRigiditySpec(
        energy_midpoints,
        fluxes_per_mev_n,
        particleMassAU=particle.atomicMass,
        particleChargeAU=particle.atomicCharge,
    )
    rigidity_dose = DAFcalc.calculate_from_rigidity_spec_array(
        rigidity_bins,
        rigidity_spec["Rigidity distribution values"].to_numpy(),
        KGO_ALTITUDES_KM,
        particleName=particle_name,
    )

    for column in DOSE_COLUMNS:
        np.testing.assert_allclose(
            energy_dose[column].to_numpy(dtype=float),
            rigidity_dose[column].to_numpy(dtype=float),
            rtol=5e-3,
        )
