from atmosphericRadiationDoseAndFlux.doseAndFluxCalculator import calculate_from_energy_spec_array
import numpy as np
import pytest
import ParticleRigidityCalculationTools as PRCT
import atmosphericRadiationDoseAndFlux.doseAndFluxCalculator as DAFcalc
from atmosphericRadiationDoseAndFlux.particle import Particle


DOSE_COLUMNS = ["adose", "edose", "dosee", "tn1", "tn2", "tn3", "SEU", "SEL"]
KGO_RTOL = 1e-7
KGO_ATOL = 1e-16
KGO_ALTITUDES_KM = [0.0, 18.5928, 60.0]
PARTICLE_NAMES = ["proton", "alpha"]

# Known-good dose rows are stored in DOSE_COLUMNS order.

# Unit flux in the first MAIRE MeV/n bin only, on the default 60-altitude grid.
SINGLE_STEP_ENERGY_KGO = {
    "proton": {
        0.0: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        60.0: [
            0.0028298404753319883,
            0.0014273066825505071,
            0.0006361271373600208,
            0.0009164120491422661,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    },
    "alpha": {
        0.0: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        60.0: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    },
}

# Unit differential flux on every MAIRE MeV/n bin.
FLAT_ENERGY_KGO = {
    "proton": {
        0.0: [
            446940602.3287374,
            436875594.0749431,
            718367866.6535243,
            20720949.640839044,
            12793093.691877475,
            8618015.43158192,
            1.2793093691877476e-06,
            0.12793093691877475,
        ],
        60.0: [
            681326991.39099,
            1037297891.9663112,
            665877743.1070614,
            130533518.5518617,
            91324127.43872547,
            59828946.97905593,
            9.132412743872547e-06,
            0.9132412743872548,
        ],
    },
    "alpha": {
        0.0: [
            2275766153.8019385,
            2258316957.1751366,
            3630429172.2581286,
            104131954.56634524,
            62257951.62991166,
            44229399.608333275,
            6.225795162991166e-06,
            0.6225795162991167,
        ],
        60.0: [
            3320606687.3991294,
            4853360773.164932,
            3402771594.818055,
            617023751.3180094,
            437125229.91899216,
            292044450.07831144,
            4.3712522991899216e-05,
            4.371252299189922,
        ],
    },
}

# j(E) = E**-7 on the default MeV/n bins.
ENERGY_POWER_LAW_KGO = {
    "proton": {
        0.0: [
            1.1581534978513966e-19,
            8.298795858089139e-20,
            1.8242000460838533e-21,
            6.17820858593639e-20,
            2.874415832386681e-20,
            1.5808375487687443e-20,
            2.874415832386681e-33,
            2.874415832386681e-28,
        ],
        18.5928: [
            7.792349521281516e-12,
            3.2724633714611035e-12,
            1.691511472185948e-12,
            1.8140506683820501e-12,
            8.64105570333197e-14,
            3.946884835935371e-15,
            8.641055703331969e-27,
            8.64105570333197e-22,
        ],
        60.0: [
            2.737488168362968e-07,
            5.82317095212612e-09,
            2.5149506406912646e-06,
            1.0846107958241945e-10,
            6.979737210028951e-12,
            1.7169661171496722e-14,
            6.979737210028951e-25,
            6.97973721002895e-20,
        ],
    },
    "alpha": {
        0.0: [
            1.0331513839110175e-11,
            1.01345138578318e-11,
            1.727023788023973e-11,
            7.945720852814965e-13,
            4.662403197230455e-13,
            3.392822919438783e-13,
            4.6624031972304555e-26,
            4.662403197230455e-21,
        ],
        18.5928: [
            3.662209847023633e-10,
            5.12032465940558e-10,
            4.4430166347112934e-10,
            5.7966741871301606e-11,
            3.756690062951024e-11,
            2.7280370711010762e-11,
            3.756690062951024e-24,
            3.756690062951024e-19,
        ],
        60.0: [
            4.033121813911825e-06,
            6.951592554539839e-07,
            4.355887455705869e-05,
            1.9460880802306542e-11,
            1.3699850722053808e-11,
            8.499312579881725e-12,
            1.3699850722053808e-24,
            1.3699850722053808e-19,
        ],
    },
}

# j(R) = R**-7 on the default total-GV bins.
RIGIDITY_POWER_LAW_KGO = {
    "proton": {
        0.0: [
            0.0002690919639564223,
            0.00018805064473081094,
            1.87311141713132e-05,
            0.00013351957369313956,
            7.492851560204757e-05,
            4.643421827288639e-05,
            7.492851560204757e-18,
            7.492851560204757e-13,
        ],
        18.5928: [
            18.882578980454603,
            12.077659888932066,
            6.12364248886505,
            6.928738423568604,
            2.882635370856761,
            1.1920443788882045,
            2.882635370856761e-13,
            2.882635370856761e-08,
        ],
        60.0: [
            216485.34708352434,
            13533.32503720353,
            1303292.2725818453,
            94.71482962245014,
            24.60732569989683,
            2.02428919852434,
            2.460732569989683e-12,
            2.460732569989683e-07,
        ],
    },
    "alpha": {
        0.0: [
            2.1117700278998495,
            2.071488275521432,
            3.5299716605141316,
            0.16243217892286055,
            0.09531206104829884,
            0.06935690517556849,
            9.531206104829885e-15,
            9.531206104829884e-10,
        ],
        18.5928: [
            75.06268394673502,
            105.06275991388294,
            90.98733799098214,
            11.910502791954261,
            7.721633634299716,
            5.6101386564689895,
            7.721633634299716e-13,
            7.721633634299716e-08,
        ],
        60.0: [
            39470.57565727135,
            17281.28619902485,
            371569.2866504834,
            3.6313667789406123,
            2.523778612354634,
            1.5983305654040947,
            2.5237786123546345e-13,
            2.5237786123546342e-08,
        ],
    },
}

# j(R) = 27593.36 * R**-2.82844 on the default total-GV bins.
GLE_LIKE_RIGIDITY_KGO = {
    "proton": {
        0.0: [
            676.0738509748024,
            617.1965548589659,
            834.4373103591985,
            136.46294642026854,
            83.93686450656759,
            56.95395755239967,
            8.393686450656759e-12,
            8.393686450656759e-07,
        ],
        18.5928: [
            406798.7808757193,
            608980.1207140976,
            328000.87638485903,
            121537.08373355807,
            75325.21184583633,
            50340.87363911171,
            7.532521184583633e-09,
            0.0007532521184583633,
        ],
        60.0: [
            15482297.027492698,
            8006614.965156902,
            45071119.91194518,
            84856.42088991875,
            48356.3777253399,
            25028.229387252595,
            4.83563777253399e-09,
            0.000483563777253399,
        ],
    },
    "alpha": {
        0.0: [
            11404.326194010302,
            11110.152297207651,
            18625.976965655955,
            1032.7003674182074,
            611.9864547650554,
            438.293320953414,
            6.119864547650554e-11,
            6.119864547650554e-06,
        ],
        18.5928: [
            1000091.5057644332,
            1462199.85258709,
            909772.7766756713,
            271775.0374883809,
            176533.53119496605,
            124398.00939512804,
            1.7653353119496605e-08,
            0.0017653353119496606,
        ],
        60.0: [
            37467763.51137213,
            58395601.01121761,
            320534147.59327877,
            128182.36610577359,
            87065.06549276336,
            51106.99099612424,
            8.706506549276337e-09,
            0.0008706506549276337,
        ],
    },
}


def _power_law_spectrum(x):
    return x ** -7


def _gle_like_rigidity_spectrum(rigidity_gv):
    return 27593.36 * (rigidity_gv ** -2.82844)


def _maire_energy_bins_MeV_n():
    return 10**(0.1*(np.array(range(1, 52))-1)+1)


def assert_dose_rows_match_kgo(output_df, kgo_by_altitude, rtol=KGO_RTOL, atol=KGO_ATOL):
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
    for column in DOSE_COLUMNS:
        np.testing.assert_allclose(
            proton_df[column].to_numpy(dtype=float) + alpha_df[column].to_numpy(dtype=float),
            both_df[column].to_numpy(dtype=float),
            rtol=rtol,
            atol=atol,
        )


def assert_seu_sel_derived_from_tn2(output_df):
    np.testing.assert_allclose(output_df["SEU"], output_df["tn2"] * 1e-13, rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(output_df["SEL"], output_df["tn2"] * 1e-8, rtol=1e-12, atol=0.0)


def getDefaultInputParameters():
    inputEnergyBins = _maire_energy_bins_MeV_n()
    inputLatitudes = np.full(60, 0.0)
    inputLongitudes = np.full(60, 0.0)
    inputAltitudes = np.linspace(0, 60, 60)
    return inputEnergyBins, inputLatitudes, inputLongitudes, inputAltitudes


def runOverSingleStepSpectrumOnly(particleSpecies: str):
    inputEnergyBins, inputLatitudes, inputLongitudes, inputAltitudes = getDefaultInputParameters()
    inputFluxes = np.append(1, np.full(49, 0))
    return calculate_from_energy_spec_array(
        inputEnergyBins, inputFluxes, inputAltitudes, particleName=particleSpecies
    )


def runOverFlatSpectrum(particleSpecies: str):
    inputEnergyBins, inputLatitudes, inputLongitudes, inputAltitudes = getDefaultInputParameters()
    inputFluxes = np.full(50, 1)
    return calculate_from_energy_spec_array(
        inputEnergyBins, inputFluxes, inputAltitudes, particleName=particleSpecies
    )


def run_energy_power_law(particle_name):
    return DAFcalc.calculate_from_energy_spec(
        _power_law_spectrum, KGO_ALTITUDES_KM, particleName=particle_name
    )


def run_rigidity_power_law(particle_name):
    return DAFcalc.calculate_from_rigidity_spec(
        _power_law_spectrum, KGO_ALTITUDES_KM, particleName=particle_name
    )


def run_gle_like_rigidity(particle_name):
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
@pytest.mark.parametrize("particle_name", PARTICLE_NAMES)
def test_spectrum_matches_kgo(particle_name, run_spectrum, kgo):
    output_df = run_spectrum(particle_name)
    assert_dose_rows_match_kgo(output_df, kgo[particle_name])
    assert_seu_sel_derived_from_tn2(output_df)


@pytest.mark.parametrize("run_spectrum, _kgo", SPECTRUM_CASES)
def test_spectrum_proton_and_alpha_sum_to_both(run_spectrum, _kgo):
    proton_df = run_spectrum("proton")
    alpha_df = run_spectrum("alpha")
    both_df = run_spectrum("both")
    assert_species_sum_to_both(proton_df, alpha_df, both_df)


def test_default_alpha_rigidity_bins_use_total_kinetic_energy():
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
