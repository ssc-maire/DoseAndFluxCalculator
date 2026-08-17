from atmosphericRadiationDoseAndFlux.doseAndFluxCalculator import calculate_from_energy_spec_array
import numpy as np
import ParticleRigidityCalculationTools as PRCT
import atmosphericRadiationDoseAndFlux.doseAndFluxCalculator as DAFcalc
from atmosphericRadiationDoseAndFlux.particle import Particle

def getDefaultInputParameters():

    inputEnergyBins = 10**(0.1*(np.array(range(1,52))-1)+1)

    inputLatitudes = np.full(60,0.0)
    inputLongitudes = np.full(60,0.0)
    inputAltitudes = np.linspace(0,60,60)

    return inputEnergyBins,inputLatitudes,inputLongitudes,inputAltitudes

def runOverSingleStepSpectrumOnly(particleSpecies:str):

    inputEnergyBins, inputLatitudes, inputLongitudes, inputAltitudes = getDefaultInputParameters()

    inputFluxes = np.append(1,np.full(49,0))

    outputDF = calculate_from_energy_spec_array(inputEnergyBins,inputFluxes,inputAltitudes,particleName = particleSpecies)

    return outputDF

def runOverFlatSpectrum(particleSpecies:str):

    inputEnergyBins, inputLatitudes, inputLongitudes, inputAltitudes = getDefaultInputParameters()

    inputFluxes = np.full(50,1)

    outputDF = calculate_from_energy_spec_array(inputEnergyBins,inputFluxes,inputAltitudes,particleName = particleSpecies)

    return outputDF

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def test_function_input_rigidity():

    outputDFproton = DAFcalc.calculate_from_rigidity_spec(lambda x:x**-7, [60.0], particleName="proton")
    outputDFalpha = DAFcalc.calculate_from_rigidity_spec(lambda x:x**-7, [60.0], particleName="alpha")
    outputDFboth = DAFcalc.calculate_from_rigidity_spec(lambda x:x**-7, [60.0], particleName="both")

    for doseType in ["adose","edose","dosee","tn1","tn2","tn3","SEU","SEL"]:
        assert outputDFproton[doseType].iloc[0] + outputDFalpha[doseType].iloc[0] == outputDFboth[doseType].iloc[0]

    print(outputDFboth)

def test_function_input_rigidity2():

    inputRigidityFunction = lambda x:27593.36 * (x**-2.82844)

    outputDFproton = DAFcalc.calculate_from_rigidity_spec(inputRigidityFunction, [60.0], particleName="proton")
    outputDFalpha = DAFcalc.calculate_from_rigidity_spec(inputRigidityFunction, [60.0], particleName="alpha")
    outputDFboth = DAFcalc.calculate_from_rigidity_spec(inputRigidityFunction, [60.0], particleName="both")

    for doseType in ["adose","edose","dosee","tn1","tn2","tn3","SEU","SEL"]:
        assert outputDFproton[doseType].iloc[0] + outputDFalpha[doseType].iloc[0] == outputDFboth[doseType].iloc[0]

    print(outputDFboth)

def test_function_input_energy():

    outputDFproton = DAFcalc.calculate_from_energy_spec(lambda x:x**-7, [60.0], particleName="proton")
    outputDFalpha = DAFcalc.calculate_from_energy_spec(lambda x:x**-7, [60.0], particleName="alpha")
    outputDFboth = DAFcalc.calculate_from_energy_spec(lambda x:x**-7, [60.0], particleName="both")

    for doseType in ["adose","edose","dosee","tn1","tn2","tn3","SEU","SEL"]:
        assert outputDFproton[doseType].iloc[0] + outputDFalpha[doseType].iloc[0] == outputDFboth[doseType].iloc[0]

    print(outputDFboth)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def test_comparison_to_original_DAF_proton_SingleStep():

    runOverSingleStepSpectrumOnly("proton")


def test_comparison_to_original_DAF_alpha_SingleStep():

    runOverSingleStepSpectrumOnly("alpha")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def test_comparison_to_original_DAF_proton_Flat():

    protonFlatSpec = runOverFlatSpectrum("proton")

    protonSlice = protonFlatSpec[protonFlatSpec["altitude (km)"] == 60.0]

    assert int(protonSlice["tn3"].iloc[0].round(-4)) == round(59828888,-4) #178580000
    assert int(protonSlice["adose"].iloc[0].round(-6)) == round(681326272,-6) #1.974e9

    print(protonFlatSpec)


def test_comparison_to_original_DAF_alpha_Flat():

    alphaFlatSpec = runOverFlatSpectrum("alpha")

    alphaSlice = alphaFlatSpec[alphaFlatSpec["altitude (km)"] == 60.0]

    assert int(alphaSlice["tn3"].iloc[0].round(-4)) == round(292044128,-4) #int(1.84192e9)
    assert int(alphaSlice["adose"].iloc[0].round(-6)) == round(3320603140,-6) #int(9.871e9)

    print(alphaFlatSpec)

def test_comparison_to_original_DAF_both_Flat():

    protonFlatSpec = runOverFlatSpectrum("proton")
    alphaFlatSpec = runOverFlatSpectrum("alpha")
    bothFlatSpec = runOverFlatSpectrum("both")

    protonSlice = protonFlatSpec[protonFlatSpec["altitude (km)"] == 60.0]
    alphaSlice = alphaFlatSpec[alphaFlatSpec["altitude (km)"] == 60.0]
    bothSlice = bothFlatSpec[bothFlatSpec["altitude (km)"] == 60.0]

    for doseType in ["adose","edose","dosee","tn1","tn2","tn3","SEU","SEL"]:
        assert protonFlatSpec[doseType].iloc[0] + alphaFlatSpec[doseType].iloc[0] == bothFlatSpec[doseType].iloc[0]

    print(protonFlatSpec)
    print(alphaFlatSpec)
    print(bothFlatSpec)


def _maire_energy_bins_MeV_n():
    return 10**(0.1*(np.array(range(1,52))-1)+1)


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


def test_alpha_energy_and_rigidity_paths_are_consistent():
    energy_bins = _maire_energy_bins_MeV_n()
    energy_midpoints = (energy_bins[1:] + energy_bins[:-1]) / 2
    fluxes_per_mev_n = energy_midpoints ** -7
    alpha = Particle("alpha")

    energy_dose = calculate_from_energy_spec_array(
        energy_bins,
        fluxes_per_mev_n,
        [60.0],
        particleName="alpha",
    )

    rigidity_bins = np.array(PRCT.convertPerNucleonEnergyToTotalRigidity(
        energy_bins,
        particleMassAU=alpha.atomicMass,
        particleChargeAU=alpha.atomicCharge,
    ))
    rigidity_spec = PRCT.convertPerNucleonEnergySpecToTotalRigiditySpec(
        energy_midpoints,
        fluxes_per_mev_n,
        particleMassAU=alpha.atomicMass,
        particleChargeAU=alpha.atomicCharge,
    )
    rigidity_dose = DAFcalc.calculate_from_rigidity_spec_array(
        rigidity_bins,
        rigidity_spec["Rigidity distribution values"].to_numpy(),
        [60.0],
        particleName="alpha",
    )

    for dose_type in ["adose", "edose", "dosee", "tn1", "tn2", "tn3"]:
        np.testing.assert_allclose(
            energy_dose[dose_type].iloc[0],
            rigidity_dose[dose_type].iloc[0],
            rtol=5e-3,
        )


def test_proton_energy_and_rigidity_paths_remain_consistent():
    energy_bins = _maire_energy_bins_MeV_n()
    energy_midpoints = (energy_bins[1:] + energy_bins[:-1]) / 2
    fluxes_per_mev = energy_midpoints ** -7
    proton = Particle("proton")

    energy_dose = calculate_from_energy_spec_array(
        energy_bins,
        fluxes_per_mev,
        [60.0],
        particleName="proton",
    )

    rigidity_bins = np.array(PRCT.convertPerNucleonEnergyToTotalRigidity(
        energy_bins,
        particleMassAU=proton.atomicMass,
        particleChargeAU=proton.atomicCharge,
    ))
    rigidity_spec = PRCT.convertPerNucleonEnergySpecToTotalRigiditySpec(
        energy_midpoints,
        fluxes_per_mev,
        particleMassAU=proton.atomicMass,
        particleChargeAU=proton.atomicCharge,
    )
    rigidity_dose = DAFcalc.calculate_from_rigidity_spec_array(
        rigidity_bins,
        rigidity_spec["Rigidity distribution values"].to_numpy(),
        [60.0],
        particleName="proton",
    )

    for dose_type in ["adose", "edose", "dosee", "tn1", "tn2", "tn3"]:
        np.testing.assert_allclose(
            energy_dose[dose_type].iloc[0],
            rigidity_dose[dose_type].iloc[0],
            rtol=5e-3,
        )


