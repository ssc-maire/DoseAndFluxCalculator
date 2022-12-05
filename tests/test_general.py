from atmosphericRadiationDoseAndFlux.doseAndFluxCalculator import calculate_from_energy_spec_array
import numpy as np
import atmosphericRadiationDoseAndFlux.doseAndFluxCalculator as DAFcalc

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

    assert int(protonSlice["tn3"].round(-4)) == 178580000
    assert int(protonSlice["adose"].round(-6).iloc[0]) == 1.974e9

    print(protonFlatSpec)


def test_comparison_to_original_DAF_alpha_Flat():

    alphaFlatSpec = runOverFlatSpectrum("alpha")

    alphaSlice = alphaFlatSpec[alphaFlatSpec["altitude (km)"] == 60.0]

    assert int(alphaSlice["tn1"].round(-4)) == int(1.84192e9)
    assert int(alphaSlice["adose"].round(-6).iloc[0]) == int(9.871e9)

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

