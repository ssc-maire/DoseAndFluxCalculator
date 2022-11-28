from doseAndFluxCalculator import calculate_from_energy_spectrum
import numpy as np

def getDefaultInputParameters():

    inputEnergyBins = 10**(0.1*(np.array(range(1,52))-1)+1)

    inputLatitudes = np.full(60,0.0)
    inputLongitudes = np.full(60,0.0)
    inputAltitudes = np.linspace(0,60,60)

    return inputEnergyBins,inputLatitudes,inputLongitudes,inputAltitudes

def runOverSingleStepSpectrumOnly(particleSpecies:str):

    inputEnergyBins, inputLatitudes, inputLongitudes, inputAltitudes = getDefaultInputParameters()

    inputFluxes = np.append(1,np.full(49,0))

    outputDF = calculate_from_energy_spectrum(inputEnergyBins,inputFluxes,inputAltitudes,particleSpecies)

    print(outputDF)

def runOverFlatSpectrum(particleSpecies:str):

    inputEnergyBins, inputLatitudes, inputLongitudes, inputAltitudes = getDefaultInputParameters()

    inputFluxes = np.full(50,1)

    outputDF = calculate_from_energy_spectrum(inputEnergyBins,inputFluxes,inputAltitudes,particleSpecies)

    print(outputDF)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def test_comparison_to_original_DAF_proton_SingleStep():

    runOverSingleStepSpectrumOnly("proton")


def test_comparison_to_original_DAF_alpha_SingleStep():

    runOverSingleStepSpectrumOnly("alpha")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def test_comparison_to_original_DAF_proton_Flat():

    runOverFlatSpectrum("proton")


def test_comparison_to_original_DAF_alpha_Flat():

    runOverFlatSpectrum("alpha")