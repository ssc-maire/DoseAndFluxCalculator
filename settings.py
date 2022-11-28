import os

import numpy as np

def convertListOrFloatToArray(value)->np.array:
    if isinstance(value,np.ndarray):
        outputValue = value
    elif isinstance(value,float) or isinstance(value,int):
        outputValue = np.array([value])
    else:
        outputValue = np.array(value)
    
    return outputValue

def allowCalculationForTotalOfParticles(doseFluxCalcFunc):

    def totalDoseFlux(*args, **kwargs):

        try:
            if args[2] == "both":

                doseFluxOutputList = []

                for particleName in ["proton", "alpha"]:

                    inputArgsForDFcalculator = (args[0],args[1],particleName)

                    outputDF = doseFluxCalcFunc(*inputArgsForDFcalculator, **kwargs)

                    doseFluxOutputList.append(outputDF)

                totalOutputDoseFlux = sum(doseFluxOutputList)
                totalOutputDoseFlux["altitude (km)"] = outputDF["altitude (km)"]

            else:
                totalOutputDoseFlux = doseFluxCalcFunc(*args, **kwargs)
        except IndexError:
            raise ValueError("No or incorrect input arguments supplied!")

        return totalOutputDoseFlux

    return totalDoseFlux

directory_of_this_file = os.path.dirname(os.path.realpath(__file__))
homeDirectory = directory_of_this_file
