# atmosphericRadiationDoseAndFlux

A calculation toolkit for determining radiation dose rates in the atmosphere due to an incoming spectrum of particles from space.

# Installation

After cloning this directory to your computer, install this package from the cloned directory using the command

```
sudo pip3 install .
```
or alternatively,
```
sudo python3 setup.py install
```

# Usage

## Calculating dose rates from an input spectrum as a function

Radiation dose rates that aircraft experience are dependent on the incoming radiation spectrum from space for each incoming particle type, 
altitude above ground level and the 'vertical cut-off rigidity' (which is a function of latitude and longitude). 

You can directly calculate radiation dose rates from an incoming rigidity spectrum using the `calculate_from_rigidity_spec` function. 
`calculate_from_rigidity_spec` takes in an input rigidity specrum as a Python function, as well as an altitude or a list of altitudes, 
and outputs the dose rates aircraft will experience at that altitude. You can optionally also include a vertical cut-off rigidity (the 
default vertical cut-off rigidity is set to 0.0 GV).

The usage of `calculate_from_rigidity_spec` is:

```
calculate_from_rigidity_spec(rigiditySpectrum:Function, 
                             altitudesInkm:float|list, 
                             particleName="proton", 
                             inputRigidityBins=None,
                             verticalCutOffRigidity=0.0)
```
where the input parameters are:

| parameter | description |
| --------- | ----------- |
| rigiditySpectrum | *callable function*: must be specified in units of **particles/cm2/sr/GV/s**, and as a function of rigidity in GV. rigiditySpectrum can be specified as any callable single argument function, representing an incoming particle spectrum with any distribution. You may wish to use Python's [lambda function](https://www.w3schools.com/python/python_lambda.asp) to do this.
| altitudesInkm | *float* or *list*: if specified as a list, calculations will be performed in a loop over each of the specified altitudes and returned in the output Pandas DataFrame.
| particleName | *str*: must be either "proton" or "alpha". Currently using the setting "alpha" actually calculates dose rates associated with all particles species with atomic numbers 2 and greater rather than just for alpha particles. Default = "proton".
| inputRigidityBins | *list*: an optional input parameter, do not use unless experienced with using this library. Setting this parameter to a value overrides this module's internal application of the `rigiditySpectrum` parameter to an internal default set of energy bins for the purpose of calculating dose rates and fluxes. |
| verticalCutOffRigidity | *float*: units of gigavolts (GV). Default = 0.0.    

For example, after importing the module using
```
from atmosphericRadiationDoseAndFlux import doseAndFluxCalculator as DAFcalc
```
and defining an input rigidity spectrum as a function;
```
testA = 2.77578789
testmu = 2.82874076

testRigiditySpectrum = lambda x:testA*(x**(-testmu))
```
running
```
DAFcalc.calculate_from_rigidity_spec(testRigiditySpectrum, 
                                    particleName="proton",altitudesInkm=18.5928, verticalCutOffRigidity=0.0)
```
will return
```
   altitude (km)      edose      adose      dosee        tn1       tn2  \
0        18.5928  62.219856  41.774429  34.169926  12.323553  7.639987   

        tn3           SEU           SEL  
0  5.108492  7.639987e-13  7.639987e-08 
```
as a [Pandas DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html). Here `testRigiditySpectrum` is a power law distribution in rigidity, with a spectral index of `testmu` and a value at x = 1 of `testA`.

The outputted dose rate (or flux) labels represent the following dose rate/flux types:

|label | dose rate/flux type|
|------|--------------------|
|adose| ambient dose equivalent in µSv/hr |
|edose| effective dose in µSv/hr |
|dosee| dose equivalent in µSv/hr |
|tn1| >1 MeV neutron flux, in n/cm2/s |
|tn2| >10 MeV neutron flux, in n/cm2/s |
|tn3| >60 MeV neutron flux, in n/cm2/s |
|SEU| single event upset rate for an SRAM device in upsets/second/bit |
|SEL| single event latch-up rate for an SRAM device in latch-ups/second/device |

Dose rates can also be generated from an energy spectrum using `calculate_from_energy_spec`, which has the same syntax as `calculate_from_rigidity_spec`, 
but where the incoming spectrum (now in units of *particles/cm2/sr/MeV/s*) must be specified in terms of energy in MeV instead of rigidity and the optional `inputRigidityBins` parameter is instead `inputEnergyBins` in MeV.

## Calculating dose rates from an input spectrum as an array

While it is recommended that the user uses functions as the method of inputting spectra, you can also use the `calculate_from_rigidity_spec_array` and `calculate_from_energy_spec_array` functions to calculate dose rates directly from an array of spectral flux values and energy bins.

The general syntax for inputs to these are:

```
calculate_from_rigidity_spec_array(
                            inputRigidityBins:list,
                            inputFluxesGV:list, 
                            altitudesInkm:list, 
                            particleName="proton",
                            verticalCutOffRigidity = 0.0)
```
                            
where `inputRigidityBins` is the full list of rigidity bin edge locations in GV that the input fluxes specified by `inputFluxesGV` are specified for. As `inputRigidityBins` represent the bin edge locations for `inputFluxesGV`, the length of `inputRigidityBins` must therefore be exactly 1 greater than the length of `inputFluxesGV`.

`calculate_from_energy_spec_array` has a nearly identical syntax to `calculate_from_rigidity_spec_array`, where the only difference is that the first argument of `calculate_from_energy_spec_array` must be specified in MeV instead of GV, and the second argument is the energy spectrum in particles/cm2/sr/MeV/s instead of the rigidity spectrum in particles/cm2/sr/GV/s.





