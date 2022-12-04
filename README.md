# atmosphericRadiationDoseAndFlux

A calculation toolkit for determining radiation dose rates in the atmosphere due to an incoming spectrum of particles from space.

# Installation

After cloning this directory to your computer, install this package using the command

```
sudo python setup.py install
```

# Usage

Radiation dose rates that aircraft experience are dependent on the incoming radiation spectrum from space for each incoming particle type, 
altitude above ground level and the 'vertical cut-off rigidity' (which is a function of latitude and longitude). 

You can directly calculate radiation dose rates from an incoming rigidity spectrum using the `calculate_from_rigidity_spec` function. 
`calculate_from_rigidity_spec` takes in an input rigidity specrum as a Python function, as well as an altitude or a list of altitudes, 
and outputs the dose rates aircraft will experience at that altitude. You can optionally also include a vertical cut-off rigidity (the 
default vertical cut-off rigidity is set to 0.0 GV).

The usage of `calculate_from_rigidity_spec` is:

**calculate_from_rigidity_spec(rigiditySpectrum:Function, particleName:str, altitudesInkm:float|list, verticalCutOffRigidity:float)**

rigiditySpectrum: must be specified in units of ______
particleName:
altitudesInkm:
verticalCutOffRigidity:
                      

For example, after importing the module using
```
import atmosphericRadiationDoseAndFlux as DAFcalc
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
as a Pandas DataFrame. 

Dose rates can also be generated from an energy spectrum using `calculate_from_energy_spec`, which has the same syntax as `calculate_from_rigidity_spec`, 
but where the incoming spectrum must be specified in terms of energy (in units of ______) instead of rigidity.

