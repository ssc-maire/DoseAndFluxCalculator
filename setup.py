from setuptools import find_packages, setup
import os

# get requirements for installation
lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = lib_folder + '/requirements.txt'
install_requires = [] # Here we'll get: ["gunicorn", "docutils>=0.3", "lxml==0.5a7"]
if os.path.isfile(requirement_path):
    with open(requirement_path) as f:
        install_requires = f.read().splitlines()

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='atmosphericRadiationDoseAndFlux',
    packages=find_packages(exclude='tests'),
    package_data={"atmosphericRadiationDoseAndFlux":["atmosphericRadiationDoseAndFlux/data/proton/*.rpf","atmosphericRadiationDoseAndFlux/data/alpha/*.rpf"]},
    include_package_data=True,
    version='1.0.12',
    description='Python library for calculating doses and fluxes at a particular altitude given an input spectrum',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Space Environment and Protection Group, University of Surrey',
    url='https://github.com/ssc-maire/DoseAndFluxCalculator',
    keywords = 'atmosphere radiation space weather cosmic rays spectrum alpha particles yield response dose rates ambient effective GLE sun solar aviation neutron monitor',
    license='CC BY-NC-SA 4.0',
    install_requires=['importlib_resources>=5.10.0',
                        'numpy>=1.21.6',
                        'pandas>=1.3.5',
                        'ParticleRigidityCalculationTools>=1.5.4',
                        'setuptools>=45.2.0'],
    #install_requires,
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    test_suite='tests',
)
