from setuptools import find_packages, setup
import os

# get requirements for installation
lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = lib_folder + '/requirements.txt'
install_requires = [] # Here we'll get: ["gunicorn", "docutils>=0.3", "lxml==0.5a7"]
if os.path.isfile(requirement_path):
    with open(requirement_path) as f:
        install_requires = f.read().splitlines()

setup(
    name='atmosphericRadiationDoseAndFlux',
    packages=find_packages(exclude='tests'),
    package_data={"atmosphericRadiationDoseAndFlux":["atmosphericRadiationDoseAndFlux/data/proton/*.rpf","atmosphericRadiationDoseAndFlux/data/alpha/*.rpf"]},
    include_package_data=True,
    version='0.2.1',
    description='Python library for calculating doses and fluxes at a particular altitude given an input spectrum',
    author='Me',
    license='MIT',
    install_requires=install_requires,
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    test_suite='tests',
)
