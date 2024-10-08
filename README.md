# SLME
Code and files for simple calculation of Spectroscopic Limited Maximum Efficiency (SLME) based on absorption or absorptance spectrum.

The SLME is an estimate of maximimum theoretical photovoltaic (PV) efficiency from [Yu and Zunger (2012)](10.1103/PhysRevLett.108.068701).

The code contains a class called Efficiency which, when instantiated, loads in the AM1.5G spectrum data from the CSV file AM15G.csv also in the repository. Then, the calculate method is used to calculate an efficiency. A full example of using this script is demonstrated in the Jupyter notebook in the repository.
