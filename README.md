# Rydberg Numerical model for N-level cascade Systems for Atomic RF Sensing applications

Code Author: Liam W Bussey

This is a collabrative piece of work between BT and University of Birmingham by the sponsorship of the code author from the Royal Commission for the Exhibition of 1851 as part of their industrial fellowship.

The files contained in the Github have been used in the publication of xxxx 

## Python code 

### Background

This Python Script contained in the file labeled Rydberg_Code. The main program file is called Rydberg_Cs_example_paper_natural_decay_multiple_lines.py and is used to solve the Master equation for an N-level cascade system were our atomic system interacts with either laser or RF fields. This file takes in user information on the system for atom species, energy level tranistion, laser and RF parameters to solve the Master equation.

The Master equation:

 $$\frac{\partial \hat{\rho}}{\partial t} = -\frac{i}{\hbar}\left[\hat{H}, \hat{\rho}\right] + \hat{\mathcal{L}}(\hat{\rho})$$

where $\hat{H}$ is the Hamiltonian of the atom-light system, $\hat{\rho}$ is the density matrix and $\hat{\mathcal{L}}$ describes the decay and dephasing in the system. 



 ### assumption used


## Python dependencies - libaries and modules required to make this code run 

- sys https://docs.python.org/3/library/sys.html
- qutip https://qutip.org/
- ARC https://arc-alkali-rydberg-calculator.readthedocs.io/en/latest/
- pandas https://pandas.pydata.org/
- numpy https://numpy.org/
- Scipy https://scipy.org/
- matlibplot https://matplotlib.org/stable/tutorials/introductory/pyplot.html


## user notes



