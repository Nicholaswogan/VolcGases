# VolcGases
This is a `Python` wrapper to a `Fortran` program that calculates the gases produced by a volcano.

Magma deep in the Earth has volatiles in it - molecules that are gases at standard temperature and pressure. As magma rises to near the surface of Earth, the overburden pressure drops. When the pressure is low enough, volatiles exsolve from the magma, and form gas bubbles (see diagram below). Chemical reactions in gas bubbles change the chemical make-up of the bubble. Eventually, the gas bubble is released from the magma to the atmosphere.

The Python function```solve_gases``` in ```VolcGases/functions.py``` calculates the gases produced by volcanoes by accounting for gases exsolving from magma to form bubbles, and the chemical reactions in those bubbles.

<p align="center">
  <img src="images/diagram.jpg" alt="Diagram" width="500"/>
</p>

## Installation
First, download or clone this repository and then navigate to the repository with a bash terminal. Finally, you can install VolcGases with with pip command. For this to work you must have `python`, and a fortran compiler. On Mac you can install a fortran compiler with the terminal command `brew install gcc`.
```bash
python -m pip install .
```

## Usage
```python
from VolcGases.functions import solve_gases

### INPUTS ###
# total mass fraction of CO2 and H2O in magma before degassing
mCO2tot = 1000e-6 # mass fraction
mH2Otot = 1000e-6 # mass fraction

# temperature and total pressure
T = 1473 # kelvin
P = 1 # bar

# Oxygen fugacity
# set to FMQ (fayalite-magnetite-quartz) mineral redox buffer.
A = 25738
B = 9
C = 0.092
log_FMQ = (-A/T+B+C*(P-1)/T)
f_O2 = 10**(log_FMQ) # units - bar

### OUTPUT ###
P_H2O,P_H2,P_CO2,P_CO,P_CH4,alphaG,x_CO2,x_H2O = solve_gases(T,P,f_O2,mCO2tot,mH2Otot)
# P_H2O = Partial pressure of H2O in gas (bar)
# P_H2 = Partial pressure of H2 in gas (bar)
# P_CO2 = Partial pressure of CO2 in gas (bar)
# P_CO = Partial pressure of CO in gas (bar)
# P_CH4 = Partial pressure of CH4 in gas (bar)
# alphaG = (mol gas)/(mol gas and magma)
# x_CO2 = mol fraction of CO2 in magma after degassing
# x_H2O = mol fraction of H2O in magma after degassing
```

Also, there is a version of `solve_gases` which is jit-ed with numba: `VolcGases.functions.solve_gases_jit`. You can use this version of solve_gases from within other numba functions.


## Other

The directory `examples` contains an example calculation. Also, the directory `article_calculations` contains calculations for a science article I wrote: https://iopscience.iop.org/article/10.3847/PSJ/abb99e/meta
