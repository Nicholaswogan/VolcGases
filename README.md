# VolcGases
This Python program calculates the gases produced by a volcano.

## Installation
First, download or clone this repository and then navigate to the repository with a bash terminal. Finally, you can install GibbsPy with with pip command
```bash
pip install .
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
