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
```

| Input/Output |               Variable               |                 Units                  | Definition                                            |
| :- | :----------------------------------: | :------------------------------------: | :---------------------------------------------------- |
| Input |                $P$                 |                  bar                   | Total pressure of degassing                           |
| Input |                \(T\)                 |                   K                    | Temperature of magma and gas of degassing             |
| Input |         \(f_{\mathrm{O_2}}\)         |                  bar                   | Oxygen fugacity of the magma                          |
| Input | \(m_{\mathrm{CO_2}}^{\mathrm{tot}}\) | (g CO\(_2\))(g gas and magma)\(^{-1}\) | mass fraction CO\(_2\) in magma before degassing      |
| Input | \(m_{\mathrm{H_2O}}^{\mathrm{tot}}\) | (g CO\(_2\))(g gas and magma)\(^{-1}\) | mass fraction H\(_2\)O in magma before degassing      |
| Output |        \(x_{\mathrm{H_2O}}\)         |  (mol H\(_2\)O) (mol magma)\(^{-1}\)   | mol fraction of H\(_2\)O in the magma after degassing |
| Output |        \(x_{\mathrm{CO_2}}\)         |  (mol CO\(_2\)) (mol magma)\(^{-1}\)   | mol fraction of CO\(_2\) in the magma after degassing |
| Output |        \(P_{\mathrm{H_2O}}\)         |                  bar                   | Partial pressure of H\(_2\)O                          |
| Output |        \(P_{\mathrm{CO_2}}\)         |                  bar                   | Partial pressure of CO\(_2\)                          |
| Output |         \(P_{\mathrm{H_2}}\)         |                  bar                   | Partial pressure of H\(_2\)                           |
| Output |         \(P_{\mathrm{CO}}\)          |                  bar                   | Partial pressure of CO                                |
| Output |        \(P_{\mathrm{CH_4}}\)         |                  bar                   | Partial pressure of CH\(_4\)                          |
| Output |      \(\alpha_{\mathrm{gas}}\)       | (mol gas)(mol gas and magma)\(^{-1}\)  | mol fraction in gas phase                             |
