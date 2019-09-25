# Volcano-Speciation
The function "solve_gases" contained within the functions.py folder calculates the speciation of gases produced by volcano by solving gas-melt and gas-gas equilibrium. 

Inputs:
T = temperature of the magma and gas in kelvin
P = pressure of the gas in bar
f_O2 = oxygen fugacity of the melt
mCO2tot = mass fraction of CO2 in the magma
mH2Otot = mass fraction of H2O in the magma

Outputs:
an array which contains PH2O, PH2, PCO2, PCO, PCH2, alphaG, xCO2, xH2O
where,
PH2O = partial pressure of H$_2$O in the gas in bar
PH2 = partial pressure of H2 in the gas in bar
PCO2 = partial pressure of CO2 in the gas in bar
PCO = partial pressure of CO in the gas in bar
PCH4 = partial pressure of CH4 in the gas in bar
alphaG = moles of gas divide by total moles in gas and magma combined
xCO2 = mol fraction of the CO2 in the magma
xH2O = mol fraction of the H2O in the magma
