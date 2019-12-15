import numpy as np
from scipy import optimize
import sys


def solve_gases(T,P,f_O2,mCO2tot,mH2Otot):
    '''
    This function solves for the speciation of gases produced by
    a volcano. This code assumes magma composition of the lava erupting at
    Mt. Etna Italy.

    Inputs:
    T = temperature of the magma and gas in kelvin
    P = pressure of the gas in bar
    f_O2 = oxygen fugacity of the melt
    mCO2tot = mass fraction of CO2 in the magma
    mH2Otot = mass fraction of H2O in the magma

    Outputs:
    an array which contains
    [PH2O, PH2, PCO2, PCO, PCH2, alphaG, xCO2, xH2O]
    where
    PH2O = partial pressure of H2O in the gas in bar
    PH2 = partial pressure of H2 in the gas in bar
    PCO2 = partial pressure of CO2 in the gas in bar
    PCO = partial pressure of CO in the gas in bar
    PCH4 = partial pressure of CH4 in the gas in bar
    alphaG = moles of gas divide by total moles in gas and magma combined
    xCO2 = mol fraction of the CO2 in the magma
    xH2O = mol fraction of the H2O in the magma
    '''

    ###### Solubility constants
    F1 = -14.234368891317805
    F2 = -5.925014547418225
    a_H2O = 0.54
    a_CO2 = 1
    d_H2O = 2.3

    # mol of magma/g of magma
    x = 0.01550152865954013

    # molar mass in g/mol
    M_H2O = 18.01528
    M_CO2 = 44.01

    # calculate mol fraction of CO2 and H2O in the magma
    xCO2tot=(mCO2tot/M_CO2)/x
    xH2Otot=(mH2Otot/M_H2O)/x

    # equilibrium constants
    # made with Nasa thermodynamic database (Burcat database)
    K1 = np.exp(-29755.11319228574/T+6.652127716162998)
    K2 = np.exp(-33979.12369002451/T+10.418882755464773)
    K3 = np.exp(-96444.47151911151/T+0.22260815074146403)

    #constants
    C1 = K1/f_O2**0.5
    C2 = K2/f_O2**0.5
    C3 = K3/f_O2**2

    # system of equations
    def system(y):
        ln_x_H2O,ln_x_CO2,ln_H2O,ln_CO2,alphaG,ln_H2,ln_CH4,ln_CO = y
        return (np.exp(ln_H2O)+np.exp(ln_CO2)+np.exp(ln_H2)+np.exp(ln_CH4)+np.exp(ln_CO)-P,\
                -ln_x_CO2+np.exp(ln_x_H2O)*d_H2O+a_CO2*ln_CO2+F1,\
                -ln_x_H2O+a_H2O*ln_H2O+F2,\
                -xH2Otot*P + (np.exp(ln_H2O)+np.exp(ln_H2)+2*np.exp(ln_CH4))*alphaG+(1-alphaG)*np.exp(ln_x_H2O)*P,\
                -xCO2tot*P + (np.exp(ln_CO2)+np.exp(ln_CO)+np.exp(ln_CH4))*alphaG+(1-alphaG)*np.exp(ln_x_CO2)*P,\
                np.log(C1)+ln_H2O-ln_H2,\
                np.log(C2)+ln_CO2-ln_CO,\
                np.log(C3)+ln_CO2+2*ln_H2O-ln_CH4)
    # the system of equtions written slightly differently
    def system1(y):
        ln_x_H2O,ln_x_CO2,ln_H2O,ln_CO2,lnalphaG,ln_H2,ln_CH4,ln_CO = y
        return (np.exp(ln_H2O)+np.exp(ln_CO2)+np.exp(ln_H2)+np.exp(ln_CH4)+np.exp(ln_CO)-P,\
                -ln_x_CO2+np.exp(ln_x_H2O)*d_H2O+a_CO2*ln_CO2+F1,\
                -ln_x_H2O+a_H2O*ln_H2O+F2,\
                -xH2Otot*P + (np.exp(ln_H2O)+np.exp(ln_H2)+2*np.exp(ln_CH4))*np.exp(lnalphaG)+(1-np.exp(lnalphaG))*np.exp(ln_x_H2O)*P,\
                -xCO2tot*P + (np.exp(ln_CO2)+np.exp(ln_CO)+np.exp(ln_CH4))*np.exp(lnalphaG)+(1-np.exp(lnalphaG))*np.exp(ln_x_CO2)*P,\
                np.log(C1)+ln_H2O-ln_H2,\
                np.log(C2)+ln_CO2-ln_CO,\
                np.log(C3)+ln_CO2+2*ln_H2O-ln_CH4)

    # simple system
    def equation(PCO2):
        return ( -1 * ( np.e )**( ( F1 + a_CO2 * np.log( PCO2 ) ) ) * ( P )**( -1 ) \
                * ( P + ( -1 * PCO2 + -1 * C2 * PCO2 ) ) * ( ( 1 + ( C1 + C3 * PCO2 ) \
                ) )**( -1 ) + ( ( P )**( -1 ) * ( P + ( -1 * PCO2 + -1 * C2 * PCO2 ) \
                ) * ( ( 1 + ( C1 + C3 * PCO2 ) ) )**( -1 ) * xCO2tot + ( -1 * ( np.e \
                )**( ( F2 + a_H2O * np.log( (P-PCO2-C2*PCO2)/(1+C1+C3*PCO2) ) ) ) * ( -1 * ( P )**( -1 ) * PCO2 + \
                xCO2tot ) + ( ( np.e )**( ( F1 + a_CO2 * np.log( PCO2 ) ) ) * xH2Otot \
                + -1 * ( P )**( -1 ) * PCO2 * xH2Otot ) ) ) )

    # Now solve simple system
    find = np.logspace(np.log10(P),-15,1000)
    find2 = np.logspace(-15,np.log10(P),1000)
    for i in range(0,len(find)):
        if np.isnan(equation(find[i]))==False:
            found_high = find[i]
            break
    for i in range(0,len(find2)):
        if np.isnan(equation(find2[i]))==False:
            found_low = find2[i]
            break
    try:
        sol = optimize.root(equation,found_high,method='lm')
        if sol['success']==False:
            sys.exit()
    except:
        sol = optimize.root(equation,found_low,method='lm')
        if sol['success']==False:
            sys.exit('Convergence issues! First optimization.')

    P_CO2 = sol['x'][0]
    # Solve for the rest of the variables in the simple system
    P_CO = C2*P_CO2
    x_CO2 = np.exp(F1+a_CO2*np.log(P_CO2))
    alphaG = P*(x_CO2-xCO2tot)/(-P_CO2+P*x_CO2)
    P_H2O = (P-P_CO2-C2*P_CO2)/(1+C1+C3*P_CO2)
    P_H2 = C1*P_H2O
    P_CH4 = C3*P_CO2*P_H2O**2
    x_H2O = np.exp(F2+a_H2O*np.log(P_H2O))
    # use different alphaG as inital guess
    alphaG = .5

    # now use the solution of the simple system to solve the
    # harder problem. I will try to solve it two different ways to
    # make sure I avoid errors.

    # error tolerance
    tol = 1e-7

    try:
        init_cond = [np.log(x_H2O),np.log(x_CO2),np.log(P_H2O),np.log(P_CO2),alphaG,np.log(P_H2),np.log(P_CH4),np.log(P_CO)]
        sol = optimize.root(system,init_cond,method='lm',options={'maxiter': 10000})
        error = np.linalg.norm(system(sol['x']))
        if error>tol or sol['success']==False:
            sys.exit('Convergence issues!')

        ln_x_H2O,ln_x_CO2,ln_P_H2O,ln_P_CO2,alphaG,ln_P_H2,ln_P_CH4,ln_P_CO = sol['x']

        if alphaG<0:
            sys.exit('alphaG is negative')
    except:
        init_cond1 = [np.log(x_H2O),np.log(x_CO2),np.log(P_H2O),np.log(P_CO2),np.log(alphaG),np.log(P_H2),np.log(P_CH4),np.log(P_CO)]
        sol1 = optimize.root(system1,init_cond1,method='lm',options={'maxiter': 10000})
        error1 = np.linalg.norm(system1(sol1['x']))
        if error1>tol or sol1['success']==False:
            print(error)
            print(sol1)
            sys.exit('Convergence issues!')
        ln_x_H2O,ln_x_CO2,ln_P_H2O,ln_P_CO2,ln_alphaG,ln_P_H2,ln_P_CH4,ln_P_CO = sol1['x']
        alphaG = np.exp(ln_alphaG)

    return (np.exp(ln_P_H2O),np.exp(ln_P_H2),np.exp(ln_P_CO2),np.exp(ln_P_CO),\
           np.exp(ln_P_CH4),alphaG,np.exp(ln_x_CO2),np.exp(ln_x_H2O))