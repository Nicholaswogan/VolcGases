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
    # H2O solubility
    # Constants from figure table 6 in Iacono-Marziano et al. 2012. Using Anhydrous constants
    a_H2O = 0.54
    b_H2O = 1.24
    B_H2O = -2.95
    C_H2O = 0.02

    # CO2 Solubility
    # Constants from table 5 in Iacono-Marziano et al. 2012. Using anhydrous
    d_H2O = 2.3
    d_AI = 3.8
    d_FeO_MgO = -16.3
    d_Na2O_K2O = 20.1
    a_CO2 = 1
    b_CO2 = 15.8
    C_CO2 = 0.14
    B_CO2 = -5.3

    # Mass fractions of different species in Mt. Etna magma.
    # From Table 1 in Iacono-Marziano et al. 2012.
    m_SiO2 = 0.4795
    m_TiO2 = 0.0167
    m_Al2O3 = 0.1732
    m_FeO = 0.1024
    m_MgO = 0.0576
    m_CaO = 0.1093
    m_Na2O = 0.034
    m_K2O = 0.0199
    m_P2O5 = 0.0051

    # molar masses in g/mol
    M_SiO2 = 60
    M_TiO2 = 79.866
    M_Al2O3 = 101.96
    M_FeO = 71.844
    M_MgO = 40.3044
    M_CaO = 56
    M_Na2O = 61.97
    M_K2O = 94.2
    M_P2O5 = 141.94
    M_H2O = 18.01528
    M_CO2 = 44.01

    # convert mass fractions to mol fractions
    x = (m_SiO2/M_SiO2)+(m_TiO2/M_TiO2)+(m_Al2O3/M_Al2O3)+(m_FeO/M_FeO)+(m_MgO/M_MgO)+\
        (m_CaO/M_CaO)+(m_Na2O/M_Na2O)+(m_K2O/M_K2O)+(m_P2O5/M_P2O5)+(mCO2tot/M_CO2)+(mH2Otot/M_H2O)
    x_SiO2 = (m_SiO2/M_SiO2)/x
    x_TiO2 = (m_TiO2/M_TiO2)/x
    x_Al2O3 = (m_Al2O3/M_Al2O3)/x
    x_FeO = (m_FeO/M_FeO)/x
    x_MgO = (m_MgO/M_MgO)/x
    x_CaO = (m_CaO/M_CaO)/x
    x_Na2O = (m_Na2O/M_Na2O)/x
    x_K2O = (m_K2O/M_K2O)/x
    x_P2O5 = (m_P2O5/M_P2O5)/x

    # calculate NBO/O anhydrous. Appendix A in Iacono-Marziano et al. 2012
    NBO_O = (2*(x_K2O+x_Na2O+x_CaO+x_MgO+x_FeO-x_Al2O3))/ \
            (2*x_SiO2+2*x_TiO2+3*x_Al2O3+x_MgO+x_FeO+x_CaO+x_Na2O+x_K2O)

    # Calculate some constants
    x_AI = x_Al2O3/(x_CaO+x_K2O+x_Na2O)
    A1 = x_AI*d_AI+(x_FeO+x_MgO)*d_FeO_MgO\
         +(x_Na2O+x_K2O)*d_Na2O_K2O+b_CO2*NBO_O\
         +B_CO2
    A2 = b_H2O*NBO_O+B_H2O

    # calculate mol fraction of CO2 and H2O in the magma
    xCO2tot=(mCO2tot/M_CO2)/x
    xH2Otot=(mH2Otot/M_H2O)/x

    # annoying constants
    F1 = np.log(10**-6/(1))+C_CO2*P/T+A1
    F2 = np.log(1/(M_H2O*x*100))+C_H2O*P/T+A2

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
    # use different alphaG
    alphaG = .5

    # now use the solution of the simple system to solve the
    # harder problem. I will try to solve it two different ways to
    # make sure I avoid errors.
    tol = 1e-8

    try:
        init_cond = [np.log(x_H2O),np.log(x_CO2),np.log(P_H2O),np.log(P_CO2),alphaG,np.log(P_H2),np.log(P_CH4),np.log(P_CO)]
        sol = optimize.root(system,init_cond,method='lm',options={'maxiter': 10000})
        error = np.linalg.norm(system(sol['x']))
        if error>tol or sol['success']==False:
            sys.exit('Convergence issues!')
        ln_x_H2O,ln_x_CO2,ln_P_H2O,ln_P_CO2,alphaG,ln_P_H2,ln_P_CH4,ln_P_CO = sol['x']
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
