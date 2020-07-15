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
    A1 = -0.4200250000201988
    A2 = -2.59560737813789
    M_H2O = 18.01528
    M_CO2 = 44.01
    C_CO2 = 0.14
    C_H2O = 0.02
    # mol of magma/g of magma
    x = 0.01550152865954013
    F1 = np.log(1/(M_CO2*x*10**6))+C_CO2*P/T+A1
    F2 = np.log(1/(M_H2O*x*100))+C_H2O*P/T+A2
    a_H2O = 0.54
    a_CO2 = 1
    d_H2O = 2.3

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

    def jacob(y):
        lnxH2O,lnxCO2,lnPH2O,lnPCO2,alphaG,lnPH2,lnPCH4,lnPCO = y
        return np.array( [np.array( [0,0,( np.e )**( lnPH2O ),( np.e )**( lnPCO2 \
                ),0,( np.e )**( lnPH2 ),( np.e )**( lnPCH4 ),( np.e )**( lnPCO ),] \
                ),np.array( [d_H2O * ( np.e )**( lnxH2O ),-1,0,a_CO2,0,0,0,0,] \
                ),np.array( [-1,0,a_H2O,0,0,0,0,0,] ),np.array( [( 1 + -1 * alphaG ) * \
                ( np.e )**( lnxH2O ) * P,0,alphaG * ( np.e )**( lnPH2O ),0,( 2 * ( \
                np.e )**( lnPCH4 ) + ( ( np.e )**( lnPH2 ) + ( ( np.e )**( lnPH2O ) + \
                -1 * ( np.e )**( lnxH2O ) * P ) ) ),alphaG * ( np.e )**( lnPH2 ),2 * \
                alphaG * ( np.e )**( lnPCH4 ),0,] ),np.array( [0,( 1 + -1 * alphaG ) \
                * ( np.e )**( lnxCO2 ) * P,0,alphaG * ( np.e )**( lnPCO2 ),( ( np.e \
                )**( lnPCH4 ) + ( ( np.e )**( lnPCO ) + ( ( np.e )**( lnPCO2 ) + -1 * \
                ( np.e )**( lnxCO2 ) * P ) ) ),0,alphaG * ( np.e )**( lnPCH4 ),alphaG \
                * ( np.e )**( lnPCO ),] ),np.array( [0,0,1,0,0,-1,0,0,] ),np.array( \
                [0,0,0,1,0,0,0,-1,] ),np.array( [0,0,2,1,0,0,-1,0,] ),] )

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

    def jacob1(y):
        ln_x_H2O,ln_x_CO2,ln_H2O,ln_CO2,lnalphaG,ln_H2,ln_CH4,ln_CO = y
        return np.array( [np.array( [0,0,( np.e )**( ln_H2O ),( np.e )**( ln_CO2 \
                ),0,( np.e )**( ln_H2 ),( np.e )**( ln_CH4 ),( np.e )**( ln_CO ),] \
                ),np.array( [d_H2O * ( np.e )**( ln_x_H2O ),-1,0,a_CO2,0,0,0,0,] \
                ),np.array( [-1,0,a_H2O,0,0,0,0,0,] ),np.array( [( np.e )**( ln_x_H2O ) \
                * ( 1 + -1 * ( np.e )**( lnalphaG ) ) * P,0,( np.e )**( ( lnalphaG + \
                ln_H2O ) ),0,( ( np.e )**( ( lnalphaG + ln_H2O ) ) + -1 * ( np.e )**( \
                ( lnalphaG + ln_x_H2O ) ) * P ),0,0,0,] ),np.array( [0,( np.e )**( \
                ln_x_CO2 ) * ( 1 + -1 * ( np.e )**( lnalphaG ) ) * P,0,( np.e )**( ( \
                lnalphaG + ln_CO2 ) ),( ( np.e )**( ( lnalphaG + ln_CO2 ) ) + -1 * ( \
                np.e )**( ( lnalphaG + ln_x_CO2 ) ) * P ),0,0,0,] ),np.array( \
                [0,0,1,0,0,-1,0,0,] ),np.array( [0,0,0,1,0,0,0,-1,] ),np.array( \
                [0,0,2,1,0,0,-1,0,] ),] )

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
    alphaG = .1

    # now use the solution of the simple system to solve the
    # harder problem. I will try to solve it two different ways to
    # make sure I avoid errors.

    # error tolerance
    tol = 1e-7

    try:
        init_cond = [np.log(x_H2O),np.log(x_CO2),np.log(P_H2O),np.log(P_CO2),alphaG,np.log(P_H2),np.log(P_CH4),np.log(P_CO)]
        sol = optimize.root(system,init_cond,method='lm',options={'maxiter': 10000},jac=jacob)
        error = np.linalg.norm(system(sol['x']))
        if error>tol or sol['success']==False:
            sys.exit('Convergence issues!')

        ln_x_H2O,ln_x_CO2,ln_P_H2O,ln_P_CO2,alphaG,ln_P_H2,ln_P_CH4,ln_P_CO = sol['x']

        if alphaG<0:
            sys.exit('alphaG is negative')
    except:
        alphaG = .1
        init_cond1 = [np.log(x_H2O),np.log(x_CO2),np.log(P_H2O),np.log(P_CO2),np.log(alphaG),np.log(P_H2),np.log(P_CH4),np.log(P_CO)]
        sol1 = optimize.root(system1,init_cond1,method='lm',options={'maxiter': 10000},jac=jacob1)
        error1 = np.linalg.norm(system1(sol1['x']))
        ln_x_H2O,ln_x_CO2,ln_P_H2O,ln_P_CO2,ln_alphaG,ln_P_H2,ln_P_CH4,ln_P_CO = sol1['x']
        alphaG = np.exp(ln_alphaG)
        if (error1>tol and alphaG>1e-4) or sol1['success']==False:
            print(error1)
            print(sol1)
            sys.exit('Convergence issues!')
        if error1>tol:
            # assume no outgassing happens here
            return (0,0,0,0,\
                    0,0,xCO2tot,xH2Otot)
            # print('warning: outgassing equations not solved to high tolerance')

    return (np.exp(ln_P_H2O),np.exp(ln_P_H2),np.exp(ln_P_CO2),np.exp(ln_P_CO),\
           np.exp(ln_P_CH4),alphaG,np.exp(ln_x_CO2),np.exp(ln_x_H2O))

def degassing_pressure(T,DFMQ,mCO2tot,mH2Otot,P_range = [1e-4,30000]):
    """
    This function determines the overburden pressure where degassing begins.

    Inputs:
    T = temperature of the magma and gas in kelvin
    f_O2 = oxygen fugacity of the melt
    mCO2tot = mass fraction of CO2 in the magma
    mH2Otot = mass fraction of H2O in the magma
    P_range = the degassing pressure must fall within this range

    Outputs:
    P = Pressure where degassing begins
    """
    A = 25738
    B = 9
    C = 0.092
    def func(P):
        log_FMQ = (-A/T+B+C*(P-1)/T)
        f_O2 = 10**(log_FMQ+DFMQ) # units - bar
        P_H2O,P_H2,P_CO2,P_CO,P_CH4,alphaG,x_CO2,x_H2O = solve_gases(T,P,f_O2,mCO2tot,mH2Otot)
        return alphaG
    P = find_first_zero(func, P_range[0], P_range[1])
    return P

def find_first_zero(func, min, max, tol=1e-5):
    min, max = float(min), float(max)
    assert (max + tol) > max
    while (max - min) > tol:
        mid = (min + max) / 2
        sol = func(mid)
        if sol == 0 or sol != sol:
            max = mid
        else:
            min = mid
    return max
