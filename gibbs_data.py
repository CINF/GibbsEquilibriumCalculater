import atom
import molecule
import math
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np 

entropy_data = {} # unit J/mol/K
entropy_data[molecule.Molecule([atom.Atom(1), atom.Atom(1)] )] = 130.679 #H2
entropy_data[molecule.Molecule([atom.Atom(1), atom.Atom(1), atom.Atom(8)])] = 188.84 #H2O
#entropy_data[molecule.Molecule([atom.Atom(1), atom.Atom(1), atom.Atom(8), atom.Atom(8)])] = 113 #H2O2
#entropy_data[molecule.Molecule([atom.Atom(8), atom.Atom(8)])] = 113 #O2

entropy_data[molecule.Molecule([atom.Atom(6), atom.Atom(8)])] = 197.7 #CO
entropy_data[molecule.Molecule([atom.Atom(6), atom.Atom(8), atom.Atom(8)])] = 213.7 #CO2

entropy_data[molecule.Molecule([atom.Atom(6), atom.Atom(1), atom.Atom(1), atom.Atom(1), atom.Atom(1)])] = 186.25 #CH4
entropy_data[molecule.Molecule([atom.Atom(6), atom.Atom(1), atom.Atom(1), atom.Atom(1), atom.Atom(8), atom.Atom(1)])] = 239.9 #CH3OH

entropy_data[molecule.Molecule([atom.Atom(7), atom.Atom(1), atom.Atom(1), atom.Atom(1)])] = 192.778 #NH3

entropy_data[molecule.Molecule([atom.Atom(7), atom.Atom(7)])] = 191.61 #N2
#entropy_data[molecule.Molecule([atom.Atom(7), atom.Atom(8)])] = 113 #NO
#entropy_data[molecule.Molecule([atom.Atom(7), atom.Atom(8), atom.Atom(8)])] = 113 #NO2

enthalpy_data = {} # unit J/mol
enthalpy_data[molecule.Molecule([atom.Atom(1), atom.Atom(1)] )] = 0.0 #H2
enthalpy_data[molecule.Molecule([atom.Atom(1), atom.Atom(1), atom.Atom(8)])] = -241830.0 #H2O
#enthalpy_data[molecule.Molecule([atom.Atom(1), atom.Atom(1), atom.Atom(8), atom.Atom(8)])] = 113 #H2O2
enthalpy_data[molecule.Molecule([atom.Atom(8), atom.Atom(8)])] = 0.0 #O2

enthalpy_data[molecule.Molecule([atom.Atom(6), atom.Atom(8)])] = -110500.0 #CO
enthalpy_data[molecule.Molecule([atom.Atom(6), atom.Atom(8), atom.Atom(8)])] = -393500.0 #CO2

enthalpy_data[molecule.Molecule([atom.Atom(6), atom.Atom(1), atom.Atom(1), atom.Atom(1), atom.Atom(1)])] = -74870.0 #CH4
enthalpy_data[molecule.Molecule([atom.Atom(6), atom.Atom(1), atom.Atom(1), atom.Atom(1), atom.Atom(8), atom.Atom(1)])] = -201300.0 #CH3OH

enthalpy_data[molecule.Molecule([atom.Atom(7), atom.Atom(1), atom.Atom(1), atom.Atom(1)])] = -45940.0 #NH3

enthalpy_data[molecule.Molecule([atom.Atom(7), atom.Atom(7)])] = 0.0 #N2
#enthalpy_data[molecule.Molecule([atom.Atom(7), atom.Atom(8)])] = 113 #NO
#enthalpy_data[molecule.Molecule([atom.Atom(7), atom.Atom(8), atom.Atom(8)])] = 113 #NO2


def assign_random_data(Set):
    for M in Set:
        try:
            M.entropy = entropy_data[M]
        except KeyError:
            print 'Molecule doesnt exist in DB'
        try:
            M.enthalpy = enthalpy_data[M]
        except KeyError:
            print 'Molecule doesnt exist in DB'

def K_equilibrium(dH,dS,T): # should not be in this file, just a check if it is working
    R=8.3144621 # Gas constant
    K = math.exp(-(dH-T*dS)/(R*T))
    return K

def Equilibrium_MeOH((MeOH, CO, H2O, CO2, H2),T,pressure): # should not be in this file, just a check if it is working
    K_MeOH = K_equilibrium(dH['MeOH'],dS['MeOH'],T)
    K_CO = K_equilibrium(dH['CO'],dS['CO'],T)
    return (
            K_MeOH-(((MeOH) * (H2O)) / ((CO2) * (H2)**3))*(pressure)**p_factor['MeOH'],
            K_CO-(((CO) * (H2O)) / ((CO2) * (H2)))*(pressure)**p_factor['CO'],
            (2*Flow['H2'])/(Flow['CO2'])-(4*MeOH+2*H2O+2*H2)/(MeOH+CO2+CO), #H/C
            (2*Flow['H2'])/(2*Flow['CO2'])-(4*MeOH+2*H2O+2*H2)/(MeOH+H2O+2*CO2+CO), #H/O
            1.0-(MeOH+CO+H2O+CO2+H2)
            )

if __name__ == '__main__':
    print 'Creating molecules'
    MeOH = molecule.Molecule([atom.Atom(6), atom.Atom(1), atom.Atom(1), atom.Atom(1), atom.Atom(8), atom.Atom(1)])
    CO = molecule.Molecule([atom.Atom(6), atom.Atom(8)])
    CO2 = molecule.Molecule([atom.Atom(6), atom.Atom(8), atom.Atom(8)])
    H2O = molecule.Molecule([atom.Atom(1), atom.Atom(1), atom.Atom(8)])
    H2 = molecule.Molecule([atom.Atom(1), atom.Atom(1)])
    
    print 'adding ethanlpy and entropy to molecules'
    assign_random_data([MeOH,CO,CO2,H2O,H2])
    
    print 'manual calculating change in enthanlpy and entropy for reactions'
    dH={}
    dS={}
    dH['CO']=(CO.enthalpy+H2O.enthalpy)-(CO2.enthalpy+H2.enthalpy)
    dS['CO']=(CO.entropy+H2O.entropy)-(CO2.entropy+H2.entropy)

    dH['MeOH']=(MeOH.enthalpy+H2O.enthalpy)-(CO2.enthalpy+3*H2.enthalpy)
    dS['MeOH']=(MeOH.entropy+H2O.entropy)-(CO2.entropy+3*H2.entropy)

    print 'Setting up flow'
    Flow = {}
    Flow['H2'] = 8
    Flow['CO2'] = 2
    p_factor = {}
    p_factor['MeOH'] = -2.0
    p_factor['CO'] = 0.0
    pressure = 2.5 #bar
    T=300.0 # K

    print 'Eq test'
    #print Equilibrium_MeOH((0.10, 0.001, 0.10, 0.14, 0.66),T,pressure)
    start_0 = fsolve(Equilibrium_MeOH, (0.01, 0.005, 0.015, 0.19, 0.78), args=(400, 1.0), xtol=1.49012e-12,maxfev=10000 )
    for pi in np.linspace(1.0,pressure,100):
        start_0 = fsolve(Equilibrium_MeOH, start_0, args=(400, pi), xtol=1.49012e-12,maxfev=10000 )
    for ti in np.linspace(400.0,T,100):
        start_0 = fsolve(Equilibrium_MeOH, start_0, args=(ti, pressure), xtol=1.49012e-12,maxfev=10000 )
    print start_0
    if True:
        fig = plt.figure()
        axis = fig.add_subplot(1,1,1)
        axis.tick_params(direction='in', length=6, width=1, colors='k',labelsize=10,axis='both',pad=3)
        axis.set_xlabel('T / [K]', fontsize=10)
        axis.set_ylabel('c', fontsize=10)
        x = np.zeros(300)
        y1 = np.zeros(300)
        y2 = np.zeros(300)
        i = 0
        cMeOH, cCO, cH2O,cCO2,cH2 =  fsolve(Equilibrium_MeOH, start_0,args=(T,pressure))
        

        cMeOH, cCO, cH2O,cCO2,cH2 =  fsolve(Equilibrium_MeOH, start_0,args=(T,pressure))
        print cMeOH, cCO, cH2O,cCO2,cH2
        for t in range(300,600):
            #print t
            cMeOH, cCO, cH2O,cCO2,cH2 =  fsolve(Equilibrium_MeOH, (cMeOH, cCO, cH2O,cCO2,cH2),args=(t*1.0,pressure))
            x[i] = t
            y1[i] = cMeOH
            y2[i] = cCO
            i+=1
        if True:
            axis.semilogy()
            axis.plot(x,y1, 'r',label='MeOH')
            axis.plot(x,y2, 'g',label='CO')
            axis.legend(loc='upper right',prop={'size':10})
            plt.xlim(300,600)
            #axis.set_ylim(1E-5,1E3)
            plt.show()
            print 'Saving'
