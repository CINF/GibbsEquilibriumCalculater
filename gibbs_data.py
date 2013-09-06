import atom
import molecule
import math
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np 

#from known_molecules import *
import known_molecules as km


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

    print 'manual calculating change in enthanlpy and entropy for reactions'
    dH={}
    dS={}
    dH['CO']=(km.CO.enthalpy + km.H2O.enthalpy)-(km.CO2.enthalpy+km.H2.enthalpy)
    dS['CO']=(km.CO.entropy + km.H2O.entropy)-(km.CO2.entropy+km.H2.entropy)

    dH['MeOH']=(km.CH3OH.enthalpy+km.H2O.enthalpy)-(km.CO2.enthalpy+3*km.H2.enthalpy)
    dS['MeOH']=(km.CH3OH.entropy+km.H2O.entropy)-(km.CO2.entropy+3*km.H2.entropy)

    print 'Setting up flow'
    Flow = {}
    Flow['H2'] = 8
    Flow['CO2'] = 2
    p_factor = {}
    p_factor['MeOH'] = -2.0
    p_factor['CO'] = 0.0
    pressure = 2.0 #bar
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
        import control_data
        CD = control_data.MeOH_control_set(pressure,0.75,0.25,0.0,0.0,0.0)
        
        if True:
            axis.semilogy()
            axis.plot(x,y1, 'r',label='MeOH')
            axis.plot(CD['T'],np.array(CD['MeOH'])*0.01, 'm-',label='MeOH_CD')
            axis.plot(x,y2, 'g',label='CO')
            axis.legend(loc='upper right',prop={'size':10})
            plt.xlim(500,600)
            plt.ylim(1E-5,1E-3)
            #axis.set_ylim(1E-5,1E3)
            plt.show()
            print 'Saving'
