import atom
import molecule
import known_molecules as km
import chemical_reaction
import control_data



#import math
#from scipy.optimize import fsolve
#import matplotlib.pyplot as plt
#import numpy as np 

#from known_molecules import *




if __name__ == '__main__':
    print 'START'
    print ''
    print 'known molecules'
    mol_1 = km.H2
    T=float('Nan')
    print'H2 entropy at unknown temp: ' + str(mol_1.standard_entropy())
    print ''
    print'H2 entropy at 298.15K: ' + str(mol_1.standard_entropy(298.15))

    print ''
    print'H2 enthalpy at unknown temp: ' + str(mol_1.standard_enthalpy())
    print ''
    print'H2 enthalpy at 298.15K: ' + str(mol_1.standard_enthalpy(298.15))

    print ''
    print'H2 gibbs at unknown temp: ' + str(mol_1.standard_gibbs())
    print ''
    print'H2 gibbs at 298.15K: ' + str(mol_1.standard_gibbs(298.15))
    

    print 'manual calculating change in enthanlpy and entropy for reactions'
    dH={}
    dS={}
    dH['CO']=(km.CO.standard_enthalpy() + km.H2O.standard_enthalpy())-(km.CO2.standard_enthalpy()+km.H2.standard_enthalpy())
    dS['CO']=(km.CO.standard_entropy() + km.H2O.standard_entropy())-(km.CO2.standard_entropy()+km.H2.standard_entropy())

    dH['MeOH']=(km.CH3OH.standard_enthalpy()+km.H2O.standard_enthalpy())-(km.CO2.standard_enthalpy()+3*km.H2.standard_enthalpy())
    dS['MeOH']=(km.CH3OH.standard_entropy()+km.H2O.standard_entropy())-(km.CO2.standard_entropy()+3*km.H2.standard_entropy())

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
            plt.xlim(300,600)
            plt.ylim(1E-5,1.0)
            #axis.set_ylim(1E-5,1E3)
            plt.show()
            print 'Saving'
