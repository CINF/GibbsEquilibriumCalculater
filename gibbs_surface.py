# -*- coding: utf-8 -*-
"""
Created on Wed Apr 02 09:45:43 2014

@author: aufn
"""

#import molecule
#import atom
import known_molecules as km
import gas
import numpy as np
import matplotlib.pyplot  as plt
import time
import scipy.interpolate as interpolate

#gas_0 = gas.Gas({km.CO: 0.01, km.H2: 0.99, km.CH4: 0.0, km.H2O: 0.0},temperature=310)
gas_0 = gas.Gas({km.H2: 0.75, 
                 km.CO: 0.0, 
                 km.CO2: 0.25, 
                 km.CH3OH: 0.0, 
                 km.H2O: 0.0},temperature=310)
#gas_1 = gas.Gas({km.H2: 0.5, km.CO: 0.0, km.CH3OH: 0.5, km.H2O: 0.0},temperature=310)
#gas_1 = gas.Gas({km.CO: 0.0, km.H2: 0.98, km.CH4: 0.01, km.H2O: 0.01},temperature=310)

P_final = 1.0
#P_final = 0.1
gas_0.set_pressure(P_final)

#reaction = {km.CO: -1.0, km.H2: -2.0, km.CH3OH: 1.0, km.H2O: 0.0}
reaction1 = {km.CO2: -1.0, km.H2: -3.0, km.CH3OH: 1.0, km.H2O: 1.0}
reaction2 = {km.CO2: -1.0, km.H2: -1.0, km.CO: 1.0, km.H2O: 1.0}
#reaction = {km.CO: -1.0, km.H2: -3.0, km.CH4: 1.0, km.H2O: 1.0}
#print gas.Gas(reaction).gas_atom_composition()
#{km.H2: -1.0, km.CO: -1.0, km.MeOH: 0.0, km.H2O: 0.0}
def Methanol(Temperature=None, Pressure=None,NO=None, gas_start=None):
    if Temperature == None:
        Temperature = 150.0+273.15 # K
    if Pressure == None:
        Pressure = 1.0 # bar
    if NO == None:
        NO = 5 # bar
    if gas_start == None:
        gas_start = {km.H2: amount_of_H2, 
                                 km.CO: amount_of_CO, 
                                 km.CO2: amount_of_CO2, 
                                 km.CH3OH: 0.0, 
                                 km.H2O: 0.0}
    if True:
        range_CO = [0.0,0.1]
        range_CO2 = [0.0,0.1]
        x = np.linspace(range_CO[0], range_CO[1], NO) # CO/H2
        y = np.linspace(range_CO2[0], range_CO2[1], NO) # CO2/H2
        X, Y = np.meshgrid(x, y)
        MeOH = np.zeros((len(x),len(y)))
        reaction0 = {km.CO2: 0.0, km.H2: -2.0, km.CO: -1.0, km.CH3OH: 1.0, km.H2O: 0.0}
        reaction1 = {km.CO2: -1.0, km.H2: -1.0, km.CO: 1.0, km.CH3OH: 0.0, km.H2O: 1.0}
        
        i=0
        i_max = len(x)*len(y)
        time_start = time.time()
        time_last = time_start
        max_point_0 = {'MeOH':0.0,'xi':[0,0],'XY':[0.0,0.0]}
        for xi in range(len(x)):
            for yi in range(len(y)):
                #print 'cordinates: ' + str(xi) + ', ' + str(yi)
                X[xi,yi] = x[xi]
                Y[xi,yi] = y[yi]
                gas_0 = gas.Gas(gas_start,temperature=Temperature)
                gas_0.set_pressure(Pressure)
                gas_test = gas_0 + X[xi,yi] *gas.Gas(reaction0,temperature=Temperature) + Y[xi,yi] *gas.Gas(reaction1,temperature=Temperature) 
                if (np.array(gas_test.partial_pressures.values()) >= 0.0).all():
                    gas_test.set_pressure(Pressure)
                    MeOH[xi,yi]=gas_test.gibbs()
                else:
                    MeOH[xi,yi]=0.0
                if time.time() - time_last > 10.0:
                    time_last = time.time()
                    print 'Procent done : ' + str(float(i)/float(i_max)*100)+ ' time: ' + str(time_last-time_start) 
                i += 1
        fig = plt.figure()
        fig.subplots_adjust(bottom=0.2, left=0.2, right=0.8,top=0.8)
        fig.set_size_inches(14/2.54,14/2.54)
        axis = fig.add_subplot(1,1,1)
        CS = axis.contour(X, Y, MeOH)
        
        axis.plot(X.flatten()[MeOH.argmin()],Y.flatten()[MeOH.argmin()],'r*')
        axis.text(X.flatten()[MeOH.argmin()],Y.flatten()[MeOH.argmin()],' x=%.3f, y=%.3f, g=%.3f' %(X.flatten()[MeOH.argmin()],Y.flatten()[MeOH.argmin()],MeOH.flatten()[MeOH.argmin()]))
        #CS = axis.contour(xnew, ynew, znew*100/Pressure,[0.125,0.25,0.5,1.0,2.0,4.0,8.0,16.0])
        #axis.plot(0.04*(1.0+np.linspace(0,1,100)),np.linspace(0,1,100),'-r')
        #axis.plot(np.linspace(0,1,100),0.1*(1.0+np.linspace(0,1,100)),'-b')
        #axis.plot(0.04/0.86,0.1/0.86,'*k')
        #axis.plot([0,1],[1,0],'-k')
        #axis.plot([0,0.9],[0.9,0],'-r')
        #axis.plot([0.7,0.7+(x[1]-x[0]),0.7+(x[1]-x[0])],[0.7,0.7,0.7+(y[1]-y[0])],'-r')
        #axis.plot(max_point_0['XY'][0],max_point_0['XY'][1],'*r')
        axis.set_xlabel('c(CO)', fontsize=8)
        axis.set_ylabel('c(CO2)', fontsize=8)
        axis.text(0.65  ,0.6,'T=%.1fK' %Temperature, fontsize=8)
        axis.text(0.65,0.5,'T=%.1f$^\circ$C' %(Temperature-273.15), fontsize=8)
        axis.text(0.65,0.4,'P=%.1fBar' %Pressure, fontsize=8)
        axis.set_xlim(0*range_CO[0], range_CO[1])
        axis.set_ylim(0*range_CO2[0], range_CO2[1])
        plt.clabel(CS,inline=1, fontsize=8)
        plt.savefig('fig/Gibbs_surface_%.1fK_%.1fbar_CO-%.3f_CO2-%.3f.'%(Temperature,Pressure,gas_start[km.CO],gas_start[km.CO2])+file_format,dpi=600)
        fig.clf()
        plt.close()
        del fig

if __name__ == '__main__':
    file_format = 'png'
    Pressure = 1.0
    Temperature = 150.0+273.15
    for x in np.linspace(0.0,0.5,20):
        print x
        amount_of_CO = x
        amount_of_CO2 = x
        amount_of_H2 = 1.0 - amount_of_CO - amount_of_CO2
        gas_start = {km.H2: amount_of_H2, 
                     km.CO: amount_of_CO, 
                     km.CO2: amount_of_CO2, 
                     km.CH3OH: 0.0, 
                     km.H2O: 0.0}
        #gas_start.set_pressure(Pressure)
        Methanol(Temperature=Temperature, Pressure=Pressure,NO=50,gas_start=gas_start)