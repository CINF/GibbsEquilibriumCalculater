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

if __name__ == '__main__':
    file_format = 'png'
    T_max = 600
    T_min = 200
    T_steps = 50
    #for ti in np.linspace(T_min,T_max,10):
    #    print gas.Gas(reaction).gibbs(T=ti)
    x = np.linspace(T_min,T_max,T_steps)
    yCO,yH2,yCO2,yH2O,yCH3OH = [],[],[],[],[]

    result = []
    if False:
        for i in range(len(x)):
            print i
            result.append(gas_0.gas_equlibrium([reaction1,reaction2],T=x[i]))
            yH2.append(result[i].partial_pressures[km.H2])
            yCO.append(result[i].partial_pressures[km.CO])
            yCO2.append(result[i].partial_pressures[km.CO2])
            yH2O.append(result[i].partial_pressures[km.H2O])
            yCH3OH.append(result[i].partial_pressures[km.CH3OH])
            #y1[i] =result.gibbs(x[i])
            #y2[i] =gas_0.partial_pressures[km.CH4]
            #y2[i] =gas_0.gibbs(x[i])
            #y3[i] =gas_1.partial_pressures[km.CH4]
            #y3[i] =gas_1.gibbs(x[i])
        sums = np.array(yCO)+np.array(yH2)+np.array(yCO2)+np.array(yH2O)+np.array(yCH3OH)
    if False:
        fig = plt.figure()
        fig.subplots_adjust(bottom=0.2, left=0.2, right=0.8,top=0.95) # Make room for x-label
        #ratio = 0.61803398             # Golden mean
        ratio = 0.4                    # This figure should be very wide to span two columns
        fig_width = 11
        fig_width = fig_width /2.54     # width in cm converted to inches
        fig_height = fig_width*ratio
        fig.set_size_inches(fig_width,fig_height)
        axis = fig.add_subplot(1,1,1)
        color_list=['b','g','r','k','b','k']
        axis.plot(x,np.array(yCH3OH)/sums*100, 'r-',linewidth=1.0)
        axis.plot(x,np.array(yCO)/sums*100, 'g-',linewidth=1.0)
        axis.plot(x,np.array(yH2O)/sums*100, 'b-',linewidth=1.0)
        #axis.plot(x,np.array(sums)/sums*100, 'c-',linewidth=1.0)
        #axis.plot(x,y2/sums, 'r-',linewidth=1.0)
        #axis.plot(x,y3/sums, 'c-')
        #axis.plot(x,y4/sums, 'k--')
        #axis.plot(x,(y1+y2+y3+y4)/sums, 'g-',linewidth=1.0)
        axis.set_xlim(250,600)
        #axis.set_ylim(-0.1,1.1)
        #axis.set_yscale('log')
        axis.tick_params(direction='in', length=6, width=1, colors='k',labelsize=8,axis='both',pad=3)
        axis.set_xlabel('Temperature / [K]', fontsize=8)
        axis.set_ylabel('CH4 / [au]', fontsize=8)
        plt.savefig('test.'+file_format,dpi=600)
        fig.clf()
        plt.close()
        del fig
    if False:
        y_final=[]
        y_ini = []
        for i in range(len(result)):
            y_final.append(result[i].gibbs(T=x[i])) 
            y_ini.append(gas_0.gibbs(T=x[i])) 
        fig = plt.figure()
        fig.subplots_adjust(bottom=0.2, left=0.2, right=0.8,top=0.95) # Make room for x-label
        #ratio = 0.61803398             # Golden mean
        ratio = 0.4                    # This figure should be very wide to span two columns
        fig_width = 11
        fig_width = fig_width /2.54     # width in cm converted to inches
        fig_height = fig_width*ratio
        fig.set_size_inches(fig_width,fig_height)
        axis = fig.add_subplot(1,1,1)
        color_list=['b','g','r','k','b','k']
        axis.plot(x,y_final, 'k-',linewidth=1.0)
        axis.plot(x,y_ini, 'g-',linewidth=1.0)
        #axis.set_xlim(250,1000)
        #axis.set_ylim(-0.1,1.2)
        #axis.set_yscale('log')
        axis.tick_params(direction='in', length=6, width=1, colors='k',labelsize=8,axis='both',pad=3)
        axis.set_xlabel('Temperature / [K]', fontsize=8)
        axis.set_ylabel('CH4 / [au]', fontsize=8)
        plt.savefig('test_gibbs.'+file_format,dpi=600)
        fig.clf()
        plt.close()
        del fig
    if True:
        Pressure = 2.5 # bar
        Temperature = 150+273.15 # K
        range_CO = [0.0,1.0]
        range_CO2 = [0.0,1.0]
        NO = 10
        #delta=MAX/30.0
        x = np.linspace(range_CO[0], range_CO[1], NO) # CO/H2
        y = np.linspace(range_CO2[0], range_CO2[1], NO) # CO2/H2
        X, Y = np.meshgrid(x, y)
        result_0 = {}#np.zeros((len(x),len(y)))
        MeOH = np.zeros((len(x),len(y)))
        reaction1 = {km.CO2: -1.0, km.H2: -3.0, km.CH3OH: 1.0, km.H2O: 1.0}
        reaction2 = {km.CO2: -1.0, km.H2: -1.0, km.CO: 1.0, km.H2O: 1.0}
        i=0
        i_max = len(x)*len(y)
        time_start = time.time()
        time_last = time_start
        max_point_0 = {'MeOH':0.0,'xi':[0,0]}
        for xi in range(len(x)):
            for yi in range(len(y)):
                amount_of_H2 = float(1.0-float(X[xi,yi]+Y[xi,yi]))
                amount_of_CO = float(X[xi,yi])
                amount_of_CO2 = float(Y[xi,yi])
                gas_0 = gas.Gas({km.H2: amount_of_H2, 
                                 km.CO: amount_of_CO, 
                                 km.CO2: amount_of_CO2, 
                                 km.CH3OH: 0.0, 
                                 km.H2O: 0.0},temperature=Temperature)
                gas_0 = gas_0.set_pressure(Pressure)
                if amount_of_H2 < 0.0 or amount_of_CO < 0.0 or amount_of_CO2 < 0.0:
                    MeOH[xi,yi] = 0.0
                    result_0[X[xi,yi],Y[xi,yi]] = gas_0
                else:
                    result_0[X[xi,yi],Y[xi,yi]]=gas_0.gas_equlibrium([reaction1,reaction2],T=Temperature)
                if float(result_0[X[xi,yi],Y[xi,yi]].partial_pressures[km.CH3OH]) < 0.0:
                    print 'Error negative partial pressure'
                MeOH[xi,yi]=float(result_0[X[xi,yi],Y[xi,yi]].partial_pressures[km.CH3OH]) # in %
                if MeOH[xi,yi] > max_point_0['MeOH']:
                    max_point_0['MeOH'] = MeOH[xi,yi]
                    max_point_0['xi'] = [xi,yi]
                    max_point_0['gas'] = result_0[X[xi,yi],Y[xi,yi]]
                i +=1
                if time.time()-time_last > 10.0:
                    time_last = time.time()
                    print 'Procent done : ' + str(float(i)/float(i_max)*100)
        print max_point_0
        fig = plt.figure()
        fig.subplots_adjust(bottom=0.2, left=0.2, right=0.8,top=0.8)
        fig.set_size_inches(14/2.54,14/2.54)
        axis = fig.add_subplot(1,1,1)
        CS = axis.contour(X, Y, MeOH*100/Pressure,[0.125,0.25,0.5,1.0,2.0,4.0,8.0,16.0])
        #axis.plot(0.04*(1.0+np.linspace(0,1,100)),np.linspace(0,1,100),'-r')
        #axis.plot(np.linspace(0,1,100),0.1*(1.0+np.linspace(0,1,100)),'-b')
        #axis.plot(0.04/0.86,0.1/0.86,'*k')
        axis.set_xlabel('c(CO)', fontsize=8)
        axis.set_ylabel('c(CO2)', fontsize=8)
        axis.text(0.55,0.3,'T=%.1fK' %Temperature, fontsize=8)
        axis.text(0.55,0.2,'T=%.1f$^\circ$C' %(Temperature-273.15), fontsize=8)
        axis.text(0.55,0.1,'P=%.1fBar' %Pressure, fontsize=8)
        axis.set_xlim(range_CO[0], range_CO[1])
        axis.set_ylim(range_CO2[0], range_CO2[1])
        plt.clabel(CS,inline=1, fontsize=8)
        plt.savefig('CS.'+file_format,dpi=600)
        fig.clf()
        plt.close()
        del fig
    if False: # fig 6.5
        Pressure = 2.5 # bar
        Temperature = 150+273.15 # K
        range_CO2 = [0.0,1.0]
        NO = 20
        #delta=MAX/30.0
        x = np.linspace(range_CO2[0], range_CO2[1], NO)
        result_1 = {}#np.zeros((len(x),len(y)))
        yMeOH = np.zeros(len(x))
        yH2O = np.zeros(len(x))
        yCO = np.zeros(len(x))
        reaction1 = {km.CO2: -1.0, km.H2: -3.0, km.CH3OH: 1.0, km.H2O: 1.0}
        reaction2 = {km.CO2: -1.0, km.H2: -1.0, km.CO: 1.0, km.H2O: 1.0}
        i=0
        i_max = len(x)
        time_start = time.time()
        time_last = time_start
        max_point_1 = {'MeOH':0.0,'xi':0}
        amount_of_CO = 0.0
        for xi in range(len(x)):
            amount_of_H2 = 1.0-float(x[xi]+amount_of_CO)
            #amount_of_CO = float(X[xi,yi])
            amount_of_CO2 = float(x[xi])
            gas_0 = gas.Gas({km.H2: amount_of_H2, 
                             km.CO: amount_of_CO, 
                             km.CO2: amount_of_CO2, 
                             km.CH3OH: 0.0, 
                             km.H2O: 0.0},temperature=Temperature)
            gas_0.set_pressure(Pressure)
            if amount_of_H2 < 0.0 or amount_of_CO < 0.0 or amount_of_CO2 < 0.0:
                result_1[x[xi]]=gas_0
            else:
                result_1[x[xi]]=gas_0.gas_equlibrium([reaction1,reaction2],T=Temperature)
            if float(result_1[x[xi]].partial_pressures[km.CH3OH]) < 0.0:
                print 'Error negative partial pressure'
            yMeOH[xi]=float(result_1[x[xi]].partial_pressures[km.CH3OH]) # in %
            yH2O[xi]=float(result_1[x[xi]].partial_pressures[km.H2O]) # in %
            yCO[xi]=float(result_1[x[xi]].partial_pressures[km.CO]) # in %
            if yMeOH[xi] > max_point_1['MeOH']:
                max_point_1['MeOH'] = yMeOH[xi]
                max_point_1['xi'] = xi
                max_point_1['gas'] = result_1[x[xi]]
            i +=1
            if time.time()-time_last > 10.0:
                time_last = time.time()
                print 'Procent done : ' + str(float(i)/float(i_max)*100)
        print max_point_1
        fig = plt.figure()
        fig.subplots_adjust(bottom=0.2, left=0.2, right=0.8,top=0.8)
        fig.set_size_inches(14/2.54,14/2.54)
        axis = fig.add_subplot(1,1,1)
        axis.plot(x,yMeOH*100.0/Pressure,'g--')
        axis.plot(x,yH2O*100.0/Pressure,'b--')
        axis.plot(x,yCO*100.0/Pressure,'r--')
        axis.set_xlabel('c(CO2)', fontsize=8)
        axis.set_ylabel('c(MeOH) at $\%$', fontsize=8)
        axis.text(0.75,4,'T=%.1fK' %Temperature, fontsize=8)
        axis.text(0.75,3.5,'T=%.1f$^\circ$C' %(Temperature-273.15), fontsize=8)
        axis.text(0.75,3.0,'P=%.1fBar' %Pressure, fontsize=8)
        axis.set_xlim(range_CO2[0], range_CO2[1])
        axis.set_ylim(0,20)
        #plt.clabel(CS, inline=1, fontsize=8)
        plt.savefig('Figure 6.5.'+file_format,dpi=600)
        fig.clf()
        plt.close()
        del fig
