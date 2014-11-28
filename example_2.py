
import molecule
import scipy.optimize
import copy
import numpy as np
import known_molecules as km
from gas import Gas
from control_data import MeOH_control_set
import matplotlib.pyplot as plt


R = 8.314462175
pi = 3.141592653589793
e = 2.718281828459045

def demo_of_methanol_equlibrium():
        fig = plt.figure()
        fig.subplots_adjust(bottom=0.2, left=0.2, right = 0.8) # Make room for x-label
        ratio = 0.61803398             # Golden mean
        ratio = 0.4                     # This figure should be very wide to span two columns
        fig_width = 14.5
        fig_width = fig_width /2.54     # width in cm converted to inches
        fig_height = fig_width*ratio
        fig.set_size_inches(fig_width,fig_height)
        axis = fig.add_subplot(1,1,1)
        axis.semilogy()
        axis2 = axis.twinx()
         # in K
        Y = {}
        for mol in [km.H2, km.CO2, km.CO, km.H2O, km.CH3OH]:
            Y[mol] = []
        Pressure = 2.5
        reaction1 = {km.CO2: 0.0, km.H2: -2.0, km.CO: -1.0, km.CH3OH: 1.0, km.H2O: 0.0}
        reaction2 = {km.CO2: -1.0, km.H2: -1.0, km.CO: 1.0, km.CH3OH: 0.0, km.H2O: 1.0}
        amount_of_CO = 0.0
        amount_of_CO2 = 0.25
        amount_of_H2 = float(1.0-float(amount_of_CO+amount_of_CO2))
        Temperature = 400.
        gas_0 = Gas({km.H2: amount_of_H2, 
                     km.CO: amount_of_CO,
                     km.CO2: amount_of_CO2,
                     km.CH3OH: 0.0,
                     km.H2O: 0.0},temperature=Temperature)
        gas_0.set_pressure(Pressure)
        eq_result_0 = gas_0.gas_equlibrium_v2([reaction1,reaction2],guess=[0.001,0.001],T=Temperature)
        guess_last = eq_result_0[1]
        guess_new = guess_last
        CS = MeOH_control_set(2.5,cH2=0.75,cCO2=0.25,cCO=0.0,cMeOH=0.0,cH2O=0.0)
        for mol in ['MeOH','CO2', 'H2', 'CO', 'H2O']:
            CS[mol] = np.array(CS[mol])/100.*Pressure
        CS['T'] = np.array(CS['T']) +273.15
        T = CS['T']
        for t_i in T:
            print('Temperature: ' + str(t_i))
            Temperature = t_i
            gas_0 = Gas({km.H2: amount_of_H2, 
                     km.CO: amount_of_CO,
                     km.CO2: amount_of_CO2,
                     km.CH3OH: 0.0,
                     km.H2O: 0.0},temperature=Temperature)
            gas_0.set_pressure(Pressure)
            eq_result_0 = gas_0.gas_equlibrium_v2([reaction1,reaction2], guess=guess_new, T=Temperature)
            guess_last = eq_result_0[1]
            guess_last = guess_last + 0.5*(guess_last - guess_last)
            gas_0.plot_guess_historic(filename = 'fig-guess/' + str(t_i))
            for mol in [km.H2, km.CO2, km.CO, km.H2O, km.CH3OH]:
                Y[mol] += [eq_result_0[0].partial_pressures[mol]]
        for mol in [km.H2, km.CO2, km.CO, km.H2O, km.CH3OH]:
            Y[mol] = np.array(Y[mol])
        print('Calculation are Done')
        #print Y[km.CH3OH]
        color_list_mol = {km.CH3OH:'m', km.H2O:'b', km.CO:'r', km.CO2:'g', km.H2:'c'}
        for mol in [km.H2, km.CO2, km.CO, km.H2O, km.CH3OH]:
            axis.plot(T,Y[mol],color_list_mol[mol]+'-')
        color_list = {'MeOH':'m', 'H2O':'b', 'CO':'r', 'CO2':'g', 'H2':'c'}
        for mol in ['MeOH','CO2', 'H2', 'CO', 'H2O']:
            axis.plot(CS['T'], CS[mol],color_list[mol]+'--')
        axis2.plot(CS['T'], (Y[km.CH3OH] - CS['MeOH'])/CS['MeOH'], 'k-')
        """axis.annotate('6nm 20% coverage', xy=(15, 4calculate_partial_pressuresE-4),  xycoords='data',
                xytext=(10, 0.01), textcoords='data',
                horizontalalignment='left', verticalalignment='top',
                arrowprops=dict(arrowstyle='->',facecolor=color_list['MR260'],edgecolor=color_list['MR260']),fontsize=8,color=color_list['MR260']
                )"""
        #axis.set_ylim(1e-13,1E-6)
        axis2.set_ylim(-0.1,0.1)
        plt.xlim(300,600)
        #axis.legend(loc='upper right',prop={'size':6})
        axis.grid(False) 
        axis.tick_params(direction='in', length=6, width=1, colors='k',labelsize=8,axis='both',pad=3)
        axis2.tick_params(direction='in', length=6, width=1, colors='k',labelsize=8,axis='both',pad=3)
        axis.set_xlabel('Temperature / [C]', fontsize=8)
        axis.set_ylabel('Pressure / [bar]', fontsize=8)
        axis2.set_ylabel('Error in fraction', fontsize=8)
        #plt.tight_layout()
        #plt.show()
        print('Saving')
        plt.savefig('demo_2 test.png',dpi=600)
        plt.clf()
        plt.close()
        #del fig
        pass
        

if __name__ == '__main__':
    print 'Start gas.py'
    demo_of_methanol_equlibrium()
    
    
    
