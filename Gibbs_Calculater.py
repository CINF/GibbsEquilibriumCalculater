import math
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

H = {} # unit J/mol
H['CO(g)'] = -110500.0
H['CO2(g)'] = -393500.0
H['H2O(g)'] = -241830.0
H['H2(g)'] = 0.0
H['O2(g)'] = 0.0
H['CH4(g)'] = -74870.0
H['CH3OH(g)'] = -201300.0
H['N2(g)'] = 0.0
H['NH3(g)'] = -45940.0

H['CO(Ru-a)'] = -110500.0-77200.0 #(0.8eV)

S = {} # unit J/mol/K
S['CO(g)'] = 197.7
S['CO2(g)'] = 213.7
S['H2O(g)'] = 188.84
S['H2(g)'] = 130.679
S['CH4(g)'] = 186.25
S['CH3OH(g)'] = 239.9
S['N2(g)'] = 191.61
S['NH3(g)'] = 192.778

S['CO(Ru-a)'] = 0.0

def K_equilibrium(dG_reaction,T):
    R=8.3144621
    K = math.exp(-dG_reaction/(R*T))
    return K

def dG(dH_reaction,dS_reaction,T):
    dG_reaction = dH_reaction-T*dS_reaction
    return dG_reaction

print K_equilibrium(-1,1)


dH = {}
dH['MeOH'] = (H['CH3OH(g)']+H['H2O(g)'])-(H['CO2(g)']+3*H['H2(g)'])
dH['CO'] = (H['CO(g)']+H['H2O(g)'])-(H['CO2(g)']+H['H2(g)'])
dH['NH3'] = (H['NH3(g)'])-(0.5*H['N2(g)']+1.5*H['H2(g)'])

dH['CO(Ru-a)'] = (H['CO(Ru-a)'])-(H['CO(g)'])

dS = {}
dS['MeOH'] = (S['CH3OH(g)']+S['H2O(g)'])-(S['CO2(g)']+3*S['H2(g)'])
dS['CO'] = (S['CO(g)']+S['H2O(g)'])-(S['CO2(g)']+S['H2(g)'])
dS['NH3'] = (S['NH3(g)'])-(0.5*S['N2(g)']+1.5*S['H2(g)'])

dS['CO(Ru-a)'] = (S['CO(Ru-a)'])-(S['CO(g)'])

print 'dG MeOH: ' + str(dG(dH['MeOH'],dS['MeOH'],500))

p_factor = {}
p_factor['MeOH'] = -2.0
p_factor['CO'] = 0.0
p_factor['NH3'] = -1.0
p_factor['CO(Ru-a)'] = -1.0

#pressure = 1.0 #unit bar

Flow = {}
Flow['H2']=3.0
Flow['CO2']=1.0

def Equilibrium_full((MeOH, CO, H2O, CO2, H2),T,pressure):
    K_MeOH = K_equilibrium(dG(dH['MeOH'],dS['MeOH'],T),T)
    K_CO = K_equilibrium(dG(dH['CO'],dS['CO'],T),T)
    return (K_MeOH-((MeOH) * (H2O)) / ((CO2) * (H2)**3)*(pressure)**p_factor['MeOH'],
            K_CO-((CO) * (H2O)) / ((CO2) * (H2))*(pressure)**p_factor['CO'],
            (2*Flow['H2'])/(Flow['CO2'])-(4*MeOH+2*H2O+2*H2)/(MeOH+CO2+CO), #H/C
            (2*Flow['H2'])/(2*Flow['CO2'])-(4*MeOH+2*H2O+2*H2)/(MeOH+H2O+2*CO2+CO), #H/O
            1-(MeOH+CO+H2O+CO2+H2)
            )

def Equilibrium_simple((MeOH, H2O, CO2, H2),T,pressure):
    K_MeOH = K_equilibrium(dG(dH['MeOH'],dS['MeOH'],T),T)
    return (K_MeOH-(((MeOH) * (H2O)) / ((CO2) * (H2)**3))*(pressure)**p_factor['MeOH'],
            (2*Flow['H2'])/(Flow['CO2'])-(4*MeOH+2*H2O+2*H2)/(MeOH+CO2), #H/C
            (2*Flow['H2'])/(2*Flow['CO2'])-(4*MeOH+2*H2O+2*H2)/(MeOH+H2O+2*CO2), #H/O
            1-(MeOH+H2O+CO2+H2)
            )

def Equilibrium_NH3((NH3, N2, H2),T,pressure):
    K_NH3 = K_equilibrium(dG(dH['NH3'],dS['NH3'],T),T)
    return (K_NH3-(((NH3)) / ((N2)**0.5 * (H2)**1.5))*(pressure)**p_factor['NH3'],
            3.0/1.0-(2*H2+3*NH3)/(2*N2+1*NH3),
            1.0 - (NH3+N2+H2)
            )
#pressure = 1.0
def Equilibrium_RuCO((CO, COa,Ru),T,pressure):
    K_RuCO = K_equilibrium(dG(dH['CO(Ru-a)'],dS['CO(Ru-a)'],T),T)
    return (K_RuCO-(((COa)) / ((CO)*(Ru)))*(pressure)**p_factor['CO(Ru-a)'],
            1.0-(COa+Ru),
            pressure-CO,
            )

print Equilibrium_RuCO((1.0, 0.1,0.9),300,1.0)
print 'RuCO'
print fsolve(Equilibrium_RuCO, (1.0,0.3,0.1 ), args=(300.0, 1.0), xtol=1.49012e-12,maxfev=10000 )
print 'hej'
print fsolve(Equilibrium_RuCO, (1.0,0.3,0.7 ), args=(400.0,1.0), xtol=1.49012e-12,maxfev=10000 )
print fsolve(Equilibrium_RuCO, (1.0,0.1,0.9 ), args=(500.0,1.0), xtol=1.49012e-12,maxfev=10000 )
#print dH['NH3']
#print dS['NH3']
temp=298.15
print 'dG NH3 ved 500K: ' + str(dG(dH['NH3'],dS['NH3'],500))
print 'END --------'
print 'K: ' + str(K_equilibrium(dG(dH['NH3'],dS['NH3'],temp),temp))
cNH3, cN2,cH2 =  fsolve(Equilibrium_NH3, (1.0,0.001 ,0.001),args=(temp,1.0),xtol=1.49012e-12,maxfev=10000 ) #500K : NH3=0.0877
print 'NH3 ved '+str(temp)+'K: '+str(cNH3)
print 'cN2 ved 500K: '+str(cN2)
print 'cH2 ved 500K: '+str(cH2)
print 'p ved 500K: '+str(cNH3+ cN2+cH2)
print fsolve(Equilibrium_NH3, (cNH3,cN2 ,cH2),args=(temp,1.0),xtol=1.49012e-10,maxfev=10000 )
print Equilibrium_NH3((cNH3,cN2 ,cH2),temp,1.0)




if False:
    x = np.zeros(300)
    y1 = np.zeros(300)
    y2 = np.zeros(300)
    i = 0
    cMeOH, cCO, cH2O,cCO2,cH2 =  fsolve(Equilibrium_full, (0.01, 0.01,0.02,0.23,0.71),args=(500,1.0))

    cMeOH, cCO, cH2O,cCO2,cH2 =  fsolve(Equilibrium_full, (0.01, 0.01,0.02,0.23,0.71),args=(300,1.0))
    
    for t in range(300,600):
        print t
        cMeOH, cCO, cH2O,cCO2,cH2 =  fsolve(Equilibrium_full, (cMeOH, cCO, cH2O,cCO2,cH2),args=(t*1.0,1.0))
        x[i] = t
        y1[i] = cMeOH
        y2[i] = cCO
        i+=1
    if True:
        fig = plt.figure()
        fig.subplots_adjust(bottom=0.2, left=0.2) # Make room for x-label

        ratio = 0.61803398             # Golden mean
        #ratio = 0.4                     # This figure should be very wide to span two columns
        fig_width = 17
        fig_width = fig_width /2.54     # width in cm converted to inches
        fig_height = fig_width*ratio
        fig.set_size_inches(fig_width,fig_height)

        axis = fig.add_subplot(1,1,1)
        axis.semilogy()
        axis.plot(x,y1, 'r',label='MeOH')
        axis.plot(x,y2, 'g',label='CO')
        axis.legend(loc='upper right',prop={'size':10})
        plt.xlim(300,600)
        #axis.set_ylim(1E-5,1E3)
        axis.tick_params(direction='in', length=6, width=1, colors='k',labelsize=10,axis='both',pad=3)
        axis.set_xlabel('T / [K]', fontsize=10)
        axis.set_ylabel('c', fontsize=10)
        #plt.tight_layout()
        plt.show()
        print 'Saving'
        #plt.savefig('fig/Arrhenius_all.png',dpi=600)
        #plt.close()
        #plt.clf()
        #del fig

if False:
    x = np.zeros(300)
    y1 = np.zeros(300)
    y2 = np.zeros(300)
    i = 0
    cNH3, cN2,cH2 =  fsolve(Equilibrium_NH3, (0.01, 0.24,0.74),args=(500.0,1.0),xtol=1.49012e-10,maxfev=10000 )
    print 'NH3 ved 500K: '+str(cNH3)
    cNH3, cN2,cH2 =  fsolve(Equilibrium_NH3, (0.99, 0.01,0.01),args=(300.0,1.0),xtol=1.49012e-10,maxfev=10000)
    for t in range(300,600):
        #print t
        cNH3, cN2,cH2 =  fsolve(Equilibrium_NH3, (cNH3, cN2,cH2),args=(t*1.0,1.0),xtol=1.49012e-10,maxfev=10000)
        x[i] = t
        y1[i] = cNH3
        y2[i] = cN2
        i+=1
    if True:
        fig = plt.figure()
        fig.subplots_adjust(bottom=0.2, left=0.2) # Make room for x-label

        ratio = 0.61803398             # Golden mean
        #ratio = 0.4                     # This figure should be very wide to span two columns
        fig_width = 17
        fig_width = fig_width /2.54     # width in cm converted to inches
        fig_height = fig_width*ratio
        fig.set_size_inches(fig_width,fig_height)

        axis = fig.add_subplot(1,1,1)
        #axis.semilogy()
        axis.plot(x,y1, 'r',label='NH3')
        axis.plot(x,y2, 'g',label='N2')
        axis.legend(loc='upper right',prop={'size':10})
        plt.xlim(300,600)
        #axis.set_ylim(1E-5,1E3)
        axis.tick_params(direction='in', length=6, width=1, colors='k',labelsize=10,axis='both',pad=3)
        axis.set_xlabel('1000/T / [1/K]', fontsize=10)
        axis.set_ylabel('TOF / [molecules/site/s]', fontsize=10)
        #plt.tight_layout()
        plt.show()
        print 'Saving'
        #plt.savefig('fig/Arrhenius_all.png',dpi=600)
        #plt.close()
        #plt.clf()
        #del fig
if True:
    x = np.zeros(300)
    y1 = np.zeros(300)
    y2 = np.zeros(300)
    i = 0
    cNH3, cN2,cH2 =  fsolve(Equilibrium_NH3, (0.01, 0.24,0.74),args=(500.0,1.0),xtol=1.49012e-10,maxfev=10000 )
    print 'NH3 ved 500K: '+str(cNH3)
    cNH3, cN2,cH2 =  fsolve(Equilibrium_NH3, (0.6, 0.2,0.3),args=(600.0,1.0),xtol=1.49012e-10,maxfev=10000)
    fig = plt.figure()
    fig.subplots_adjust(bottom=0.2, left=0.2) # Make room for x-label

    ratio = 0.61803398             # Golden mean
    #ratio = 0.4                     # This figure should be very wide to span two columns
    fig_width = 17
    fig_width = fig_width /2.54     # width in cm converted to inches
    fig_height = fig_width*ratio
    fig.set_size_inches(fig_width,fig_height)

    axis = fig.add_subplot(1,1,1)
    for pi in [1.0,2.0,4.0]:
        x = np.zeros(300)
        y1 = np.zeros(300)
        y2 = np.zeros(300)
        i = 0
        for t in range(600,300,-1):
            #print t
            cNH3, cN2,cH2 =  fsolve(Equilibrium_NH3, (cNH3, cN2,cH2),args=(t*1.0,pi),xtol=1.49012e-10,maxfev=10000)
            x[i] = t
            y1[i] = cNH3
            y2[i] = cN2
            i+=1
        if True:

            #axis.semilogy()
            axis.plot(x,y1,label='p='+str(pi))
    axis.legend(loc='upper right',prop={'size':10})
    plt.xlim(300,600)
    #axis.set_ylim(1E-5,1E3)
    axis.tick_params(direction='in', length=6, width=1, colors='k',labelsize=10,axis='both',pad=3)
    axis.set_xlabel('1000/T / [1/K]', fontsize=10)
    axis.set_ylabel('TOF / [molecules/site/s]', fontsize=10)
    #plt.tight_layout()
    plt.show()
    print 'Saving'
    #plt.savefig('fig/Arrhenius_all.png',dpi=600)
    #plt.close()
    #plt.clf()
    #del fig
