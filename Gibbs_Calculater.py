import numpy as np
from random import shuffle
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math


###
# list of functions:
#def K_equilibrium(dG_reaction,T):
#    return K
#def dG(dH_reaction,dS_reaction,T):
#    return dG_reaction
#def dGv2(C_value, T):
#    return dG_reaction
#def Equilibrium_full((MeOH, CO, H2O, CO2, H2),T,pressure):
#    return [X,X,X,X]
#def Equilibrium_simple((MeOH, H2O, CO2, H2),T,pressure):
#    return [X,X,X,X]
#def Equilibrium_NH3((NH3, N2, H2),T,pressure):
#    return [X,X,X]
#def matrix_problem(reduced_list_of_elements,element):
#    return C_value
#def compound_to_be_removed(reduced_list_of_elements,element,grade):
#    return elements_to_be_removed
#def simple_reaction(reduced_list_of_elements,element,list_of_elements_to_be_removed):
#    return return_string, list_of_Cvalue
#def reaction_eq(list_of_elements,grade):
#    return reaction_simple_string, list_of_elements
#def full_reaction(list_of_elements):
#    return return_string, list_of_elements


#Enthalpy
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

#Entropy
S = {} # unit J/mol/K
S['CO(g)'] = 197.7
S['CO2(g)'] = 213.7
S['H2O(g)'] = 188.84
S['H2(g)'] = 130.679
S['CH4(g)'] = 186.25
S['CH3OH(g)'] = 239.9
S['N2(g)'] = 191.61
S['NH3(g)'] = 192.778


composition = {}
composition_size = 10
composition['H']=np.zeros(composition_size)
composition['H'][1]=1

composition['He']=np.zeros(composition_size)
composition['He'][2]=1

composition['Li']=np.zeros(composition_size)
composition['Li'][3]=1

composition['Be']=np.zeros(composition_size)
composition['Be'][4]=1

composition['B']=np.zeros(composition_size)
composition['B'][5]=1

composition['C']=np.zeros(composition_size)
composition['C'][6]=1

composition['N']=np.zeros(composition_size)
composition['N'][7]=1

composition['O']=np.zeros(composition_size)
composition['O'][8]=1

composition['H2O']=np.zeros(composition_size)
composition['H2O'][1]=2
composition['H2O'][8]=1

composition['H2O2']=np.zeros(composition_size)
composition['H2O2'][1]=2
composition['H2O2'][8]=2

composition['NH3']=np.zeros(composition_size)
composition['NH3'][1]=3
composition['NH3'][7]=1

composition['CO']=np.zeros(composition_size)
composition['CO'][6]=1
composition['CO'][8]=1

composition['CO2']=np.zeros(composition_size)
composition['CO2'][6]=1
composition['CO2'][8]=2

composition['CH4']=np.zeros(composition_size)
composition['CH4'][1]=4
composition['CH4'][6]=1

composition['CH3OH']=np.zeros(composition_size)
composition['CH3OH'][1]=4
composition['CH3OH'][6]=1
composition['CH3OH'][8]=1

composition['H2']=np.zeros(composition_size)
composition['H2'][1]=2

composition['N2']=np.zeros(composition_size)
composition['N2'][7]=2

composition['O2']=np.zeros(composition_size)
composition['O2'][8]=2

composition['NO']=np.zeros(composition_size)
composition['NO'][7]=1
composition['NO'][8]=1

composition['NO2']=np.zeros(composition_size)
composition['NO2'][7]=1
composition['NO2'][8]=2



#Reaction equlibrium constant
def K_equilibrium(dG_reaction,T):
    R=8.3144621 # Gas constant
    K = math.exp(-dG_reaction/(R*T))
    return K

#Change in Gibbs free energy
def dG(dH_reaction,dS_reaction,T):
    dG_reaction = dH_reaction-T*dS_reaction
    return dG_reaction

"""
#Change in Gibbs free energy
def dGv2(C_value, T):
    dG_reaction = 0
    for i in range(len(C_value)):
        dG_reaction += C_value[i][1] * ( H(C_value[i][0]) - T * S(C_value[i][0]))
    return dG_reaction
"""

def Equilibrium_full((MeOH, CO, H2O, CO2, H2),T,pressure):
    K_MeOH = K_equilibrium(dG(dH['MeOH'],dS['MeOH'],T),T)
    K_CO = K_equilibrium(dG(dH['CO'],dS['CO'],T),T)
    return (K_MeOH-((MeOH) * (H2O)) / ((CO2) * (H2)**3)*(pressure)**p_factor['MeOH'],
            K_CO-((CO) * (H2O)) / ((CO2) * (H2))*(pressure)**p_factor['CO'],
            (2*Flow['H2'])/(Flow['CO2'])-(4*MeOH+2*H2O+2*H2)/(MeOH+CO2+CO), #H/C
            (2*Flow['H2'])/(2*Flow['CO2'])-(4*MeOH+2*H2O+2*H2)/(MeOH+H2O+2*CO2+CO), #H/O
            1-(MeOH+CO+H2O+CO2+H2)
            )

#CO2 +3H2 -> CH3OH + H2O
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
def Equilibrium_CH4((CH4, H2O, CO, H2),T,pressure):
    K_CH4 = K_equilibrium(dG(dH['CH4'],dS['CH4'],T),T)
    return (K_CH4-(((CH4)*(H2O)) / ((CO) * (H2)**3))*(pressure)**p_factor['CH4'],
            3*2.0/1.0-(2*H2+2*H2O+4*CH4)/(CO+CH4), #H/C
            3*2.0/1.0-(2*H2+2*H2O+4*CH4)/(CO+H2O), #H/O
            1.0 - (CH4 + H2O + CO + H2)
            )
def Equilibrium_Ru_CO((cCO),T,pressure):
    K_Ru_CO = K_equilibrium(dG(dH['Ru_CO'],dS['Ru_CO'],T),T)
    return (K_Ru_CO-(cCO)*(pressure)**p_factor['Ru_CO']
            )

def matrix_problem(reduced_list_of_elements,element):
    #print 'matrix_problem'
    #print 'Long list: ' + str(reduced_list_of_elements)
    #print 'Element: ' + str(element)
    A = np.zeros((composition_size,composition_size))
    i = 0
    for reduced_element in reduced_list_of_elements:
        A[:,i]=composition[reduced_element]
        i += 1
    B=np.zeros((composition_size,1))
    B[:,0] = composition[element]
    x = np.linalg.lstsq(A,B,rcond=0.0001)
    C=np.zeros((composition_size,1))
    C[:,0] = [ elem for elem in x[0] ]#[ round(elem, 9) for elem in x[0] ]
    #print 'Sum: ' + str(((np.dot(A,C)-B)**2).sum())
    if ((np.dot(A,C)-B)**2).sum() > 1e-10:
        #print 'No solution found: ' + str(((np.dot(A,C)-B)**2).sum())
        #print str(reduced_list_of_elements) + str(element)
        #print A[:,0:4]
        #print C[0:4,:]
        #print B
        return 'NaN'
    else:
        return C

def compound_to_be_removed(reduced_list_of_elements,element,grade):
    #print 'compound_to_be_removed'
    C = matrix_problem(reduced_list_of_elements,element)
    #print C
    if C != 'NaN':
        product_element = []
        reactans_element = []
        products = 1
        product_element += [element]
        reactans = 0
        i = 0
        for reduced_element in reduced_list_of_elements:
            if C[i][0] < 0.001 and C[i][0]**2>0.001:
                products += 1
                product_element += [reduced_element]
            if C[i][0] > 0.001 and C[i][0]**2>0.001:
                reactans +=1
                reactans_element += [reduced_element]
            i+=1
        #print product_element, reactans_element
        #print products, reactans
        if products == grade:
            element_to_be_removed = product_element
            #return (element_to_be_removed)
        elif reactans == grade:
            element_to_be_removed = reactans_element
            #return (element_to_be_removed)
        else:
            element_to_be_removed = 'NaN'
        #print element_to_be_removed
        return element_to_be_removed
    else:
        return 'NaN'

def simple_reaction(reduced_list_of_elements,element,list_of_elements_to_be_removed):
    #print 'simple_reaction'
    #print reduced_list_of_elements, element, list_of_elements_to_be_removed
    C = matrix_problem(reduced_list_of_elements,element)
    #print C
    if C == 'NaN':
        return 'NaN'
    else:
        normalisationfactor = 1.0/(1.0*1.0)
##        if element in list_of_elements_to_be_removed:
##            normalisationfactor = 1.0
##        else:
##            i = 0
##            for reduced_element in reduced_list_of_elements:
##                if reduced_element in list_of_elements_to_be_removed[0]:
##                    if C[i][0]**2>0.001:
##                        normalisationfactor = C[i][0]
##                        #normalisationfactor = abs(C[i][0])
##                        #print normalisationfactor
##                        break
##                    else:
##                        return 'NaN'
##                i +=1
        #print 'normalisationfactor: ' + str(normalisationfactor)
        #print 'products: ' + str(list_of_elements_to_be_removed) + ' Element: ' + str(element)
        re_string = ''
        reactans = 0
        products = 0
        list_of_Cvalue = []
        composition_check=np.zeros(composition_size)
        #print 'Element: ' + str(element)
        #print composition_check
        re_string += str(-1.0/normalisationfactor)+'*'+str(element) + ' + '
        composition_check += -1.0/normalisationfactor * composition[element]
        list_of_Cvalue.append([element,-1.0/normalisationfactor])
        #print composition_check
        if -1.0/normalisationfactor < 0:
            products += 1
        elif -1.0/normalisationfactor < 0:
            reactans +=1
        i = 0
        for reduced_element in reduced_list_of_elements:
            C_value = C[i][0]/normalisationfactor
            if C_value**2>0.001:
                re_string += str(C_value)+'*'+str(reduced_element) + ' + '
                composition_check += (C_value) * composition[reduced_element]
                list_of_Cvalue.append([reduced_element,C_value])
                if C_value < 0.0:
                    products += 1
                elif C_value > 0.0:
                    reactans +=1
                #print reduced_element
                #print composition_check
            i += 1

        bonus_string = ''#'   BONUS: ' + str(products) +' '+ str(element) +' '+ str(element_to_be_removed)
        return_string = str(re_string[0:len(re_string)-3]) + ' = 0 ' + str(bonus_string)
        
        #print 'String result: ' + str(return_string)
        if products == len(list_of_elements_to_be_removed) or reactans == len(list_of_elements_to_be_removed):
            if sum(composition_check**2) > 0.001:
                print 'composition_check Error: ' + str(composition_check)
                print reduced_list_of_elements, element, list_of_elements_to_be_removed
                print 'Return string: ' + str(return_string)
            #print return_string
            return return_string, list_of_Cvalue
        else:
            print 'Error: ' + str(products) +' ' + str(reactans) + ' ' + str(len(list_of_elements_to_be_removed))
            return 'NaN'


#print list_of_elements
def reaction_eq(list_of_elements,grade):
    reaction_simple_string = ''
    #print 'start'
    #print 'Total list: ' + str(list_of_elements)
    L_end = len(list_of_elements)
    #print 'Length of list: ' + str(L_end)
    for L in range(2,L_end+1):
        for rec in range(len(list_of_elements)**2): #range(len(list_of_elements)**2):
            shuffle(list_of_elements)
            #print list_of_elements
            #print 'new# iteration'
            reaction_string = ''
            C_value_full = []
            for element in list_of_elements[0:L]:
                temp = list_of_elements[0:L]
                #print temp
                temp.remove(element)
                reduced_list_of_elements = temp
                #print 'Long list: ' + str(reduced_list_of_elements)
                #print 'Element: ' + str(element)
                list_of_elements_to_be_removed = compound_to_be_removed(reduced_list_of_elements,element,grade)
                #print 'list_of_elements_to_be_removed: ' + str(list_of_elements_to_be_removed)
                if list_of_elements_to_be_removed != 'NaN':
                    #print 'list_of_elements_to_be_removed: ' + str(list_of_elements_to_be_removed)
                    #for i in range(len(list_of_elements_to_be_removed)):
                        #print 'remove: '
                        #print list_of_elements_to_be_removed[i]                    
                    reaction_string, C_value = simple_reaction(reduced_list_of_elements,element,list_of_elements_to_be_removed)
                    #print 'reaction_string: ' + str(reaction_string)
                    if reaction_string != 'NaN':
                        reaction_simple_string += str(reaction_string) + '\n'
                        C_value_full.append(C_value)
                        for i in range(len(list_of_elements_to_be_removed)):
                            #print 'remove: '
                            #print list_of_elements_to_be_removed[i]
                            list_of_elements.remove(list_of_elements_to_be_removed[i])
                        #print list_of_elements
                        break

    #print 'remaining list: '+str(list_of_elements)
    #print 'iteration Done'
    #print reaction_simple_string         
    #print 'end'
    if len(reaction_simple_string) != 0:
        return reaction_simple_string, list_of_elements
    else:
        return 'NaN', list_of_elements

def full_reaction(list_of_elements):
    return_string = ''
    grade = 1#1
    while grade <= 0.5*len(list_of_elements):
        for ii in range(10):
            #print 'Grade: ' + str(grade)
            part_return_string, list_of_elements = reaction_eq(list_of_elements,grade)
            #print list_of_elements
            #print part_return_string
            if part_return_string != 'NaN':
                return_string += part_return_string
                
                #grade +=1
        grade +=1
    return return_string, list_of_elements

dH = {}
dH['MeOH'] = (H['CH3OH(g)']+H['H2O(g)'])-(H['CO2(g)']+3*H['H2(g)'])
dH['CO'] = (H['CO(g)']+H['H2O(g)'])-(H['CO2(g)']+H['H2(g)'])
dH['NH3'] = (H['NH3(g)'])-(0.5*H['N2(g)']+1.5*H['H2(g)'])
dH['CH4'] = (H['CH4(g)']+H['H2O(g)'])-(3*H['H2(g)']+H['CO(g)'])
dH['Ru_CO'] = (H['CO(g)']-67500)-(H['CO(g)'])#0.7ev/molecule

dS = {}
dS['MeOH'] = (S['CH3OH(g)']+S['H2O(g)'])-(S['CO2(g)']+3*S['H2(g)'])
dS['CO'] = (S['CO(g)']+S['H2O(g)'])-(S['CO2(g)']+S['H2(g)'])
dS['NH3'] = (S['NH3(g)'])-(0.5*S['N2(g)']+1.5*S['H2(g)'])
dS['CH4'] = (S['CH4(g)']+S['H2O(g)'])-(3*S['H2(g)']+S['CO(g)'])
dS['Ru_CO'] = (0.0)-(S['CO(g)'])

print 'dG MeOH: ' + str(dG(dH['MeOH'],dS['MeOH'],500))

p_factor = {}
p_factor['MeOH'] = -2.0
p_factor['CO'] = 0.0
p_factor['NH3'] = -1.0
p_factor['CH4'] = -2.0
p_factor['CO(Ru-a)'] = -1.0
p_factor['Ru_CO'] = -1.0

#pressure = 1.0 #unit bar

Flow = {}
Flow['H2']=3.0
Flow['CO2']=1.0
    

#list_of_elements = ['H2','H2O','CO','CO2','CH3OH']
for i in range(1):
    list_of_elements = ['H2','N2','NH3','H2O','O2','H2O2','CO','CO2','CH4','CH3OH','NO','NO2']
    #list_of_elements = ['H2','N2','O2','CO','CH4']
    total_reaction_list, list_of_elements = full_reaction(sorted(list_of_elements))
    print total_reaction_list
    print len(total_reaction_list)
    print sorted(list_of_elements)
    print len(list_of_elements)

#print Equilibrium_RuCO((1.0, 0.1,0.9),300,1.0)
#print 'RuCO'
#print fsolve(Equilibrium_RuCO, (1.0,0.3,0.1 ), args=(300.0, 1.0), xtol=1.49012e-12,maxfev=10000 )
#print 'hej'
#print fsolve(Equilibrium_RuCO, (1.0,0.3,0.7 ), args=(400.0,1.0), xtol=1.49012e-12,maxfev=10000 )
#print fsolve(Equilibrium_RuCO, (1.0,0.1,0.9 ), args=(500.0,1.0), xtol=1.49012e-12,maxfev=10000 )
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


fig = plt.figure()
axis = fig.add_subplot(1,1,1)
axis.tick_params(direction='in', length=6, width=1, colors='k',labelsize=10,axis='both',pad=3)
axis.set_xlabel('T / [K]', fontsize=10)
axis.set_ylabel('c', fontsize=10)


if True:
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
        axis.semilogy()
        axis.plot(x,y1, 'r',label='MeOH')
        axis.plot(x,y2, 'g',label='CO')
        axis.legend(loc='upper right',prop={'size':10})
        plt.xlim(300,600)
        #axis.set_ylim(1E-5,1E3)
        plt.show()
        print 'Saving'
if False:
    x = np.zeros(300)
    y1 = np.zeros(300)
    y2 = np.zeros(300)
    i = 0
    cCH4, cH2O, cCO, cH2 =  fsolve(Equilibrium_CH4, (0.01, 0.01,0.01,0.97),args=(500,1.0))

    cCH4, cH2O, cCO, cH2 =  fsolve(Equilibrium_CH4, (0.01, 0.01,0.01,0.97),args=(500,1.0))
    
    for t in range(300,600):
        print t
        cCH4, cH2O, cCO, cH2 =  fsolve(Equilibrium_CH4, (cCH4, cH2O, cCO, cH2),args=(t*1.0,1.0))
        x[i] = t
        y1[i] = cCH4
        y2[i] = cCO
        i+=1
    if True:
        axis.semilogy()
        axis.plot(x,y1, 'r',label='CH4')
        axis.plot(x,y2, 'g',label='CO')
        axis.legend(loc='upper right',prop={'size':10})
        plt.xlim(300,600)
        #axis.set_ylim(1E-5,1E3)
        plt.show()
        print 'Saving'


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
        axis.plot(x,y1, 'r',label='NH3')
        axis.plot(x,y2, 'g',label='N2')
        axis.legend(loc='upper right',prop={'size':10})
        plt.xlim(300,600)
        #axis.set_ylim(1E-5,1E3)
        plt.show()
        print 'Saving'

if False:
    x = np.zeros(300)
    y1 = np.zeros(300)
    y2 = np.zeros(300)
    i = 0
    cNH3, cN2,cH2 =  fsolve(Equilibrium_NH3, (0.01, 0.24,0.74),args=(500.0,1.0),xtol=1.49012e-10,maxfev=10000 )
    print 'NH3 ved 500K: '+str(cNH3)
    cNH3, cN2,cH2 =  fsolve(Equilibrium_NH3, (0.6, 0.2,0.3),args=(600.0,1.0),xtol=1.49012e-10,maxfev=10000)

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


