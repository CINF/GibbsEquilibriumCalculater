import numpy as np
from random import shuffle
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math
import known_molecules as km

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


def matrix_problem(reduced_list_of_elements,element):
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
    if ((np.dot(A,C)-B)**2).sum() > 1e-10:
        return None
    else:
        return C

def compound_to_be_removed(reduced_list_of_elements,element,grade):
    C = matrix_problem(reduced_list_of_elements,element)
    #print C
    if C != None:
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
        if products == grade:
            element_to_be_removed = product_element
        elif reactans == grade:
            element_to_be_removed = reactans_element
        else:
            element_to_be_removed = None
        return element_to_be_removed
    else:
        return None
    

def simple_reaction(reduced_list_of_elements,element,list_of_elements_to_be_removed):
    C = matrix_problem(reduced_list_of_elements,element)
    if C == None:
        return None
    else:
        normalisationfactor = 1.0/(1.0*1.0)
        re_string = ''
        reactans = 0
        products = 0
        list_of_Cvalue = []
        composition_check=np.zeros(composition_size)
        re_string += str(-1.0/normalisationfactor)+'*'+str(element) + ' + '
        composition_check += -1.0/normalisationfactor * composition[element]
        list_of_Cvalue.append([element,-1.0/normalisationfactor])
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
            i += 1
        bonus_string = ''#'   BONUS: ' + str(products) +' '+ str(element) +' '+ str(element_to_be_removed)
        return_string = str(re_string[0:len(re_string)-3]) + ' = 0 ' + str(bonus_string)
        if products == len(list_of_elements_to_be_removed) or reactans == len(list_of_elements_to_be_removed):
            if sum(composition_check**2) > 0.001:
                print 'composition_check Error: ' + str(composition_check)
                print reduced_list_of_elements, element, list_of_elements_to_be_removed
                print 'Return string: ' + str(return_string)
            return return_string, list_of_Cvalue
        else:
            print 'Error: ' + str(products) +' ' + str(reactans) + ' ' + str(len(list_of_elements_to_be_removed))
            return None
        


def reaction_eq(list_of_elements,grade):
    reaction_simple_string = ''
    L_end = len(list_of_elements)
    for L in range(2,L_end+1):
        for rec in range(L_end**2):
            shuffle(list_of_elements)
            reaction_string = ''
            C_value_full = []
            for element in list_of_elements[0:L]:
                temp = list_of_elements[0:L]
                temp.remove(element)
                reduced_list_of_elements = temp
                list_of_elements_to_be_removed = compound_to_be_removed(reduced_list_of_elements,element,grade)
                if list_of_elements_to_be_removed != None:
                    reaction_string, C_value = simple_reaction(reduced_list_of_elements,element,list_of_elements_to_be_removed)
                    if reaction_string != None:
                        reaction_simple_string += str(reaction_string) + '\n'
                        C_value_full.append(C_value)
                        for i in range(len(list_of_elements_to_be_removed)):
                            list_of_elements.remove(list_of_elements_to_be_removed[i])
                        break
    if len(reaction_simple_string) != 0:
        return reaction_simple_string, list_of_elements
    else:
        return None, list_of_elements
    

def full_reaction(list_of_elements):
    return_string = ''
    grade = 1#1
    while grade <= 0.5*len(list_of_elements):
        for ii in range(10):
            #print 'Grade: ' + str(grade)
            part_return_string, list_of_elements = reaction_eq(list_of_elements,grade)
            #print list_of_elements
            #print part_return_string
            if part_return_string != None:
                return_string += part_return_string
                
                #grade +=1
        grade +=1
    return return_string, list_of_elements

for i in range(1):
    list_of_elements = ['H2','N2','NH3','H2O','O2','H2O2','CO','CO2','CH4','CH3OH','NO','NO2']
    #list_of_elements = ['H2','N2','O2','CO','CH4']
    total_reaction_list, list_of_elements = full_reaction(sorted(list_of_elements))
    print total_reaction_list
    print len(total_reaction_list)
    print sorted(list_of_elements)
    print len(list_of_elements)
