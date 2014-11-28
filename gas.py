
import molecule
import scipy.optimize
import copy
import numpy as np
import known_molecules as km
import time
import matplotlib.pyplot as plt
R = 8.314462175
pi = 3.141592653589793
e = 2.718281828459045


class Gas():
    def __init__(self,partial_pressures,temperature=None):
        self.partial_pressures = partial_pressures
        self.temperature = temperature
        self.pressure = sum(self.partial_pressures.values())
        self.volume = 1.0
        if self.temperature == None:
            self.temperature = 300 #K
        else:
            self.temperature = temperature
    #@profile
    def __add__(self,other):
        partial_pressures = {k: self.partial_pressures[k] + other.partial_pressures[k] for k in self.partial_pressures.iterkeys()}
        #else:
        #    partial_pressures = {}
        #    for mol in self.partial_pressures.keys():
        #        partial_pressures[mol] = 0.0
        #    for mol in other.partial_pressures.keys():
        #        partial_pressures[mol] = 0.0
        #    for mol in partial_pressures.keys():
        #        if mol in self.partial_pressures.keys():
        #            partial_pressures[mol] += self.partial_pressures[mol]
        #        if mol in other.partial_pressures.keys():
        #            partial_pressures[mol] += other.partial_pressures[mol]
        #        #del other.partial_pressures[mol]
        temperature = 0.5*(self.temperature + other.temperature)
        return Gas(partial_pressures,temperature)

    #@profile
    def __rmul__(self, other):
        #print '__rmul__ ' + str(other)
        partial_pressures = dict(self.partial_pressures)
        for key, value in partial_pressures.iteritems():
            partial_pressures[key] = value*other
        return Gas(partial_pressures,self.temperature)
        
    def __eq__(self,other):
        result = False
        if self.temperature == other.temperature:
            if self.partial_pressures == other.partial_pressures: # doesnt include 0 pressure
                result = True
        return result

    def __hash__(self):
        return None
        
    def __lt__(self,other):
        result = False
        if self.gibbs() < other.gibbs():
            result = True
        return result
    
    #def __getitem__(self,key):
    #    return self.partial_pressures[key]

    #@profile
    def gas_equlibrium_v2(self,reactions,guess = None,T=None):
        tic = time.time()
        if T == None:
            temperature=self.temperature
        else:
            temperature = T
        dimension = len(reactions)
        Pressure = self.get_pressure()
        gas_result = Gas(self.partial_pressures,temperature) # Equilibrated gas composition
        gibbs_start = gas_result.gibbs(T=temperature) # Gibbs energy of gas
        gibbs_last = gibbs_start
        gas_reaction = []        
        for reaction in reactions:
            gas_reaction.append(Gas(reaction,temperature=temperature))
        if guess == None:
            guess = np.array([0.0,0.0])
            guess_last = guess # The latest guess
            guess_test = guess # The current guess (random)
        else: # Check that the provided guess is valid
            guess = np.array(guess)
            guess_test = np.array(guess)
            gas_test = self
            for i in range(dimension):
                gas_test += guess_test[i]*gas_reaction[i]
            if gas_test.is_gas_valid():#(np.array(gas_test.partial_pressures.values()) >= 0.0 ).all(): # Check that we no negetive partial pressures
                gas_test.set_pressure(Pressure) # Scale partial pressures to overall pressure
                #gibbs_test = gas_test.gibbs(T=temperature)
                if  (self < gas_test) == False:
                    guess_last = guess_test
                    print guess_last
                    #gibbs_last = gibbs_test
                    gas_result = gas_test
                else:
                    print 'guess was NOT better than initial gas'
                    print '%.3f < %.3f is False' %(gas_test.gibbs(),  gas_result.gibbs())
                    guess = np.array([0.0,0.0]) # Consider to return an error if input guess is not accepted
                    guess_last = guess
            else:
                print 'guess was not valid' + str(guess)
                guess = np.array([0.0,0.0])
                guess_last = guess
        list_of_guess = []
        list_of_guess += [[guess_last, gas_result.gibbs()]]
        number_of_steps = 0
        number_of_succes = 0
        step = 0.01*(max((gas_reaction[0].partial_pressures.values())))/(max(self.partial_pressures.values()))
        """for re_no in range(dimension): # optimize in each direction seperately to avoid boundary problems
            step = 0.01*(max((gas_reaction[re_no].partial_pressures.values())))/(max(self.partial_pressures.values()))
            i = 0
            n = 0
            m = 0
            crit = True
            guess_is_valid = False
            while crit == True: #Optimization for 1 reaction
                n += 1
                i+=1
                guess_test = np.array(guess_last) # to prevent point reference
                if guess_is_valid == False:
                    random_step = step*(np.random.randn(1)-0.5*np.ones(1))
                guess_test[re_no] += random_step
                gas_test = self
                for re_no in range(dimension):
                    gas_test += guess_test[re_no]*gas_reaction[re_no]
                if gas_test.is_gas_valid():
                    gas_test.set_pressure(Pressure)
                    #gibbs_test =  gas_test.gibbs(T=temperature)
                    if gas_test <  gas_result:
                        guess_last = np.array(guess_test)
                        #gibbs_last = gibbs_test
                        gas_result = gas_test
                        number_of_succes +=1
                        n = 0
                        m = 0
                        guess_is_valid = True
                        list_of_guess += [[guess_last, gas_result.gibbs()]]
                    else:
                        guess_is_valid = False
                else:
                    guess_i4.18917703629s_valid = False
                if n > 10:
                    n = 0
                    m +=1
                    step *=0.5
                if m > 20 or step < 1E-4:
                    #print i, number_of_succes, step, guess_last
                    crit = False
            number_of_steps+=i"""
        i = 0
        step *= 100
        n = 0
        m = 0
        crit = True
        guess_is_valid = False
        while crit == True: # This is the only correct loop, must check if any partial pressures is 0.0
            if guess_is_valid == False:
                random_step = step*(np.random.randn(dimension) - 0.5*np.ones(dimension))
            guess_test = np.array(guess_last) + random_step
            gas_test = self
            for re_no in range(dimension):
                gas_test += guess_test[re_no]*gas_reaction[re_no]
            #gas_test = self
            #for re_no_i in range(len(reactions)):
            #    gas_test = gas_test + guess_test[re_no_i]*gas_reaction[re_no_i]
            if gas_test.is_gas_valid():
                gas_test.set_pressure(Pressure)
                #gibbs_test =  gas_test.gibbs(T=temperature)
                if gas_test <  gas_result: #comparing gibbs energy
                    guess_last = guess_test
                    #gibbs_last = gibbs_test
                    gas_result = gas_test
                    number_of_succes +=1
                    n=0
                    m=0
                    guess_is_valid = True
                    list_of_guess += [[guess_last, gas_result.gibbs()]]
                else:
                    guess_is_valid = False
            else:
                guess_is_valid = False
            if n > 10:
                n = 0
                m +=1
                step *= 0.7
            if m > 50 or step < 1E-9 or i > 5000:
                print i, m, number_of_succes, step, guess_last
                crit = False
            i += 1
            n += 1
        number_of_steps+=i
        self.guess_historic = list_of_guess
        print('Total number of steps: ' + str(number_of_succes) + ' / ' + str(number_of_steps))
        print(' step size: ' + str(step) + ', Time: ' + str(time.time()-tic) )
        return gas_result, guess_last
            

    def list_of_atoms(self):
        atom_list = {}
        for molecule in self.partial_pressures.keys():
            for element in molecule.atoms:
                if element.Z in atom_list:
                    atom_list[element.Z] += 1
                else:
                    atom_list[element.Z] = 1
        return atom_list

    def gas_atom_composition(self, option=None):
        atom_list = {}
        for molecule in self.partial_pressures.keys():
            for element in molecule.atoms:
                if element.Z in atom_list:
                    atom_list[element.Z] += self.partial_pressures[molecule]
                else:
                    atom_list[element.Z] = self.partial_pressures[molecule]
        if option == 'nomarlized':
            NOM = sum(atom_list.values())
            for element in atom_list:
                atom_list[element] /= NOM
        return atom_list

    def entropy(self,T=None):
        #do we need to introduce a entropy for the entire gas?
        if T == None:
            T=self.temperature
        S = 0.0
        for mol, value in self.partial_pressures.iteritems():
            S += value * mol.entropy(T=T)
        return S

    def enthalpy(self,T=None):
        #do we need to introduce a enthalpy for the entire gas?
        if T == None:
            T=self.temperature
        H = 0.0
        for mol, value in self.partial_pressures.iteritems():
            H += value * mol.enthalpy(T=T)
        return H

    def gibbs(self,T=None):
        #do we need to introduce a gibbs for the entire gas?
        if T == None:
            T=self.temperature
        G = 0.0
        for mol, value in self.partial_pressures.iteritems():
            if value > 0.0:
                G += self.volume * value * (mol.standard_gibbs(T=T) + R*T*np.log(value))
        return G

    def get_partial_pressure(self,test_molecule):
        result = 0.0
        if test_molecule in self.partial_pressures.keys():
            result = float(self.partial_pressures[test_molecule])
        else:
            result = 0.0
        return result #Bar
    
    def get_pressure(self):
        return sum(self.partial_pressures.values())
    
    def set_pressure(self,P_final):
        if self.is_gas_valid() == False:
            print 'Pressure is not defined for a non valid gas'
            return self            
        #P_init = sum([abs(self.partial_pressures[mol]) for mol in self.partial_pressures.keys()])
        #print 'Pressure : ' + str(P_init) + ' ' + str(self.pressure) + ' ' + str(sum(self.partial_pressures.values()))
        P_init = self.get_pressure()
        #if (np.array(self.partial_pressures.values()) < 0.0).any():
        #    print 'Negative partial_pressures'
        if P_init < 0.01:
            print 'Pressure : ' + str(P_init)
            P_init = 0.01
        for mol, value in self.partial_pressures.iteritems():
            self.partial_pressures[mol] = value * P_final / P_init
        self.volume *= P_init / P_final
        #gas_sum(result.x).partial_pressures
        self.pressure = P_final
        return self#.get_pressure()
        
    def equilibrium_constants(self,T=None):
        if T == None:
            T=self.temperature
        Q = exp(-self.gibbs(T=T)/(R*T))
        return Q
        
    def is_gas_valid(self):
        if (np.array(self.partial_pressures.values()) < 0.0 ).any():
            gas_validity = False
        else:
            gas_validity = True
        return gas_validity

    def plot_guess_historic(self,filename = None):
        fig = plt.figure()
        fig.subplots_adjust(bottom=0.2, left=0.2, right = 0.9, wspace=0.4) # Make room for x-label
        ratio = 0.61803398             # Golden mean
        ratio = 0.4                     # This figure should be very wide to span two columns
        fig_width = 14.5
        fig_width = fig_width /2.54     # width in cm converted to inches
        fig_height = fig_width*ratio
        fig.set_size_inches(fig_width,fig_height)
        axis = fig.add_subplot(1,2,1)
        Y = {}
        dimension = len(self.guess_historic[0][0])
        for i in range(dimension):
            Y[i] =  [re[0][i] for re in self.guess_historic]
        axis.plot(Y[0],Y[1],'k.-')
        axis.axhline(y=Y[1][-1],color='r')
        axis.axvline(x=Y[0][-1],color='r')
        #axis.legend(loc='upper right',prop={'size':6})
        axis.grid(False) 
        axis.tick_params(direction='in', length=6, width=1, colors='k',labelsize=8,axis='both',pad=3)
        axis.set_xlabel('reaction 1', fontsize=8)
        axis.set_ylabel('reaction 2', fontsize=8)
        
        axis2 = fig.add_subplot(1,2,2)
        axis2.semilogy()
        Y = {}
        dimension = len(self.guess_historic[0][0])
        for i in ['gibbs',]:
            Y[i] =  [re[1] for re in self.guess_historic]
        axis2.plot(range(len(Y['gibbs'])),Y['gibbs']-Y['gibbs'][-1],'k.')
        #axis.legend(loc='upper right',prop={'size':6})
        axis2.grid(False) 
        axis2.tick_params(direction='in', length=6, width=1, colors='k',labelsize=8,axis='both',pad=3)
        axis2.set_xlabel('step no', fontsize=8)
        axis2.set_ylabel('gibbs energy', fontsize=8)
        #plt.tight_layout()
        #plt.show()
        print('Saving')
        if filename == None:
            filename = 'demo_2 guess'
        plt.savefig( str(filename) +'.png',dpi=600)
        plt.clf()
        plt.close()
        del fig
        pass

if __name__ == '__main__':
    print 'Start gas.py'
    import known_molecules as km
    #gas_1 = Gas({km.CO: 0.5, km.O2: 0.5, km.CO2: 0.0},temperature=310)
    #gas_2 = Gas({km.CO: 0.0, km.O2: 0.25, km.CO2: 0.5},temperature=310)
    #gas_3 = Gas({km.CO: 0.0, km.CO2: 0.5, km.O2: 0.25},temperature=310)
    
    #gas_2 = Gas([km.CO,km.O2,km.H2,km.H2O], [0.2,0.2,0.2,0.4],pressure=2.0,temperature=310)
    #gas_3 = Gas([km.CO,km.CO2,km.O2,km.H2,km.H2O],[0.1,0.125,0.1,0.475,0.2],pressure=4.0,temperature=310)

    """gas_test = {}
    gas_test[km.CO] = 0.5
    gas_test[km.O2] = 0.5
    gas_test[km.CO2] = 0.0

    print 'Test:__eq__()'
    print 'gas_1 == gas_3: ' + str(gas_1 ==gas_3)
    print 'gas_2 == gas_3: ' + str(gas_2 ==gas_3)

    print 'Test:__add__()'
    print 'gas_1 + gas_3: ' + str(gas_1 + gas_3)
    print 'gas_2 + gas_3: ' + str(gas_2 + gas_3)
    print '(gas_1 + gas_2) + (gas_1 + gas_3): ' + str((gas_1 + gas_2) == (gas_1 + gas_3))

    print 'Test:get_pressure()'
    print 'P = ' + str(gas_1.get_pressure())
    print gas_1.partial_pressures.values()

    print 'Test:set_pressure()'
    print 'P = ' + str(gas_1.set_pressure(2.0))
    print gas_1.partial_pressures.values()
    #gas_1.gas_equlibrium(reaction,T=300)"""
    #Temperature = 1000
    Temperature = 150.0+273.15
    Pressure = 1.0
    reaction1 = {km.CO2: 0.0, km.H2: -2.0, km.CO: -1.0, km.CH3OH: 1.0, km.H2O: 0.0}
    reaction2 = {km.CO2: -1.0, km.H2: -1.0, km.CO: 1.0, km.CH3OH: 0.0, km.H2O: 1.0}
    amount_of_CO = 0.3
    amount_of_CO2 = 0.3
    amount_of_H2 = float(1.0-float(amount_of_CO+amount_of_CO2))
    gas_0 = Gas({km.H2: amount_of_H2, 
                     km.CO: amount_of_CO,
                     km.CO2: amount_of_CO2,
                     km.CH3OH: 0.0,
                     km.H2O: 0.0},temperature=Temperature)
    gas_0.set_pressure(Pressure)
    eq_result_0 = gas_0.gas_equlibrium_v2([reaction1,reaction2],guess=[0.001,0.001],T=Temperature)
    print eq_result_0
    for i in range(5):
        print('---- ' + str(i) +' ----')
        gas_0 = Gas({km.H2: amount_of_H2, 
                         km.CO: amount_of_CO,
                         km.CO2: amount_of_CO2,
                         km.CH3OH: 0.0,
                         km.H2O: 0.0},temperature=Temperature)
        gas_0.set_pressure(Pressure)
        eq_result_2 = gas_0.gas_equlibrium_v2([reaction1,reaction2],T=Temperature)
        print eq_result_2
        gas_0.plot_guess_historic()
    #gas_0.gas_equlibrium(reaction,T=310)
    #gas_reaction = Gas(reaction,temperature=Temperature)
    #gas_sum = lambda delta_P : gas_1 + delta_P*gas_reaction
    #fun = lambda delta_P: ( gas_sum(delta_P)).gibbs(T=Temperature)
    #result =  scipy.optimize.minimize_scalar(fun,bounds=(0.0, 0.01), method='bounded')
    
    
    
