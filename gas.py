
import molecule
import scipy.optimize
import copy
import numpy as np
import known_molecules as km
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

    def __add__(self,other):
        partial_pressures = copy.deepcopy(self.partial_pressures)
        for mol in other.partial_pressures.keys():
            if mol in self.partial_pressures.keys():
                partial_pressures[mol] += other.partial_pressures[mol]
            else:
                partial_pressures[mol] = other.partial_pressures[mol]
            #del other.partial_pressures[mol]
        temperature = 0.5*(self.temperature + other.temperature)
        return Gas(partial_pressures,temperature)

    def __rmul__(self, other):
        #print '__rmul__ ' + str(other)
        partial_pressures = copy.deepcopy(self.partial_pressures)
        for mol in self.partial_pressures.keys():
            partial_pressures[mol] *= other
        return Gas(partial_pressures,self.temperature)
        
    def __eq__(self,other):
        if self.temperature == other.temperature:
            if self.partial_pressures == other.partial_pressures: # doesnt include 0 pressure
                return True
        return False

    def __hash__(self):
        return None

        
    def gas_equlibrium_v2(self,reactions,guess = None,T=None):
        if T == None:
            temperature=self.temperature
        else:
            temperature = T
        Pressure = self.get_pressure()
        gas_0 = Gas(self.partial_pressures,temperature) # Initial gas composition
        gas_result = gas_0 # Equilibrated gas composition
        gibbs_start = gas_result.gibbs(T=temperature) # Gibbs energy of gas
        gibbs_last = gibbs_start
        gas_reaction = []
        for reaction in reactions:
            gas_reaction.append(Gas(reaction,temperature=temperature))
        if guess == None:
            guess = [0.0,0.0]
            guess_last = guess # The latest guess
            guess_test = guess # The current guess (random)
        else: # Check that the provided guess is valid
            guess_test = guess
            gas_temp = (gas_0 + guess[0]*gas_reaction[0] + guess[1]*gas_reaction[1])
            if (np.array(gas_temp.partial_pressures.values()) >= 0.0 ).all(): # Check that we no negetive partial pressures
                gas_test = gas_temp.set_pressure(Pressure) # Scale partial pressures to overall pressure
                gibbs_test = gas_test.gibbs(T=temperature)
                if gibbs_test <  gibbs_last:
                    guess_last = guess
                    gibbs_last = gibbs_test
                    gas_result = gas_test
                else:
                    guess = [0.0,0.0] # Consider to return an error if input guess is not accepted
                    guess_last = guess
            else:
                guess = [0.0,0.0]
                guess_last = guess
        number_of_steps = 0
        number_of_succes = 0
        step = 0.01
        n = 0
        m = 0
        for i in range(1000):
            guess_test = np.array(guess_last) + (step*(np.random.randn(1))) # Make a random pertubation
            gas_test = (gas_0 + guess_test[0]*gas_reaction[0])
            n += 1
            if (np.array(gas_test.partial_pressures.values()) >= 0.0 ).all(): # Check that we no negetive partial pressures
                gas_test.set_pressure(Pressure)
                gibbs_test =  gas_test.gibbs(T=temperature)
                if gibbs_test <  gibbs_last:
                    guess_last = guess_test
                    gibbs_last = gibbs_test
                    gas_result = gas_test
                    number_of_succes +=1
                    n=0
                    m=0
            if n > 10:
                n = 0
                m +=1
                step *=0.5
            if m > 5 or step < 1E-9:
                break
        number_of_steps+=i
        step = 0.01
        n = 0
        m = 0
        for i in range(1000):
            guess_test = np.array(guess_last) + (step*(np.random.randn(1)))
            gas_test = (gas_0 + guess_test[1]*gas_reaction[1])
            n += 1
            if (np.array(gas_test.partial_pressures.values()) >= 0.0 ).all():
                gas_test.set_pressure(Pressure)
                gibbs_test =  gas_test.gibbs(T=temperature)
                if gibbs_test <  gibbs_last:
                    guess_last = guess_test
                    gibbs_last = gibbs_test
                    gas_result = gas_test
                    number_of_succes +=1
                    n=0
                    m=0
            if n > 10:
                n = 0
                m +=1
                step *=0.5
            if m > 5 or step < 1E-9:
                break
        number_of_steps+=i
        step = 0.1
        n = 0
        m = 0
        for i in range(5000):
            guess_test = np.array(guess_last) + (step*(np.random.randn(2)))
            gas_test = (gas_0 + guess_test[0]*gas_reaction[0] + guess_test[1]*gas_reaction[1])
            n += 1
            if (np.array(gas_test.partial_pressures.values()) >= 0.0 ).all():
                gas_test.set_pressure(Pressure)
                gibbs_test =  gas_test.gibbs(T=temperature)
                if gibbs_test <  gibbs_last:
                    guess_last = guess_test
                    gibbs_last = gibbs_test
                    gas_result = gas_test
                    number_of_succes +=1
                    n=0
                    m=0
            if n > 10:
                n = 0
                m +=1
                step *=0.5
            if m > 10 or step < 1E-9:
                break
        number_of_steps+=i
        print 'Total number of steps: ' + str(number_of_succes) + ' / ' + str(number_of_steps)
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
        for mol in self.partial_pressures.keys():
            S += self.partial_pressures[mol] * mol.entropy(T=T)
        return S

    def enthalpy(self,T=None):
        #do we need to introduce a enthalpy for the entire gas?
        if T == None:
            T=self.temperature
        H = 0.0
        for mol in self.partial_pressures.keys():
            H += self.partial_pressures[mol] * mol.enthalpy(T=T)
        return H

    def gibbs(self,T=None):
        #do we need to introduce a gibbs for the entire gas?
        if T == None:
            T=self.temperature
        G = 0.0
        for mol in self.partial_pressures.keys():
            #if self.partial_pressures[mol] < 0.0:
            #    G = float('nan')
            G += self.volume*self.partial_pressures[mol] * (mol.standard_gibbs(T=T) + 
                R*T * np.log(abs(self.partial_pressures[mol]+1e-15))
                )
        return G

    def get_partial_pressure(self,test_molecule):
        result = 0.0
        try:
            result = float(self.partial_pressures[test_molecule])
        except KeyError:
            result = 0.0
        return result #Bar
    
    def get_pressure(self):
        return sum(self.partial_pressures.values())
    
    def set_pressure(self,P_final):
        #P_init = sum([abs(self.partial_pressures[mol]) for mol in self.partial_pressures.keys()])
        #print 'Pressure : ' + str(P_init) + ' ' + str(self.pressure) + ' ' + str(sum(self.partial_pressures.values()))
        P_init = sum(self.partial_pressures.values())
        #if (np.array(self.partial_pressures.values()) < 0.0).any():
        #    print 'Negative partial_pressures'
        if P_init < 0.01:
            print 'Pressure : ' + str(P_init)
            P_init = 0.01
        for mol in self.partial_pressures.keys():
            self.partial_pressures[mol] *= P_final / P_init
        self.volume *= P_init / P_final
        #gas_sum(result.x).partial_pressures
        self.pressure = P_final
        return self#.get_pressure()
    def equilibrium_constants(self,T=None):
        if T == None:
            T=self.temperature
        Q = exp(-self.gibbs(T=T)/(R*T))
        return Q
        

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
    Pressure = 2.5
    reaction1 = {km.CO2: 0.0, km.H2: -2.0, km.CO: -1.0, km.CH3OH: 1.0, km.H2O: 0.0}
    reaction2 = {km.CO2: -1.0, km.H2: -1.0, km.CO: 1.0, km.CH3OH: 0.0, km.H2O: 1.0}
    amount_of_CO = 0.1
    amount_of_CO2 = 0.1
    amount_of_H2 = float(1.0-float(amount_of_CO+amount_of_CO2))
    gas_0 = Gas({km.H2: amount_of_H2, 
                     km.CO: amount_of_CO,
                     km.CO2: amount_of_CO2,
                     km.CH3OH: 0.0,
                     km.H2O: 0.0},temperature=Temperature)
    gas_0 = gas_0.set_pressure(Pressure)
    gas_1 = gas_0.gas_equlibrium([reaction1,reaction2],T=Temperature)
    print gas_1.partial_pressures.values()
    gas_2 = gas_0.gas_equlibrium_v2([reaction1,reaction2],T=Temperature)
    print gas_2.partial_pressures.values()
    #gas_0.gas_equlibrium(reaction,T=310)
    #gas_reaction = Gas(reaction,temperature=Temperature)
    #gas_sum = lambda delta_P : gas_1 + delta_P*gas_reaction
    #fun = lambda delta_P: ( gas_sum(delta_P)).gibbs(T=Temperature)
    #result =  scipy.optimize.minimize_scalar(fun,bounds=(0.0, 0.01), method='bounded')
    
    
    
