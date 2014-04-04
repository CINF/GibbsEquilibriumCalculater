
import molecule
import scipy.optimize
import copy
import numpy as np
R = 8.314462175
pi = 3.141592653589793
e = 2.718281828459045

def gas_mixing(main,delta,):
     1.0*gas_reaction
     gas_0 = Gas(self.partial_pressures,temperature)
     gas_sum = lambda delta_P : ((gas_0 + delta_P*gas_reaction))#.set_pressure(Pressure))
     validate = lambda delta_P: gas_sum(delta_P) #(np.array(gas_sum(delta_P).partial_pressures.values()) < 0.0).any()*gas_0 + (np.array(gas_sum(delta_P).partial_pressures.values()) > 0.0).all()*gas_sum(delta_P)
     fun = lambda delta_P: ( validate(delta_P)).gibbs(T=temperature)
     result =  scipy.optimize.minimize_scalar(fun,bounds=(0.0, 0.001), method='Bounded')

class Gas():
    def __init__(self,partial_pressures,temperature=None):
        self.partial_pressures = partial_pressures
        self.temperature = temperature
        self.volume = 1.0
        if self.temperature == None:
            self.temperature = 300 #K

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

    def gas_equlibrium(self,reactions,T=None):
        if T == None:
            temperature=self.temperature
        else:
            temperature = T
        Pressure = self.get_pressure()
        gas_reaction = []
        for reaction in reactions:
            gas_reaction.append(Gas(reaction,temperature=temperature))
        #if (np.array(gas_reaction.gas_atom_composition().values()) != 0.0).any():
        #    print 'Error: ' + str(gas_reaction.gas_atom_composition())
        #    print 'Error specs: ' + str(np.array((gas_reaction.gas_atom_composition().values()) != 0.0))
        #1.0*gas_reaction
        gas_0 = Gas(self.partial_pressures,temperature)
        gas_sum = lambda delta_P : ((gas_0 + delta_P[0]*gas_reaction[0] + delta_P[1]*gas_reaction[1]).set_pressure(Pressure))
        validate = lambda delta_P: gas_0 if (np.array(gas_sum(delta_P).partial_pressures.values()) < 1e-12).any() else gas_sum(delta_P)
        fun = lambda delta_P: ( validate(delta_P)).gibbs(T=temperature)
        #result =  scipy.optimize.minimize_scalar(fun,bounds=(0.0, 0.25), method='Bounded')
        x0 = [0.001,0.001]
        result =  scipy.optimize.minimize(fun,x0, method='L-BFGS-B',options={'gtol': 1e-9, 'disp': False})
        #gx = 0.0
        #G = fun(gx)
        #for i in np.linspace(0.0,0.25,10):
        #    if fun(i) < G:
        #        G=fun(i)
        #        gx = i
        #print 'fit params: ' + str(result.x) + ' i: ' + str(gx) + ' at T: ' + str(temperature) + ' P: ' + str(validate(gx).get_pressure())
        #print 'fit params: ' + str(result.x)
        return validate(result.x)

    def list_of_atoms(self):
        atom_list = {}
        for molecule in self.partial_pressures.keys():
            for element in molecule.atoms:
                if element.Z in atom_list:
                    atom_list[element.Z] += 1
                else:
                    atom_list[element.Z] = 1
        return atom_list

    def gas_atom_composition(self):
        atom_list = {}
        for molecule in self.partial_pressures.keys():
            for element in molecule.atoms:
                if element.Z in atom_list:
                    atom_list[element.Z] += self.partial_pressures[molecule]
                else:
                    atom_list[element.Z] = self.partial_pressures[molecule]
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
        #total_pressure = sum(self.partial_pressures.values())
        for mol in self.partial_pressures.keys():
            G += self.volume*self.partial_pressures[mol] * (mol.standard_gibbs(T=T) + 
                R*T * np.log(abs(self.partial_pressures[mol]+1e-15))
                )
        #G0 = self.enthalpy(T=T) - T * self.entropy(T=T)
        #if abs(G-G0) > 1e-10:
        #    print 'Error: ' + str(G) + ' ' + str(G0)
        #    print 'Error specs: ' + str((G-G0)/G)
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
        P_init = abs(sum(self.partial_pressures.values()))
        for mol in self.partial_pressures.keys():
            self.partial_pressures[mol] *= P_final / P_init
        self.volume *= P_init / P_final
        return self#.get_pressure()
    def equilibrium_constants(self,T=None):
        if T == None:
            T=self.temperature
        Q = exp(-self.gibbs(T=T)/(R*T))
        return Q
        

if __name__ == '__main__':
    print 'Start gas.py'
    import known_molecules as km
    gas_1 = Gas({km.CO: 0.5, km.O2: 0.5, km.CO2: 0.0},temperature=310)
    gas_2 = Gas({km.CO: 0.0, km.O2: 0.25, km.CO2: 0.5},temperature=310)
    gas_3 = Gas({km.CO: 0.0, km.CO2: 0.5, km.O2: 0.25},temperature=310)
    
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
    Temperature = 1000
    gas_0 = Gas({km.H2: 0.99, km.CO: 0.01, km.H2O: 0.0,km.CH4: 0.0},temperature=Temperature)
    gas_0.equilibrium_constants(T=300)
    print gas_0.gibbs(T=300)
    print gas_0.gas_atom_composition()
    reaction = {km.H2: -3.0, km.CO: -1.0, km.H2O: 1.0,km.CH4: 1.0}
    gas_0.gas_equlibrium(reaction,T=310)
    #gas_reaction = Gas(reaction,temperature=Temperature)
    #gas_sum = lambda delta_P : gas_1 + delta_P*gas_reaction
    #fun = lambda delta_P: ( gas_sum(delta_P)).gibbs(T=Temperature)
    #result =  scipy.optimize.minimize_scalar(fun,bounds=(0.0, 0.01), method='bounded')
    
    
    
