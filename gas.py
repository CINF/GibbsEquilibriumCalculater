
import molecule
class Gas():
    def __init__(self,partial_pressures,temperature=None):
        self.partial_pressures = partial_pressures
        self.temperature = temperature
        if self.temperature == None:
            self.temperature = 300 #K

    def __add__(self,other):
        partial_pressures = self.partial_pressures
        for mol in other.partial_pressures.keys():
            try:
                partial_pressures[mol] += other.partial_pressures[mol]
            except KeyError:
                partial_pressures[mol] = other.partial_pressures[mol]
            #del other.partial_pressures[mol]
        temperature = 0.5*(self.temperature + other.temperature)
        return Gas(partial_pressures,temperature)

    def __eq__(self,other):
        if self.temperature == other.temperature:
            if self.partial_pressures == other.partial_pressures: # doesnt include 0 pressure
                return True
        return False
    def __hash__(self):
        return None

    def gas_equlibrium(self):
        # calculate the equlibrium gas composition given from the initial gas composition, temperature and pressure
        # minimize self.gibbs() by adjusting partial pressure, without changing self.gas_atom_composition()
        # line1: equlibrium constant - equlibrium state=0
        # line2: conservation self.gas_atom_composition() = initial_partial_pressure
        # line3: --||--
        # line4: sum(self.partial_pressure)=1.0
        # self.assign_pressure_to_molecule()
        
        return Gas(self.partial_pressures,self.temperature)

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
                    atom_list[element.Z] += molecule.partial_pressure
                else:
                    atom_list[element.Z] = molecule.partial_pressure
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
            G += self.partial_pressures[mol] * mol.gibbs(T=T)
        return G

    def get_partial_pressure(self,test_molecule):
        result = 0.0
        try:
            result = self.partial_pressures[test_molecule]
        except KeyError:
            result = 0.0
        return result #Bar

if __name__ == '__main__':
    print 'Start gas.py'
    import known_molecules as km
    gas_1 = Gas({km.CO: 0.5, km.O2: 0.5, km.CO2: 0.0},temperature=310)
    gas_2 = Gas({km.CO: 0.0, km.O2: 0.25, km.CO2: 0.5},temperature=310)
    gas_3 = Gas({km.CO: 0.0, km.CO2: 0.5, km.O2: 0.25},temperature=310)
    #gas_2 = Gas([km.CO,km.O2,km.H2,km.H2O], [0.2,0.2,0.2,0.4],pressure=2.0,temperature=310)
    #gas_3 = Gas([km.CO,km.CO2,km.O2,km.H2,km.H2O],[0.1,0.125,0.1,0.475,0.2],pressure=4.0,temperature=310)

    gas_test = {}
    gas_test[km.CO] = 0.5
    gas_test[km.O2] = 0.5
    gas_test[km.CO2] = 0.0

    print 'Test:__eq__()'
    print 'gas_1 == gas_3: ' + str(gas_1 ==gas_3)
    print 'gas_2 == gas_3: ' + str(gas_2 ==gas_3)\

    print 'Test:__add__()'
    print 'gas_1 + gas_3: ' + str(gas_1 + gas_3)
    print 'gas_2 + gas_3: ' + str(gas_2 + gas_3)
    print '(gas_1 + gas_2) + (gas_1 + gas_3): ' + str((gas_1 + gas_2) == (gas_1 + gas_3))
    #print gas_1.list_of_atoms()
    #print gas_1.gas_atom_composition()

    #print gas_1.gas_equlibrium()
    #print 'Gibbs: ' + str(gas_1.gibbs())

    #print 'add test'
    #gas_4 = gas_1 + gas_2
    #print gas_4.get_partial_pressure(km.O2)
    
    
