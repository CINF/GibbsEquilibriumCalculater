
import molecule
class Gas():
    def __init__(self,molecules,partial_pressure_ratio=None,pressure=None,temperature=None):
        self.molecules = molecules
        self.partial_pressure_ratio = partial_pressure_ratio
        self.pressure = pressure
        self.temperature = temperature
        if self.partial_pressure_ratio == None:
            print 'ERROR - Partial pressure needed'
        if self.pressure == None:
            self.pressure = 1.0 #bar
        if self.temperature == None:
            self.temperature = 300 #K
        for i in range(len(self.molecules)):
            self.molecules[i].partial_pressure = partial_pressure_ratio[i]

    def gas_equlibrium(self):
        # calculate the equlibrium gas composition given from the initial gas composition, temperature and pressure
        # minimize self.gibbs() by adjusting partial pressure, without changing self.gas_atom_composition()
        # line1: equlibrium constant - equlibrium state=0
        # line2: conservation self.gas_atom_composition() = initial_partial_pressure
        # line3: --||--
        # line4: sum(self.partial_pressure)=1.0
        # self.assign_pressure_to_molecule()
        
        return self.partial_pressure_ratio
    
    def set_partial_pressure_ratio(self,ppr):
        if sum(ppr) == 1.0:
            self.partial_pressure_ratio = ppr
        else:
            print 'ERROR - sum of pressure ratio should be 1.0'
            self.partial_pressure_ratio = ppr / sum(ppr)
        self.assign_pressure_to_molecule()
        return self.partial_pressure_ratio

    def list_of_atoms(self):
        atom_list = {}
        for molecule in self.molecules:
            for element in molecule.atoms:
                if element.Z in atom_list:
                    atom_list[element.Z] += 1
                else:
                    atom_list[element.Z] = 1
        return atom_list

    def gas_atom_composition(self):
        self.assign_pressure_to_molecule()
        atom_list = {}
        for molecule in self.molecules:
            for element in molecule.atoms:
                if element.Z in atom_list:
                    atom_list[element.Z] += molecule.partial_pressure
                else:
                    atom_list[element.Z] = molecule.partial_pressure
        return atom_list

    def assign_pressure_to_molecule(self):
        for i in range(len(self.molecules)):
            self.molecules[i].partial_pressure = self.partial_pressure_ratio[i]
        return True

    def entropy(self,T=None):
        #do we need to introduce a entropy for the entire gas?
        if T == None:
            T=self.temperature
        S = 0.0
        for molecule in self.molecules:
            S += molecule.partial_pressure * molecule.entropy(T=T)
        return S

    def enthalpy(self,T=None):
        #do we need to introduce a enthalpy for the entire gas?
        if T == None:
            T=self.temperature
        H = 0.0
        for molecule in self.molecules:
            H += molecule.partial_pressure * molecule.enthalpy(T=T)
        return H

    def gibbs(self,T=None):
        #do we need to introduce a gibbs for the entire gas?
        if T == None:
            T=self.temperature
        G = 0.0
        for molecule in self.molecules:
            G += molecule.partial_pressure * molecule.gibbs(T=T)
        return G

if __name__ == '__main__':
    print 'Start gas.py'
    import known_molecules as km
    gas_1 = Gas([km.CO,km.CO2,km.H2,km.H2O],[0.0,0.25,0.75,0.0],pressure=2.0,temperature=310)
    
    print gas_1.list_of_atoms()
    print gas_1.gas_atom_composition()

    print gas_1.gas_equlibrium()
    print 'Gibbs: ' + str(gas_1.gibbs())
    
