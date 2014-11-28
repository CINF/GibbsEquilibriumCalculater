import hashlib

import atom

import math

class Molecule():
    def __init__(self, atoms):
        self.atoms = atoms
        self.entropy_t0  = float('NaN')
        self.enthalpy_t0 = float('NaN')
        self.coeff    = []

    def __eq__(self, other):
        if self.list_of_atoms() == other.list_of_atoms():
            return True
        else:
            return False

    def __hash__(self):
        hash_value = ''
        for a in self.atoms:
            hash_value = hashlib.sha224(a.symbol() + hash_value).hexdigest()
        return int(hash_value,16)

    def list_of_atoms(self):
        atom_list = {}
        for element in self.atoms:
            if element.Z in atom_list:
                atom_list[element.Z] += 1
            else:
                atom_list[element.Z] = 1
        return atom_list
        
    def standard_entropy(self, T=None):
        if T == None:
            S = self.entropy_t0#J/mol/K
        else:
            Tk=float(T)/1000.0
            S = ((self.coeff[0]*math.log(Tk)+
                self.coeff[1]*Tk+
                self.coeff[2]*(Tk**2)/2.0+
                self.coeff[3]*(Tk**3)/3.0-
                self.coeff[4]/(2.0*(Tk**2))+
                self.coeff[6]))
            #S = S#J/mol/K
        return S #self.entropy_t0#S #J/mol/K

    def standard_enthalpy(self, T=None):
        if T == None:
            H = self.enthalpy_t0
        else:
            Tk=float(T)/1000.0
            H = ((self.coeff[0]*(Tk) + 
                self.coeff[1]*(Tk**2)/2.0+
                self.coeff[2]*(Tk**3)/3.0+
                self.coeff[3]*(Tk**4)/4.0-
                self.coeff[4]/(Tk)+
                self.coeff[5]-
                0*self.coeff[7]))
        return H*1000.0 + 0*self.enthalpy_t0 #J/mol/K
    
    def standard_gibbs(self, T=None):
        if T == None:
            G = self.standard_enthalpy()-298.15*self.standard_entropy()
        else:
            Tk=float(T)/1000
            G = self.standard_enthalpy(T)-T*self.standard_entropy(T)
        return G #J/mol/K

if __name__ == '__main__':
    import known_molecules as km
    
    m = Molecule([atom.Atom(5), atom.Atom(2), atom.Atom(3), atom.Atom(2)])
    l = Molecule([atom.Atom(2), atom.Atom(5), atom.Atom(3), atom.Atom(2)])
    n = Molecule([atom.Atom(6), atom.Atom(8)])

    print km.H2.standard_gibbs(300)

    print n == m
    print n == l
    print l == m
 
