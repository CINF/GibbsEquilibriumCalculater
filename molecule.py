import hashlib

import atom

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
        
    def entropy(self, T=float('Nan')):
        if T == float('Nan'):
            S = self.entropy_t0
        else:
            Tk=T/1000
            S = coef[0]*math.log(Tk)+coef[1]*Tk+coef[3]*Tk**2/2+coef[4]*Tk**3/3-coef[5]/(2*Tk**2)+coef[7]
        return S

    def enthalpy(self, T=float('Nan')):
        if T == float('Nan'):
            H = self.entropy_t0
        else:
            Tk=T/1000
            H = coef[0]*Tk+coef[1]*Tk**2/2+coef[3]*Tk**3/3+coef[4]*Tk**4/4-coef[5]/Tk+coef[6]-coef[8]+self.M.enthalpy_t0
        return H

if __name__ == '__main__':
    m = Molecule([atom.Atom(5), atom.Atom(2), atom.Atom(3), atom.Atom(2)])
    l = Molecule([atom.Atom(2), atom.Atom(5), atom.Atom(3), atom.Atom(2)])
    n = Molecule([atom.Atom(6), atom.Atom(8)])

    print n == m
    print n == l
    print l == m
 
