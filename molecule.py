import hashlib

import atom

class SimpleMolecule():
    def __init__(self, atoms):
        self.atoms = atoms

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


class Molecule(SimpleMolecule):
    
    def gibbs_free_energy(self):
        gibbs = -1
        if self == Molecule([atom.Atom(6), atom.Atom(8)]):
            gibbs = -110500.0
        return gibbs


if __name__ == '__main__':
    m = Molecule([atom.Atom(5), atom.Atom(2), atom.Atom(3), atom.Atom(2)])
    l = Molecule([atom.Atom(2), atom.Atom(5), atom.Atom(3), atom.Atom(2)])
    n = Molecule([atom.Atom(6), atom.Atom(8)])
    
    CO = SimpleMolecule([atom.Atom(6), atom.Atom(8)])
    
    print n == m
    print n == l
    print l == m
    print m.gibbs_free_energy()
    print n.gibbs_free_energy()
 
