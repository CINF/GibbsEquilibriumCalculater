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
        pass

    def list_of_atoms(self):
        atom_list = {}
        for element in self.atoms:
            if element.z in atom_list:
                atom_list[element.z] += 1
            else:
                atom_list[element.z] = 1
        return atom_list


class Molecule(SimpleMolecule):
    


    def gibbs_free_energy(self):
        gibbs = -1
        if self == Molecule([atom.Atom(6), atom.Atom(8)]):
            gibbs = -110500.0


if __name__ == '__main__':
    m = Molecule([atom.Atom(5), atom.Atom(2), atom.Atom(3), atom.Atom(2)])
    l = Molecule([atom.Atom(5), atom.Atom(2), atom.Atom(3), atom.Atom(2)])
    n = Molecule([atom.Atom(6), atom.Atom(8)])
    
    CO = SimpleMolecule([atom.Atom(6), atom.Atom(8)])
    
    a = {}
    a[CO]['gibbs'] = -110500.0
    
    print n == m
    print n == l
    print l == m
    m.gibbs_free_energy()
    n.gibbs_free_energy()
    
    #print m.list_of_atoms()
