import atom

class Molecule():
    def __init__(self, atoms):
        print type(atoms)
        for a in atoms:
            print isinstance(a, atom.Atom)
            print a.symbol()
            
            
if __name__ == '__main__':
    m = Molecule([atom.Atom(2), atom.Atom(3)])
