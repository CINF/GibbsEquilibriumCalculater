import atom
import molecule

class ChemicalReaction():
    def __init__(self, left, right):
        self.left_side  = left
        self.right_side = right

if __name__ == '__main__':

    CO = molecule.Molecule([atom.Atom(6), atom.Atom(8)])
    CO2 = molecule.Molecule([atom.Atom(6), atom.Atom(8), atom.Atom(8)])
    
    reaction = ChemicalReaction([CO], [CO2])
