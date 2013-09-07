import atom
import molecule
import known_molecules as km

class ChemicalReaction():
    def __init__(self, left, right):
        self.left_side  = left
        self.right_side = right
    
    def list_of_atoms(self, side):
        atoms = {}
        molecules = {}
        if side == 'left':
            molecules = self.left_side
        if side == 'right':
            molecules = self.right_side

        atoms = {}
        for m in molecules:
            elements = m.list_of_atoms()
            for element in elements.keys():
                if element in atoms:
                    atoms[element] += elements[element]
                else:
                    atoms[element] = elements[element]
        return atoms
    
    def is_balanced(self):
        left  = self.list_of_atoms('left')
        right = self.list_of_atoms('right')
        return left == right
        
if __name__ == '__main__':

    reaction = ChemicalReaction([km.CO,km.O2,km.CO2], [km.CO2])
    print reaction.list_of_atoms('left')
    print reaction.list_of_atoms('right')

    print reaction.is_balanced()
    """
    left_atoms = {}
    for mol in reaction.left_side:
        elements = mol.list_of_atoms()
        for element in elements.keys():
            if element in left_atoms:
                left_atoms[element] = left_atoms[element] + elements[element]
            else:
                left_atoms[element] = elements[element]

    print left_atoms
    """
