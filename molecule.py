import atom
from simple_molecule import SimpleMolecule
from known_molecules import *


class Molecule(SimpleMolecule):
    
    def enthalpy(self):
        # As this table gets larger, the actual data should properly go in a file
        gibbs = float('NaN')
        if self == CO:
            gibbs = -110500.0
        if self == H2:
            gibbs = 0.0
        if self == H2O:
            gibbs = -241830.0
        if self == O2:
            gibbs = 0.0
        if self == CO2:
            gibbs = -393500.0
        if self == CH4:
            gibbs = -74870.0
        if self == CH3OH:
            gibbs = -201300.0
        if self == NH3:
            gibbs = -45940.0
        if self == N2:
            gibbs = 0.0            
        return gibbs

    def entropy(self):
        # As this table gets larger, the actual data should properly go in a file
        entropy = float('NaN')
        if self == CO:
            entropy = 197.7
        if self == H2:
            entropy = 130.679
        if self == H2O:
            entropy = 188.84
        if self == CO2:
            entropy = 213.7
        if self == CH4:
            entropy = 186.25
        if self == CH3OH:
            entropy = 239.9
        if self == NH3:
            entropy = 192.778
        if self == N2:
            entropy = 191.61           
        return gibbs

if __name__ == '__main__':
    m = Molecule([atom.Atom(5), atom.Atom(2), atom.Atom(3), atom.Atom(2)])
    l = Molecule([atom.Atom(2), atom.Atom(5), atom.Atom(3), atom.Atom(2)])
    n = Molecule([atom.Atom(6), atom.Atom(8)])

    print n == m
    print n == l
    print l == m
    print m.enthalpy()
    print n.enthalpy()
 
