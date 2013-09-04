import atom
import molecule

H2    = molecule.Molecule([atom.Atom(1), atom.Atom(1)])
H2O   = molecule.Molecule([atom.Atom(1), atom.Atom(1), atom.Atom(8)])
H2O2  = molecule.Molecule([atom.Atom(1), atom.Atom(1), atom.Atom(8), atom.Atom(8)])
O2    = molecule.Molecule([atom.Atom(8), atom.Atom(8)])
CO    = molecule.Molecule([atom.Atom(6), atom.Atom(8)])
CO2   = molecule.Molecule([atom.Atom(6), atom.Atom(8), atom.Atom(8)])
CH4   = molecule.Molecule([atom.Atom(6), atom.Atom(1), atom.Atom(1), atom.Atom(1), atom.Atom(1)])
CH3OH = molecule.Molecule([atom.Atom(6), atom.Atom(1), atom.Atom(1), atom.Atom(1), atom.Atom(8), atom.Atom(1)])
NH3   = molecule.Molecule([atom.Atom(7), atom.Atom(1), atom.Atom(1), atom.Atom(1)])
N2    = molecule.Molecule([atom.Atom(7), atom.Atom(7)])
NO    = molecule.Molecule([atom.Atom(7), atom.Atom(8)])
NO2   = molecule.Molecule([atom.Atom(7), atom.Atom(8), atom.Atom(8)])
