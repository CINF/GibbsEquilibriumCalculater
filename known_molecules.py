import atom
import simple_molecule

H2    = simple_molecule.SimpleMolecule([atom.Atom(1), atom.Atom(1)])
H2O   = simple_molecule.SimpleMolecule([atom.Atom(1), atom.Atom(1), atom.Atom(8)])
H2O2  = simple_molecule.SimpleMolecule([atom.Atom(1), atom.Atom(1), atom.Atom(8), atom.Atom(8)])
O2    = simple_molecule.SimpleMolecule([atom.Atom(8), atom.Atom(8)])
CO    = simple_molecule.SimpleMolecule([atom.Atom(6), atom.Atom(8)])
CO2   = simple_molecule.SimpleMolecule([atom.Atom(6), atom.Atom(8), atom.Atom(8)])
CH4   = simple_molecule.SimpleMolecule([atom.Atom(6), atom.Atom(1), atom.Atom(1), atom.Atom(1), atom.Atom(1)])
CH3OH = simple_molecule.SimpleMolecule([atom.Atom(6), atom.Atom(1), atom.Atom(1), atom.Atom(1), atom.Atom(8), atom.Atom(1)])
NH3   = simple_molecule.SimpleMolecule([atom.Atom(7), atom.Atom(1), atom.Atom(1), atom.Atom(1)])
N2    = simple_molecule.SimpleMolecule([atom.Atom(7), atom.Atom(7)])
NO    = simple_molecule.SimpleMolecule([atom.Atom(7), atom.Atom(8)])
NO2   = simple_molecule.SimpleMolecule([atom.Atom(7), atom.Atom(8), atom.Atom(8)])
