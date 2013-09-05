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


H2.enthalpy    = 0.0
O2.enthalpy    = 0.0
N2.enthalpy    = 0.0
CO.enthalpy    = -110500.0
H2O.enthalpy   = -241830.0
CO2.enthalpy   = -393500.0
CH4.enthalpy   = -74870.0
CH3OH.enthalpy = -201300.0
NH3.enthalpy   = -45940.0

H2.entropy    = 130.679
CO.entropy    = 197.7
H2O.entropy   = 188.84
CO2.entropy   = 213.7
CH4.entropy   = 186.25
CH3OH.entropy = 239.9
N2.entropy    = 191.61
NH3.entropy   = 192.778

