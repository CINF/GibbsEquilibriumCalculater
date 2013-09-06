import atom
import molecule
import math

def Gibbs(M,T):
  Tk=T/1000
  coef = None
  if T<298:
    return False
  if M==H2 and T>=298 and <1000:
    coef = [33.066178,-11.363417,11.432816,-2.772874,-0.158558,-9.980797,172.707974,0.0]
  H=coef[0]*Tk+coef[1]*Tk**2/2+coef[3]*Tk**3/3+coef[4]*Tk**4/4-coef[5]/Tk+coef[6]-coef[8]+M.enthalpy_0
  S=coef[0]*math.log(Tk)+coef[1]*Tk+coef[3]*Tk**2/2+coef[4]*Tk**3/3-coef[5]/(2*Tk**2)+coef[7]
  return H-T*S

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




H2.enthalpy_0    = 0.0
O2.enthalpy_0    = 0.0
N2.enthalpy_0    = 0.0
CO.enthalpy_0    = -110500.0
H2O.enthalpy_0   = -241830.0
CO2.enthalpy_0   = -393500.0
CH4.enthalpy_0   = -74870.0
CH3OH.enthalpy_0 = -201300.0
NH3.enthalpy_0   = -45940.0

H2.entropy_0    = 130.679
CO.entropy_0    = 197.7
H2O.entropy_0   = 188.84
CO2.entropy_0   = 213.7
CH4.entropy_0   = 186.25
CH3OH.entropy_0 = 239.9
N2.entropy_0    = 191.61
NH3.entropy_0   = 192.778

