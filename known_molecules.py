import atom
import molecule
import math

def Gibbs(M,T):
  #Shoomate equation
  # solutions found for EtOH, CO2, CO, H2O, N2, O2
  # link to NIST
  # #http://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI&Mask=1#Thermo-Gas
  # this page might be of interest:
  # http://jkitchin.github.io/blog/tag/thermodynamics/
  Tk=T/1000
  coef = None
  if T<298:
    return False
  if M==H2 and T>=298 and T<1000:
    coef = [33.066178,-11.363417,11.432816,-2.772874,-0.158558,-9.980797,172.707974,0.0] #H2
    
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

H2.coeff = [33.066178,-11.363417,11.432816,-2.772874,-0.158558,-9.980797,172.707974,0.0]
H2O.coeff = [30.09200,  6.832514,     6.793435,  -2.534480,  0.082139, -250.8810, 223.3967,   -241.8264]  # H2O
CO2.coeff = [24.99735,  55.18696,   -33.69137,    7.948387, -0.136638, -403.6075, 228.2431,   -393.5224]  # CO2
CO.coeff = [25.56759,  6.096130,     4.054656,  -2.671301,  0.131021, -118.0089, 227.3665,   -110.5271]  # CO

O2.coeff = [31.32234,  -20.23531,    57.86644,  -36.50624,  -0.007374, -8.903471, 246.7945,   0.0]  # O2
CH4.coeff = [-0.703029,  108.4773,    -42.52157,  5.862788,  0.678565, -76.84376, 158.7163,   -74.87310]  # CH4
NH3.coeff = [19.99563,  49.77119,   -15.37599,  1.921168, 0.189174, -53.30667, 203.8591,  -45.89806]  # NH3

CH3OH.coeff = [14.1952,  97.7218,  -9.73279,  -12.8461,0.15819 , -209.037, 0.0, -201.102]  # NH3, coeff[6] is missing

H2.enthalpy_t0    = 0.0
O2.enthalpy_t0    = 0.0
N2.enthalpy_t0    = 0.0
CO.enthalpy_t0    = -110500.0
H2O.enthalpy_t0   = -241830.0
CO2.enthalpy_t0   = -393500.0
CH4.enthalpy_t0   = -74870.0
CH3OH.enthalpy_t0 = -201300.0
NH3.enthalpy_t0   = -45940.0

H2.entropy_t0    = 130.679
CO.entropy_t0    = 197.7
H2O.entropy_t0   = 188.84
CO2.entropy_t0   = 213.7
CH4.entropy_t0   = 186.25
CH3OH.entropy_t0 = 239.9
N2.entropy_t0    = 191.61
NH3.entropy_t0   = 192.778

if __name__ == '__main__':
  print 'Main'
  print H2.enthalpy_t0

