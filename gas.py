
import molecule
class Gas():
    def __init__(self,molecules,partial_pressure_ratio=None,pressure=None,temperature=None):
        self.molecules = molecules
        self.partial_pressure_ratio = partial_pressure_ratio
        self.p = pressure
        self.temperature = temperature
        if self.partial_pressure_ratio == None:
            print 'ERROR - Partial pressure needed'
        if self.p == None:
            self.p = 1.0 #bar
        if self.temperature == None:
            self.temperature = 300 #K

    def gas_equlibrium(self):
        # calculate the equlibrium gas composition given from the initial gas composition, temperature and pressure
        return self.partial_pressure_ratio
    
    def set_partial_pressure_ratio(self,ppr):
        if sum(ppr) == 1.0:
            self.partial_pressure_ratio = ppr
        else:
            print 'ERROR - sum of pressure ratio should be 1.0'
        return self.partial_pressure_ratio
    

            
    

if __name__ == '__main__':
    print 'Start gas.py'
    import known_molecules as km
    gas_1 = Gas([km.CO,km.CO2,km.H2,km.H2O],[0.0,0.25,0.75,0.0],pressure=2.0,temperature=310)

    print gas_1.gas_equlibrium()
    
