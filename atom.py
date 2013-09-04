def elements(z):
    periodic_table = {}
    periodic_table[1]  = {'symbol':'H'}
    periodic_table[2]  = {'symbol':'He'}
    periodic_table[3]  = {'symbol':'Li'}
    periodic_table[4]  = {'symbol':'Be'}
    periodic_table[5]  = {'symbol':'B'}
    periodic_table[6]  = {'symbol':'C'}
    periodic_table[7]  = {'symbol':'N'}
    periodic_table[8]  = {'symbol':'O'}
    periodic_table[9]  = {'symbol':'F'}
    periodic_table[10] = {'symbol':'Ne'}
    periodic_table[11] = {'symbol':'Na'}
    periodic_table[12] = {'symbol':'Mg'}
    return periodic_table[z]

def iso(z):
    list_of_iso = {}
    list_of_iso['H'] = {'H':[1.007825,0.99985],'D':[2.01410178,0.00015]}
    list_of_iso['He'] = {'3He':[3.0160293,0.00000137],'4He':[4.002602,0.99999863]}
    list_of_iso['C'] = {'12C':[12.0,0.9893],'13C':[13.0033548378,0.0107]}
    list_of_iso['O'] = {'16O':[15.9949,0.9976],'17O':[16.9991315,0.00039],'18O':[18,0.00201]} # mass of 18O is not correct
    return list_of_iso[z]
    
class Atom(): # should input be numebr or string for Carbon, 'C' or '6'
    def __init__(self, Z=0, name=''):
        #if not name == '':
        #    self.Z = 
        
        if Z > 0:
            self.Z = Z      

    def __eq__(self, other):
        if other.z == self.z:
            equal = True
        else:
            equal = False
        return equal

    def symbol(self):
        info = elements(self.z)
        return info['symbol']
        
    def iso(self):
        info = iso(elements(self.z)['symbol'])
        return info
        
    def mass(self):
        info = iso(elements(self.z)['symbol'])
        amu = 0
        for i in info.keys():
            amu += (info[i][0]*info[i][1])
        return amu
        

if __name__ == '__main__':
    atom = Atom(1)
    print atom.symbol()
    print atom.mass()
    atom.mass()
