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
    
class Atom():
    def __init__(self, z):
        self.z = z
        
    def symbol(self):
        info = elements(self.z)
        return info['symbol']
        

if __name__ == '__main__':
    atom = Atom(2)
    print atom.symbol()
