#!/bin/env python
# pylint: disable=invalid-name

"""This module contains the Atom class"""

PERIODIC_TABLE = {
    1: {'symbol': 'H'},
    2: {'symbol': 'He'},
    3: {'symbol': 'Li'},
    4: {'symbol': 'Be'},
    5: {'symbol': 'B'},
    6: {'symbol': 'C'},
    7: {'symbol': 'N'},
    8: {'symbol': 'O'},
    9: {'symbol': 'F'},
    10: {'symbol': 'Ne'},
    11: {'symbol': 'Na'},
    12: {'symbol': 'Mg'},
}
LIST_OF_ISO = {
    'H': {'H': [1.007825, 0.99985], 'D': [2.01410178, 0.00015]},
    'He': {'3He': [3.0160293, 0.00000137], '4He': [4.002602, 0.99999863]},
    'C': {'12C': [12.0, 0.9893], '13C': [13.0033548378, 0.0107]},
    'O': {'16O': [15.9949, 0.9976], '17O': [16.9991315, 0.00039],
          '18O': [18, 0.00201]},  # mass of 18O is not correct
}


def elements(z):
    """Returns an element based on its atomic number"""
    return PERIODIC_TABLE[z]


def iso(z):
    """Returns isotopes for an atomic number"""
    return LIST_OF_ISO[z]


class Atom(object):  # should input be number or string for Carbon, 'C' or '6'
    """Class that represents an atom"""

    def __init__(self, Z=0, name=''):
        """Initilize internal variables"""
        self.name = name
        if Z > 0:
            self.Z = Z

    def __eq__(self, other):
        """Returns the equals value"""
        if other.Z == self.Z:
            equal = True
        else:
            equal = False
        return equal

    def symbol(self):
        """Returns the atomic symbol"""
        info = elements(self.Z)
        return info['symbol']

    def iso(self):
        """Returns the atoms isotopes"""
        info = iso(elements(self.Z)['symbol'])
        return info

    def mass(self):
        """Returns the atoms mass"""
        info = iso(elements(self.Z)['symbol'])
        amu = 0
        for i in info.keys():
            amu += (info[i][0] * info[i][1])
        return amu


if __name__ == '__main__':
    atom = Atom(1)
    print atom.symbol()
    print atom.iso()
    print atom.mass()
    atom.mass()
