import atom
import molecule
import math
import known_molecules as km
import gas
from scipy.optimize import minimize
from scipy.optimize import fmin
import math


H2O = gas.Gas({km.H2: -1.0, km.O2: -0.5, km.H2O: 1.0})
H2O2 = gas.Gas({km.H2: -1.0, km.O2: -1.0, km.H2O2: 1.0})
#gas.Gas({km.H2O: -1.0, km.O2: -0.5, km.H2O2: 1.0})

CO2 = gas.Gas({km.CO: -1.0, km.O2: -0.5, km.CO2: 1.0})

CH3OH = gas.Gas({km.CO: -1.0, km.H2: -2.0, km.CH3OH: 1.0})
CH3OH_2 = gas.Gas({km.CO2: -1.0, km.H2: -3.0, km.CH3OH: 1.0,km.H2O: 1.0})
#gas.Gas({km.CO2: -1.0, km.H2: -3.0, km.CH3OH: 1.0, km.H2O: 1.0})

CH4 = gas.Gas({km.CO: -1.0, km.H2: -3.0, km.CH4: 1.0, km.H2O: 1.0})
#gas.Gas({km.CO2: -1.0, km.H2: -4.0, km.CH4: 1.0, km.H2O: 2.0})

NO = gas.Gas({km.N2: -0.5, km.O2: -0.5, km.NO: 1.0})
NO2 = gas.Gas({km.N2: -0.5, km.O2: -1.0, km.NO2: 1.0})

NH3 = gas.Gas({km.N2: -0.5, km.H2: -1.5, km.NH3: 1.0})


if __name__ == '__main__':
    print 'Known gas'
    check_gas = True
    for empty_gas in [H2O,H2O2,CO2,CH3OH,CH4,NO,NO2,NH3]:
        if sum([x**2 for x in empty_gas.gas_atom_composition().values()]) != 0.0:
            check_gas *= False
    print 'Are all gases empty: ' + str(bool(check_gas))

    def eq_fun(xvalues):
        P = 1.0
        print xvalues
        [xCO2,xMeOH] = xvalues
        #xMeOH = 0.0
        g0 = gas.Gas({km.CO: 0.00, km.CO2: 0.25, km.H2: 0.75, km.CH3OH: 0.0})
        g =  g0 + xCO2 * CO2 + xMeOH * CH3OH_2
        g.set_pressure(P)
        g.temperature = 750.0
        #print g.partial_pressures.values()
        
        result = g.gibbs()
        print [xCO2,xMeOH,result]
        if True in [x < 0.0 for x in g.partial_pressures.values()]:
            #print [x < 0.0 for x in g.partial_pressures.values()]
            result = 0.0
        return result
    def eq_fun2(xvalues):
        P = 1.0
        print xvalues
        [xMeOH] = xvalues
        g0 = gas.Gas({km.CO: 0.0, km.CO2: 0.25, km.H2: 0.75, km.CH3OH: 0.0})
        g =  g0 + xMeOH * CH3OH_2
        g.set_pressure(P)
        g.temperature = 750.0
        #print g.partial_pressures.values()
        
        result = g.gibbs()
        print [xMeOH,result]
        if True in [x < 0.0 for x in g.partial_pressures.values()]:
            #print [x < 0.0 for x in g.partial_pressures.values()]
            result = 0.0
        return result
    def eq_fun3(xvalues):
        P = 1.0
        T = 273.15+0.0
        #print xvalues
        [xCH4] = xvalues
        g0 = gas.Gas({km.CO: 0.1, km.CH4: 0.0, km.H2: 0.9,km.H2O:0.0})
        g =  g0 + xCH4 * CH4
        g.set_pressure(P)
        g.temperature = T
        #print g.partial_pressures.values()
        
        result = g.gibbs(T)
        #print [xCH4,result]
        if True in [x < 0.0 for x in g.partial_pressures.values()]:
            #print [x < 0.0 for x in g.partial_pressures.values()]
            result = 0.0
        return result

    #res = fmin(eq_fun3, [0.001],xtol=0.00000001 )
    print 'Results'
    #print res
    #print eq_fun2(res.x)
    for i in range(10):
        P = 1.0
        T = 273.15+1500.0
        g0 = gas.Gas({km.CO: 0.1, km.CH4: 0.0, km.H2: 0.9,km.H2O:0.0})
        g =  g0 + 0.1/10*i * CH4
        g.set_pressure(P)
        g.temperature = T
        result = g.gibbs(T)
        print result

    
