import numpy as np
from scipy.optimize import fsolve


#------------------------------------------------------------------
#                           Problem 2
#------------------------------------------------------------------
mu = 3.9860047e5
Re = 6378.14
J2 = 0.0010826

e = 0.25
i = np.radians(116.57)

def func(a):
    func = 0.9856*(np.pi/180/86400) + 3/2 * np.sqrt(mu/a**3) * J2 * (Re / (a*(1-e**2)))**2*np.cos(i)
    return func

a = fsolve(func, 2*Re)

a = a[0]
a_n = a-1000

dOmg_dt = -3/2 * np.sqrt(mu/a_n**3) * J2 * (Re/(a_n*(1-e**2)))**2 * np.cos(i) * 180/np.pi *86400

ha = a_n * (1+e) - Re
hp = a_n * (1-e) - Re

h_bar = (ha+hp)/2
ecc = (ha-hp) / (ha+hp+2*Re)

dw_dt = -9.9358/(1-ecc**2)**2 * (Re / (Re+h_bar))**3.5 * np.cos(i)

print(dOmg_dt)
print(dw_dt)


#------------------------------------------------------------------
#                           Problem 3
#------------------------------------------------------------------

e = 0.17

a1 = fsolve(func, 2*Re)
a1 = a1[0]

print(a1)

