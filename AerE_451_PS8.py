import numpy as np
# from sympy import solveset, Symbol, Eq, S, sqrt, cos, sin
from scipy import sin, cos, sqrt
from scipy.optimize import fsolve

#Problem 2

mu = 3.9860047e5
J2 = 0.0010826
Re = 6378.140
a_dot = 0.9856*np.pi/180/86400

i = np.radians(180 - 63.43)
e = 0.25
dw_dt = 0


def func2(SMA):
    func = 0.9856*np.pi/180/86400 + 3/2 * np.sqrt(mu/SMA**3) * J2 * (Re / (SMA*(1-e**2)))**2*np.cos(np.radians(i))
    return func
# e = 0.25
# i = 116.57
SMA = fsolve(func2, 3700+Re)
SMA = SMA[0] #km

print(SMA)

SMA_new = (SMA-1000)  #meters
print(SMA_new)
RAAN_Rate = -3/2 * np.sqrt(mu/SMA_new**3) * J2 * (Re/(SMA_new*(1-e**2)))**2 * np.cos(np.radians(i)) * 180/np.pi *86400 #deg/day
print(RAAN_Rate)




# def f(a):
#
#     func = a_dot - 3 / 2 * sqrt(mu / a ** 3) * J2 * (Re / (a * (1 - e ** 2))) ** 2 * cos(i)
#
#     return func
#
#
# a = fsolve(f, 2e4)
# a = a[0]
#
# print(a)
#
# a = Symbol('a')
# func = Eq(-3/2 * sqrt(mu/a**3) * J2 * (Re/(a*(1-e**2)))**2 * cos(i), a_dot)
# aa = solveset(func, a, domain=S.Reals)
#
# for i in aa:
#     a2 = float(i) - 1000
#     print(a2)
#
# dOmg_dt = -3/2 * sqrt(mu/a2**3) * J2 * (Re/(a2*(1-e**2)))**2 * cos(i) * (180/np.pi*86400)
# RAAN_Rate = -3/2 * sqrt(mu/a2**3) * J2 * (Re/(a2*(1-e**2)))**2 * cos(i) * 180/np.pi *86400 #deg/day
# print(dOmg_dt)
# print(RAAN_Rate)










