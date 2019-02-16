import numpy as np
from numpy.linalg import norm
import Astrodynamics as A
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

#################################################################
#                           Problem 1                           #
#################################################################
#Constants
mu = 3.986e5

#Determine properties of hyperbolic orbit
a_h = -3.0e4
e_h = np.sqrt(2)
p_h = a_h*(1 - e_h**2)

Em_h = -mu/(2*a_h)
v_h = np.sqrt(2*(Em_h + mu/-p_h))
nu_h = np.arccos(p_h/(p_h*e_h) - 1/e_h)

v_at_1 = [np.sqrt(mu/p_h)*-np.sin(nu_h), np.sqrt(mu/p_h)*(e_h + np.cos(nu_h)), 0]

if np.dot([0,-p_h,0], v_at_1) < 0:
    nu_h = 2*np.pi - nu_h
    v_at_1 = [np.sqrt(mu / p_h) * -np.sin(nu_h), np.sqrt(mu / p_h) * (e_h + np.cos(nu_h)), 0]

#Flip the coordinate system
v_at_1 = [-1*i for i in v_at_1]

#Determine properties of Elliptical orbit

a_e = 15000
e_e = 0.5
p_e = a_e*(1 - e_e**2)
ra_e = a_e*(1 + e_e)

Em_e = -mu/(2*a_e)
va_e = np.sqrt(2*(Em_e + mu/ra_e))
nu_e = np.arccos(p_e/(ra_e*e_e) - 1/e_e)

v_at_2 = [np.sqrt(mu/p_e)*-np.sin(nu_e), np.sqrt(mu/p_e)*(e_e + np.cos(nu_e)), 0]

TOF = 2*np.pi*np.sqrt(a_e**3/mu)

#Gauss Problem

dt = TOF/2
r_1 = [0, p_h, 0]
r_2 = [-ra_e, 0, 0]

v_1,v_2,direction = A.Lambert_UV(r_1,r_2,dt)

short_index = direction.index('Short')

v_1 = v_1[short_index]
v_2 = v_2[short_index]

#Calculate properties of transfer orbit
[a,e,i,RAAN,w,nu] = A.rv2oe(r_1,v_1)

#Calculate delta Vs
dv_1 = [i - j for i,j in zip(v_1,v_at_1)]
dv1 = norm(dv_1)

dv_2 = [i - j for i,j in zip(v_at_2,v_2)]
dv2 = norm(dv_2)

#################################################################
#                           Problem 2                           #
#################################################################

#Constants
mu = 3.986e5

r1 = 6700
r2 = 42164
delta_i = 28.5

a_0 = (r1 + r2)/2
e_0 = 1 - r1/a_0

ToF_0 = A.ToF(a_0,e_0,180,0)

alpha1_0, dV1_0, dV2_0, dV_0 = A.Split_PC(r1,r2,delta_i)

print('Optimal Split Plane Change Hohmann Transfer')
print('--------------------------------------------')
print('\u0394V = ' + str(round(dV_0,10)) + ' km/s')
print('\u03B11 = ' + str(round(alpha1_0,10)) + '\u00b0')
print('\u03B12 = ' + str(round(delta_i -alpha1_0,10)) + '\u00b0')
print('Time of Flight = ' + str(round(ToF_0/3600,10)) + ' hours')
print('Cost Function value = ' + str(round(dV_0 +2.0*ToF_0/3600,10)))


nu = np.linspace(np.pi/2,np.pi,1000)

ToF_list = []
dV_list = []
ee = []
Cost = []
alpha1_list = []
for TA in nu:

    def f(e):
        func = (r1/(1-e))*(1 - e**2)/(1 + e*np.cos(TA)) - r2
        return func

    e = fsolve(f,0.5)[0]
    ee.append(e)
    a = r1/(1-e)
    p = a*(1 - e**2)

    v_t1 = np.sqrt(2*mu/r1 - mu/a)
    v_t2 = np.sqrt(2*mu/r2 - mu/a)

    h = np.sqrt(mu*p)
    phi_2 = np.degrees(np.arccos(h/(r2*v_t2)))

    alpha1, dV1, dV2, dV = A.Split_PC(r1, r2, delta_i,phi_2,v_t1,v_t2)
    ToF = A.ToF(a,e,np.degrees(TA),0)
    ToF_h = ToF/3600

    # print(TA,dV,ToF_h)

    C = dV + 2.0 * ToF_h
    Cost.append(C)

    dV_list.append(dV)
    ToF_list.append(ToF_h)
    alpha1_list.append(alpha1)

Opt_C = min(Cost)
Opt_index = Cost.index(Opt_C)
Opt_nu = nu[Opt_index]

print('\nMinimal Cost Split Plane Change Transfer')
print('--------------------------------------------')
print('\u0394V = ' + str(round(dV_list[Opt_index],10)) + ' km/s')
print('v1 = ' + str(round(np.degrees(nu[Opt_index]),10)) + '\u00b0')
print('\u03B11 = ' + str(round(alpha1_list[Opt_index],10)) + '\u00b0')
print('\u03B12 = ' + str(round(delta_i -alpha1_list[Opt_index],10)) + '\u00b0')
print('Time of Flight = ' + str(round(ToF_list[Opt_index],10)) + ' hours')
print('Cost Function value = ' + str(round(dV_list[Opt_index] + 2.0*ToF_list[Opt_index],10)))

plt.plot(np.degrees(nu),Cost)

plt.xlim([90,180])
# plt.ylim([-4,4])
plt.xlabel('True anomaly (deg)')
plt.ylabel('Cost Function')
# plt.legend()
plt.grid(True)
plt.title('Cost Function vs True Anomaly')
plt.show()


#################################################################
#                           Problem 4                           #
#################################################################

#Constants
mu = 3.986e5
r_0 = [1131.340, -2282.343, 6672.423]
v_0 = [-5.64305, 4.30333, 2.42879]
dt = 40*60 #sec

r,v = A.Laguerre(r_0,v_0,dt)
