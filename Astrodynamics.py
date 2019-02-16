#Import Libraries
import numpy as np

def alt():

    # Import libraries
    import numpy as np
    import matplotlib.pyplot as plt

    ###################################
    #            Constants            #
    ###################################
    g0 = 9.81

    ###################################
    #           Calculations          #
    ###################################
    T0 = 288.16
    P0 = 101325
    rho0 = 1.225
    R = 287

    TT = [288.16]
    PP = [101325]
    den = [1.225]
    den2 = [1.225]

    for i in range(1,11001):
        TT.append(T0 + -6.5e-3*i)
        PP.append(P0*(TT[-1]/T0) ** ((-g0/(R*-6.5e-3))))
        den.append(rho0*(TT[-1]/T0) ** -((g0/(R*-6.5e-3)) +1))
        den2.append(PP[i] / (R * TT[i]))

    T0 = TT[-1]; P0 = PP[-1]; rho0 = den[-1]
    for i in range(1,(25000-11000)+1):
        TT.append(TT[-1])
        PP.append(P0 * (np.exp(-(g0 / (R * TT[-1]))))**(i))
        den.append(rho0 * (np.exp(-(g0 / (R * TT[-1])))) ** (i))
        den2.append(PP[i+11000] / (R * TT[i+11000]))

    T0 = TT[-1]; P0 = PP[-1]; rho0 = den[-1]
    for i in range(1,(47000-25000)+1):
        T = TT[-1] + 3e-3
        TT.append(T)
        PP.append(P0 * (TT[-1] / T0) ** ((-g0 / (R * 3e-3))))
        den.append(rho0 * (TT[-1] / T0) ** -((g0 / (R * 3e-3)) + 1))
        den2.append(PP[i+25000] / (R * TT[i+25000]))

    T0 = TT[-1]; P0 = PP[-1]; rho0 = den[-1]
    for i in range(1,(53000-47000)+1):
        TT.append(TT[-1])
        PP.append(P0 * np.exp((-g0 / (R * TT[-1])))**(i))
        den.append(rho0 * (np.exp(-(g0 / (R * TT[-1])))) ** (i))
        den2.append(PP[i+47000] / (R * TT[i+47000]))

    T0 = TT[-1]; P0 = PP[-1]; rho0 = den[-1]
    for i in range(1,(79000-53000)+1):
        T = TT[-1] + -4.5e-3
        TT.append(T)
        PP.append(P0 * (TT[-1] / T0) ** ((-g0 / (R * -4.5e-3))))
        den.append(rho0 * (TT[-1] / T0) ** -((g0 / (R * -4.5e-3)) + 1))
        den2.append(PP[i+53000] / (R * TT[i+53000]))

    T0 = TT[-1]; P0 = PP[-1]; rho0 = den[-1]
    for i in range(1,(90000-79000)+1):
        TT.append(TT[-1])
        PP.append(P0 * np.exp((-g0 / (R * TT[-1])))**(i))
        den.append(rho0 * (np.exp(-(g0 / (R * TT[-1])))) ** (i))
        den2.append(PP[i+79000] / (R * TT[i+79000]))

    T0 = TT[-1]; P0 = PP[-1]; rho0 = den[-1]
    for i in range(1,(105000-90000)+1):
        T = TT[-1] + 4e-3
        TT.append(T)
        PP.append(P0 * (TT[-1] / T0) ** ((-g0 / (R * 4e-3))))
        den.append(rho0 * (TT[-1] / T0) ** ((-g0 / (R * 4e-3)) + 1))
        den2.append(PP[i+90000] / (R * TT[i+90000]))


    # import numpy as np
    # import matplotlib.pyplot as plt


    # h = list(np.linspace(0,105000,105001))
    #
    # plt.plot(TT, h)
    #
    # #plt.figure()
    # plt.xlim([160,320])
    # plt.ylim([0,110000])
    # plt.xlabel('Temperature (K)')
    # plt.ylabel('Altitude (m)')
    # # plt.legend()
    # plt.grid(True)
    # plt.title('Standard Atmosphere')
    # plt.show()
    #
    #
    # plt.plot(PP, h)
    #
    # # plt.figure()
    # plt.xlim([-5000,101325])
    # plt.ylim([0,110000])
    # plt.xlabel('Pressure (Pa)')
    # plt.ylabel('Altitude (m)')
    # # plt.legend()
    # plt.grid(True)
    # plt.title('Standard Atmosphere')
    # plt.show()
    #
    # DEN = np.zeros(len(h))
    # for i in range(len(h)):
    #     DEN[i] = 1.225*np.exp(-h[i]/8000)

    # plt.plot(den, h)
    #
    # # plt.figure()
    # plt.xlim([-.05,1.225])
    # plt.ylim([0,110000])
    # plt.xlabel('Density (kg/m^3)')
    # plt.ylabel('Altitude (m)')
    # # plt.legend()
    # plt.grid(True)
    # plt.title('Standard Atmosphere')
    # plt.show()

    return TT, PP, den

def Laguerre(r0,v0,dt):

    # Import libraries
    import numpy as np
    from numpy.linalg import norm
    from numpy import sqrt
    import math

    #Constants
    if norm(r0) > 500:
        mu = 3.986012e5
    else:
        mu = 1.0

    #Calculations
    r_mag = norm(r0)
    v_mag = norm(v0)
    h0 = np.cross(r0,v0)
    h_mag = norm(h0)
    energy = (v_mag**2)/2 - mu/r_mag
    a = -mu/(2*energy)
    alpha = -(2*energy)/mu
    p = h_mag**2/mu
    e = sqrt(1 - (p*alpha))

    nu0 = np.arccos(round( ((p/r_mag) - 1) / e, 12)) #; print(nu0)

    t = dt

    Period = 2*np.pi*sqrt(a**3 / mu)

    if dt > 0:
        sign = 1
    else:
        sign = -1

    if alpha < 0:
        dt = dt - sign * math.floor(abs(dt)/Period) * p

    #Loop Setup
    x0 = sqrt(mu)* dt*abs(alpha)
    x = x0
    k = 0
    error = 1
    tol = 1e-12

    data = []
    xx = []
    CC = []
    SS = []

    while (abs(error) > tol) and (k < 1000):
        k += 1

        z = x ** 2 * alpha

        SC = Stumpff(z, 20)
        S = float(SC[0])
        C = float(SC[1])

        F = (1 - (r_mag*alpha))*S*x**3 + np.dot(r0,v0)/sqrt(mu)*C*x**2 + r_mag*x - sqrt(mu)*t
        Fp = C*x**2 + np.dot(r0,v0)/sqrt(mu)*(1 - S*z)*x + r_mag*(1 - C*z)
        Fpp = (1 - (r_mag*alpha))*(1 - S*z)*x + np.dot(r0,v0)/sqrt(mu)*(1-C*z)

        delta = 2*sqrt(4*Fp**2 - 5*F*Fpp)

        if Fp > 0:
            sign = 1
        else:
            sign = -1

        dx = 5*F / (Fp + sign*delta)
        error = dx**2*alpha

        x = x - dx


        data.append(dx)
        xx.append(x)
        CC.append(C)
        SS.append(S)

    #print('k = ' + str(k))

    # pulling data from loop
    x = xx[-1]
    C = CC[-1]
    S = SS[-1]
    z = x ** 2 / a

    # Calculations
    f = 1 - (x**2/r_mag)*C
    g = t - (x ** 3) * S / sqrt(mu)

    r_a = [i * f for i in r0]
    r_b = [i * g for i in v0]

    r = [i + j for i, j in zip(r_a, r_b)]

    r_new = norm(r)

    f_dot = (sqrt(mu) / (r_new * r_mag)) * x * (z * S - 1)
    g_dot = 1 - (x ** 2) * C / r_new

    v_a = [i * f_dot for i in r0]
    v_b = [i * g_dot for i in v0]

    v = [i + j for i, j in zip(v_a, v_b)]
    v_new = norm(v)

    #import pdb; pdb.set_trace()

    # print('\n'+ 'New r = ' + str(r))

    # print('\n'+ 'New v = ' + str(v))

    # print('\n' + str(f*g_dot - f_dot*g) + '  =  1?')

    return r,v
def vector_UV(r0,v0,del_t):
    # Import libraries
    import numpy as np
    from math import log
    from numpy.linalg import norm
    from math import sqrt

    #Constants
    if norm(r0) > 500:
        mu = 3.986012e5
    else:
        mu = 1.0

    #Calculations
    r_mag = norm(r0)
    v_mag = norm(v0)
    h0 = np.cross(r0,v0)
    h_mag = norm(h0)
    energy = (v_mag**2)/2 - mu/(r_mag)
    a = -mu/(2*energy)
    p = h_mag**2/mu
    e = sqrt(1 - (p/a))

    nu0 = np.arccos( ((p/r_mag) - 1) / e ); #print(nu0)

    t = del_t

    # Approximation for x1
    if a > 0:
        x1 = sqrt(mu) * (del_t) / a
    else:
        if del_t>0:
            sign = 1

        else:
            sign = -1

        num = -2 * mu * (del_t)
        den = a * (np.dot(r0, v0) + sign * sqrt(-mu*a) * (1 - (r_mag/a)))

        x1 = sign * sqrt(-a)*log(num/den)

    x = x1
    dt = 1
    i = 0
    k = 0

    data = []
    xx = []
    CC = []
    SS = []

    while (np.abs(dt) > 1e-12) and k<1000:
        k += 1

        z = x ** 2 / a

        SC = Stumpff(z, 20)   #; print(SC); print('\n')
        S = SC[0]
        C = SC[1]

        t1 = ((np.dot(r0, v0) / sqrt(mu)) * x ** 2 * C + (1 - (r_mag / a)) * x ** 3 * S + r_mag * x) / sqrt(mu)
        dt_dx = (x ** 2 * C + (np.dot(r0, v0) / sqrt(mu)) * x * (1 - z * S) + r_mag * (1 - z * C)) / sqrt(mu)

        x2 = x + (t - t1) / dt_dx

        dt = t - t1

        x = x2

        #print(dt)

        data.append(dt)
        xx.append(x)
        CC.append(C)
        SS.append(S)

    # print('k = ' + str(k))

    # pulling data from loop
    x = xx[-1]
    C = CC[-1]
    S = SS[-1]
    z = x ** 2 / a

    #import pdb; pdb.set_trace()

    # Calculations
    f = 1 - (x**2/r_mag)*C
    g = t - (x ** 3) * S / sqrt(mu)

    r_a = [i * f for i in r0]
    r_b = [i * g for i in v0]

    r = [i + j for i, j in zip(r_a, r_b)]

    r_new = norm(r)

    f_dot = (sqrt(mu) / (r_new * r_mag)) * x * (z * S - 1)
    g_dot = 1 - (x ** 2) * C / r_new

    v_a = [i * f_dot for i in r0]
    v_b = [i * g_dot for i in v0]

    v = [i + j for i, j in zip(v_a, v_b)]
    v_new = norm(v)

    # print('\n'+ 'New r = ' + str(r))

    # print('\n'+ 'New v = ' + str(v))

    # print('\n' + str(f*g_dot - f_dot*g) + '  =  1?')

    return r,v
def scalar_UV(r_mag, v_mag, nu0_rad, del_t):
    # Import libraries
    import numpy as np
    from math import log,sqrt
    from sympy import Symbol, Eq, solve
    from numpy.linalg import norm

    if r_mag > 500:
        mu = 3.986012e5
    else:
        mu = 1.0

    nu0 = nu0_rad

    # Calculations
    energy = (v_mag ** 2) / 2 - mu / (r_mag)
    a = - mu / (2 * energy)  # print(a)

    rp = r_mag * np.cos(nu0)  # print(rp)
    rq = r_mag * np.sin(nu0)  # print(rq)
    r0 = [rp, rq, 0]
    # print('r0 = ' + str(r0))

    e = Symbol('e')
    f_e = Eq(r_mag * (1 + e * np.cos(nu0)), a * (1 - e ** 2))
    ee = solve(f_e)

    if type(ee) is list:
        e = [i for i in ee if i > 0]
        e = float(e[0])
    # print('e = ' + str(e))

    p = a * (1 - e ** 2)
    coeff = sqrt(mu / p)

    vp = - coeff * np.sin(nu0)
    vq = coeff * (e + np.cos(nu0))
    v0 = [vp, vq, 0]
    # print('v0 = ' + str(v0) + '\n')

    t = del_t

    r,v = vector_UV(r0,v0,del_t)

    return r,v
def Lambert_UV(r_1, r_2, t0):
    import numpy as np
    from numpy import dot, sin, cos, sqrt
    from numpy.linalg import norm

    # Constants
    if norm(r_1) > 500:
        mu = 3.986012e5
    else:
        mu = 1.0

    # r_1 = np.array(r_1)
    # r_2 = np.array(r_2)

    r1 = norm(r_1)
    r2 = norm(r_2)
    dnu_1 = np.arccos(dot(r_1, r_2) / (r1 * r2))
    dnu_2 = (2 * np.pi - dnu_1)

    nu = [dnu_1, dnu_2]

    v_1 = [0, 0]
    v_2 = [0, 0]
    k = 0;
    n = 0
    for j in nu:

        tol = 1e-7
        i_max = 50;
        i = 0
        dt = 1
        z = 0.001

        while abs(dt) > tol and i <= i_max:
            A = (sqrt(r1 * r2) * sin(j)) / sqrt(1 - cos(j))

            S, C = Stumpff(z, 10)

            y = r1 + r2 - A * (1 - z * S) / sqrt(C)
            x = sqrt(y / C)

            t = (x ** 3 * S + A * sqrt(y)) / sqrt(mu)

            dt = t0 - t

            dSdz = 1 / (2 * z) * (C - 3 * S)
            dCdz = 1 / (2 * z) * (1 - z * S - 2 * C)

            dtdz = (x ** 3 * (dSdz - (3 * S * dCdz) / (2 * C)) + A / 8 * ((3 * S * sqrt(y)) / C + A / x)) / sqrt(mu)

            z = z + dt / dtdz
            # print(i)
            # print(str(S)+'\t'+str(C)+'\t'+str(x)+'\t'+str(y)+'\t'+str(z)+'\t'+str(t)+'\n')
            i += 1

        f = 1 - y / r1
        g = A * sqrt(y / mu)
        g_dot = 1 - y / r2

        fr1 = [f * k for k in r_1]
        v1 = [(i - j) / g for i, j in zip(r_2, fr1)]

        gdr2 = [g_dot * k for k in r_2]
        v2 = [(i - j) / g for i, j in zip(gdr2, r_1)]

        v_1[k] = v1
        v_2[k] = v2
        k += 1

        label = ['Short', 'Long']
        # print(label[n])
        # print(v_1[n])
        # print(v_2[n])
        # n += 1

    return v_1, v_2, label

def v_to_E(v, e):
    import numpy as np

    if v == np.pi:
        E = np.pi
    else:
        q = np.sqrt((1 - e) / (1 + e))
        E = 2 * np.arctan(q * np.tan(v / 2))

    return E
def v_to_F(v, e):
    import numpy as np

    if v == np.pi:
        F = np.pi
    else:
        # F = np.arccosh((e + np.cos(v)) / (1 + e * np.cos(v)))
        q = np.sqrt((e - 1) / (e + 1))
        F = 2 * np.arctanh(q * np.tan(v/2))

    return F

def Stumpff(z, n_terms):

    def factorial(n):
        if n <= 0:
            return 1
        else:
            return n * factorial(n - 1)

    S = 0
    C = 0

    for i in range(0,n_terms):

        S = S + ((-z)**i/factorial(2*i + 3))
        C = C + ((-z)**i/factorial(2*i + 2))
        #print(S); print(C); print(i); print('\n')
    return S, C

def rv2oe(r_2,v_2):
    import numpy as np
    from numpy.linalg import norm
    import math

    if norm(r_2) <= 500:
        mu = 1
    else:
        mu = 3.986012e5

    r2 = norm(r_2)
    v2 = norm(v_2)

    h0 = np.cross(r_2, v_2)
    h_mag = norm(h0)
    energy = (v2**2)/2 - mu/r2
    a = - mu / (2 * energy)
    p = h_mag**2 / mu
    e = np.sqrt(1 - (p / a))

    vva = [(v2**2 - mu/r2)/mu*i for i in r_2]
    rdvv = [np.dot(r_2,v_2)/mu*i for i in v_2]

    e_vec = [i-j for i,j in zip(vva,rdvv)]#; print(e_vec)

    n = np.cross([0,0,1],h0)

    i = math.degrees(np.arccos(np.dot([0,0,1],h0)/h_mag))

    w = math.degrees(np.arccos(np.dot(n, e_vec) / (e * norm(n))))

    RAAN = math.degrees(np.arccos(np.dot([1, 0, 0], n) / norm(n)))

    nu2 = math.degrees(np.arccos(np.dot(e_vec, r_2) / (e * r2)))

    if n[1] < 0:
        RAAN = 360 - math.degrees(np.arccos(np.dot([1,0,0],n)/norm(n)))

    if e_vec[2] < 0:
        w = 360 - math.degrees(np.arccos(np.dot(n,e_vec)/(e*norm(n))))

    if np.dot(r_2,v_2) < 0:
        nu2 = 360 - nu2


    # print('\na = '+str(round(a,5)) + '\n')
    # print('e = '+str(round(e,5))+'\n')
    # print('i = '+str(round(i,5))+u'\xb0 \n')
    # print('\u03A9 = ' + str(round(RAAN, 5)) + u'\xb0 \n')
    # print('\u03C9 = '+str(round(w,5))+u'\xb0 \n')
    # print('\u03BD = '+str(round(nu2,5))+u'\xb0 \n')
    print('State vector = ' + str([a,e,i,RAAN,w,nu2]))

    return a,e,i,RAAN,w,nu2

def ECI(H,Lat,LST):

    import numpy as np

    a_e = 6378.15  # km
    f = 1/298.26  #0.003353

    x = (a_e/(np.sqrt(1 - (2*f - f**2)*(np.sin(Lat))**2)) + H)*np.cos(Lat)
    z = (((a_e*(1-f)**2)/np.sqrt(1 - (2*f - f**2)*(np.sin(Lat))**2)) + H)*np.sin(Lat)

    R = [x * np.cos(LST), x * np.sin(LST), z]

    return R

def TimeToUTC(year, month, day, hour, min, sec, UTC_offset):
    import math
    """Takes a given date and time and adjusts it to UTC time, making corrections to the date if necessary."""
    #Check that the UTC Offset does not exceed 24 hours, this is non-sensical
    if abs(UTC_offset) > 24:
        print("Error: UTC offset specification is greater than 24 hours. Please re-adjust to a sensical value. "
              "Exiting...")
        return None
    def month_limit(month, year = 0):
        import math
        if month == 2 and year == 0:
            print("ERROR! If the month is February, the year must also be input. Exiting...")
            return None
        #Determine the day limit of the month based on the input year and month
        if month == 1 or month == 3 or month == 5 or month == 7 or month == 8 or month == 10 or month == 12:
            day_limit = 31
        elif month == 4 or month == 6 or month == 9 or month == 11:
            day_limit = 30
        elif month == 2:
            remainder = year % 4
            if remainder == 0:
                day_limit = 29
            else:
                day_limit = 28
        return day_limit

    #Check if UTC offset input is an integer or some fractional hour
    frac_check = UTC_offset - math.floor(UTC_offset)
    if frac_check != 0: #Hour fraction! Must add some minutes
        min_add = frac_check * 60 #convert from hours to min
        min = min + math.floor(min_add)
        frac_check = min_add - math.floor(min_add)
        if frac_check != 0: #Minute fraction! Must add some seconds
            sec_add = frac_check * 60 #convert from min to seconds
            sec = sec + sec_add

    hour = hour + math.floor(UTC_offset)
    #If the hour exceeds the bounds for a day, adjust the day and hour numbers
    if hour > 24:
        day = day + 1
        hour = hour % 24
    if hour < 0:
        day = day - 1
        hour = 24 + hour
    #Check Day Bounds for the month
    day_limit = month_limit(month, year)
    if day < 1:
        month = month - 1
        if month < 1:
            month = 12
            year = year - 1
            day = month_limit(month, year)
        else:
            day = month_limit(month, year)
    elif day > day_limit:
        month = month + 1
        if month > 12:
            month = 1
            year = year + 1
        day = 1

    return year, month, day, hour, min, sec
def JulianDate(year,month,day,hour,min,sec,UTC_offset=0):
    import math
    """Calculates the Julian Date of a given epoch."""
    #Check that the inputs fit within the valid bounds
    if year >= 2099 or year <=1901:
        print("Error: Year exceeds acceptable range of 1901-2099. Exiting...")
        return None
    if month < 1 or month > 12:
        print("Error: Month falls outside of acceptable range 1-12. Exiting...")
        return None
    if day < 1 or day > 31:
        print("Error: Day falls outside of acceptable range 1-31. Exiting...")
        return None

    #Adjust the input time to UTC time. Capture any potential day changes that can occur.
    if UTC_offset != 0:
        year, month, day, hour, min, sec = TimeToUTC(year, month, day, hour, min, sec, UTC_offset)
    J0 = 367*year - math.floor(7/4*(year+math.floor((month+9)/12)))+math.floor(275*month/9)+day+1721013.5
    DayFrac = (hour + min/60 + sec/3600) / 24

    JulianDate = J0 + DayFrac

    return JulianDate, J0

def Sidereal(Lat,year,month,day,hour,min,sec,UTC_offset=0):

    JD,J0 = JulianDate(year,month,day,hour,min,sec,UTC_offset)
    JD2000 = 2451545.0

    DayFrac = (hour + min / 60 + sec / 3600) / 24

    T0 = (J0 - JD2000)/36525

    Theta_g0 = 100.4606184 + 36000.77004*T0 + 0.0003879333*T0**2 - 2.583e-8*T0**3
    Theta_g0 = np.mod(Theta_g0,360)

    Theta = Theta_g0 + 360.98564724*DayFrac

    LST = Theta + Lat
    LST = np.mod(LST,360)

    return LST,JD

def Gauss_Imp(Lat_deg,H,time_array,RA_array, DEC_array, LST_array):
    # Import libraries
    import numpy as np
    from numpy import cross,dot
    from numpy.linalg import norm
    from sympy import solveset, Symbol, Eq, S
    import matplotlib.pyplot as plt

    ###################################
    #            Inputs               #
    ###################################
    Lat = np.radians(Lat_deg)
    time = time_array
    RA = RA_array
    DEC = DEC_array
    LST = LST_array

    mu = 3.9860e5

    tau1 = time[0] - time[1]
    tau3 = time[2] - time[1]
    tau = time[2] - time[0]

    ###################################
    #             R vectors           #
    ###################################

    R_list = []
    for i in range(0,len(LST)):
        R_list.append(ECI(H,Lat,LST[i]))

    R1 = np.array(R_list[0]); R2 = np.array(R_list[1]); R3 = np.array(R_list[2])

    ###################################
    #         Rho_hat vectors         #
    ###################################

    rho_list = []
    for i in range(0,len(LST)):
        Lx = np.cos(RA[i])*np.cos(DEC[i])
        Ly = np.sin(RA[i])*np.cos(DEC[i])
        Lz = np.sin(DEC[i])
        #print('\n' + str([Lx,Ly,Lz]) + '\n')

        rho_list.append([Lx,Ly,Lz])

    rho1 = np.array(rho_list[0]); rho2 = np.array(rho_list[1]); rho3 = np.array(rho_list[2])

    p1 = cross(rho2,rho3)
    p2 = cross(rho1,rho3)
    p3 = cross(rho1,rho2)

    ###################################
    #             D Matrix            #
    ###################################

    D0 = dot(rho1,p1)

    D11 = dot(R1,p1)
    D12 = dot(R1,p2)
    D13 = dot(R1,p3)

    D21 = dot(R2,p1)
    D22 = dot(R2,p2)
    D23 = dot(R2,p3)

    D31 = dot(R3,p1)
    D32 = dot(R3,p2)
    D33 = dot(R3,p3)

    ###################################
    #           r_2 and v_2           #
    ###################################

    A = (-D12*(tau3/tau) + D22 + D32*(tau1/tau))/D0
    B = (D12*(tau3**2 - tau**2)*(tau3/tau) + D32*(tau**2 - tau1**2)*(tau1/tau))/(6*D0)
    E = dot(R2,rho2)
    R2_mag2 = norm(R2)**2

    a = -(A**2 + 2*A*E + R2_mag2)
    b = -2*mu*B*(A+E)
    c = -(mu*B)**2

    x = Symbol('x')
    func = Eq(x**8 + a*x**6 + b*x**3 + c, 0)
    xx = solveset(func,x,domain=S.Reals);

    R_2 = 0;
    for i in xx:
        if i > 6378:
            R_2 = float(i)
            #print(R_2)

    rho1_mag = (((6*(D31*tau1/tau3 + D21*tau/tau3)*R_2**3 + mu*D31*(tau**2 - tau1**2)*tau1/tau3))/(6*R_2**3 + mu*(tau**2 - tau3**2)) - D11)/D0 # print(rho1_mag)
    rho2_mag = A + mu*B/R_2**3 #; print(rho2_mag)
    rho3_mag = (((6*(D13*tau3/tau1 - D23*tau/tau1)*R_2**3 + mu*D13*(tau**2 - tau3**2)*tau3/tau1))/(6*R_2**3 + mu*(tau**2 - tau1**2)) - D33)/D0 #; print(rho3_mag)

    c1 = tau3/tau * (1 + (mu/(6*R_2**3)*(tau**2 - tau3**2))) #C1
    c3 = -tau1/tau * (1 + (mu/(6*R_2**3)*(tau**2 - tau1**2))) #C3

    f1 = 1 - mu*tau1**2/(2*R_2**3)
    f3 = 1 - mu*tau3**2/(2*R_2**3)

    g1 = tau1 - mu*tau1**3/(6*R_2**3)
    g3 = tau3 - mu*tau3**3/(6*R_2**3)

    r_1 = R1 + rho1_mag*rho1 ; print('r_1 = ' + str(r_1))
    r_2 = R2 + rho2_mag*rho2 ; print('r_2 = ' + str(r_2))
    r_3 = R3 + rho3_mag*rho3 ; print('r_3 = ' + str(r_3))

    v_2 = (-f3*r_1 + f1*r_3)/(f1*g3 - f3*g1) ; print('v_2 = ' + str(v_2))

    ###################################
    #           Improvement           #
    ###################################
    r_list = [norm(r_2)]
    v_list = [norm(v_2)]
    c1_list = [c1]
    c3_list = [c3]
    rho1_list = [rho1_mag]
    rho2_list = [rho2_mag]
    rho3_list = [rho3_mag]
    d_rho = 1
    k = 0
    while d_rho > 1e-5:
        k+=1; print('\nN = ' + str(k))

        r2 = norm(r_2)
        v2 = norm(v_2)

        alpha = 2/r2 - v2**2/mu
        X = Symbol('X')
        SC = Stumpff(alpha*X**2,5); SS = SC[0]; CC = SC[1]

        ###################################
        #            X1 and X3            #
        ###################################

        func2 = Eq(np.sqrt(mu)*tau1,(dot(r_2,v_2)*X**2*CC)/np.sqrt(mu) + (1 - alpha*r2)*X**3*SS + r2*X)
        XX = solveset(func2,X,domain=S.Reals)
        for i in XX:
            X1 = i
            #print(X1)

        func3 = Eq(np.sqrt(mu)*tau3,(dot(r_2,v_2)*X**2*CC)/np.sqrt(mu) + (1 - alpha*r2)*X**3*SS + r2*X)
        XX = solveset(func3,X,domain=S.Reals);
        for i in XX:
            X3 = i
            #print(X3)

        SC1 = Stumpff(alpha*X1**2,5); S1 = SC1[0]; C1 = SC1[1]
        SC3 = Stumpff(alpha*X3**2,5); S3 = SC3[0]; C3 = SC3[1]

        f1n = 1 - (X1**2*C1)/r2 #; print(f1n)
        f3n = 1 - (X3**2*C3)/r2 #; print(f3n)
        g1n = tau1 - (X1**3*S1)/np.sqrt(mu) #; print(g1n)
        g3n = tau3 - (X3**3*S3)/np.sqrt(mu)# ; print(g3n)

        f1_avg = (f1+f1n)/2
        f3_avg = (f3+f3n)/2
        g1_avg = (g1+g1n)/2
        g3_avg = (g3+g3n)/2

        c1n = g3_avg/(f1_avg*g3_avg - f3_avg*g1_avg); c1_list.append(c1n)
        c3n = -g1_avg/(f1_avg*g3_avg - f3_avg*g1_avg); c3_list.append(c3n)

        rho1n = (-D11 + D21/c1n - c3n*D31/c1n)/D0; rho1_list.append(rho1n)
        rho2n = (-c1n*D12 + D22 - c3n*D32)/D0; rho2_list.append(rho2n) #; print(rho2n)
        rho3n = (-c1n*D13/c3n + D23/c3n - D33)/D0; rho3_list.append(rho3n)

        r_1n = R1 + rho1n*rho1 #; print(r_1n)
        r_2n = R2 + rho2n*rho2 ; print('r2 = ' + str(r_2n))
        r_3n = R3 + rho3n*rho3 #; print(r_3n)

        v_2n = (-f3_avg*r_1n + f1_avg*r_3n)/(f1_avg*g3_avg - f3_avg*g1_avg); print('v2 = ' + str(v_2n))

        d_rho = abs(rho2n - rho2_list[k-1])
        d_c1 = abs(c1n - c1_list[k-1])
        d_c3 = abs(c3n - c3_list[k-1])

        r_2 = np.array(r_2n).astype(np.float64); r_list.append(norm(r_2))
        v_2 = np.array(v_2n).astype(np.float64); v_list.append(norm(v_2))

        #print(X1)
        #print(X3)
        print('r2 = ' + str(norm(r_2)))
        print('v2 = ' + str(norm(v_2)))
        print('\u03c12 = ' + str(rho2n))
        print('\u0394\u03c1 = ' + str(d_rho))
        print('\u0394c1 = ' + str(d_c1))
        print('\u0394c3 = ' + str(d_c3))


    print('----------------------------------------------------------------------')
    print('\nr_2 = ' + str(r_2))
    print('v_2 = ' + str(v_2))

    print('\nr2 = ' + str(r2))
    print('v2 = ' + str(v2))
    print('\n----------------------------------------------------------------------\n\n')


    ###################################
    #         Orbital Elements        #
    ###################################

    rv2oe(r_2,v_2)

    ###################################
    #         Plot Convergence        #
    ###################################

    iter = []
    for i in range(0,len(rho2_list)):
        iter.append(i)
    #print(iter)

    for i in iter:
        plt.plot(iter,rho2_list)

    plt.xlim([0,k])
    plt.ylim([570,572])
    plt.xlabel('Iteration')
    plt.ylabel('Magnitude of \u03c12')
    #plt.legend()
    plt.grid(True)
    plt.title('Convergence of \u03c12')
    plt.show()

def Gibbs_Method(r_1,r_2,r_3):
    import numpy as np
    import math
    from numpy.linalg import norm

    # r_1 = [1.41422511, 0, 1.414202]
    # r_2 = [1.81065659, 1.06066883, 0.3106515]
    # r_3 = [1.35353995, 1.41422511, -0.6464495]

    if norm(r_1) <= 500 and norm(r_2) <= 500:
        mu = 1
    else:
        mu = 3.986012e5

    r1 = norm(r_1)
    r2 = norm(r_2)
    r3 = norm(r_3)

    #print(np.dot(r_1, np.cross(r_2, r_3)))
    #print(mu)

    C23 = np.array(np.cross(r_2, r_3))/(r2*r3) #; print(C23)
    n = np.array(r_1)/r1 #; print(n)

    # print(np.dot(C23, n))

    if np.dot(C23,n) < 1e-5:

        N = r1*np.array(np.cross(r_2,r_3)) + r2*np.array(np.cross(r_3,r_1)) + r3*np.array(np.cross(r_1,r_2))
        D = np.array(np.cross(r_1,r_2)) + np.array(np.cross(r_2,r_3)) + np.array(np.cross(r_3,r_1))
        S = np.array(r_1)*(r2-r3) + np.array(r_2)*(r3-r1) + np.array(r_3)*(r1-r2)

        N_mag = norm(N)
        D_mag = norm(D)

        v_2 = np.sqrt(mu/(N_mag*D_mag)) * ((np.array(np.cross(D,r_2))/r2) + S); #print(v_2)
        v2 = norm(v_2); #print(v2)

    else:
        print('r1, r2, and r3 are not coplanar vectors')
        exit()

    rv2oe(r_2,v_2)
def Gibbs(Lat_deg, H, Az_array_rad, El_array_rad, rho_array, LST_array_rad):
    # Import libraries
    import numpy as np

    #Givens
    mu = 3.986012e5

    Lat = np.radians(Lat_deg)
    LST = LST_array_rad
    Az = Az_array_rad
    El = El_array_rad
    rho = rho_array

    R_list = []
    for i in range(0, len(LST)):
        R_list.append(ECI(H, Lat, LST[i]))

    # Right ascension and declination

    dec = []
    h = []
    alpha = []
    for i in range(0, len(LST)):

        dec.append(np.arcsin(np.cos(Lat) * np.cos(Az[i]) * np.cos(El[i]) + np.sin(Lat) * np.sin(El[i])))

        if 0 < Az[i] < np.pi:
            h.append(2 * np.pi - np.arccos(
                (np.cos(Lat) * np.sin(El[i]) - np.sin(Lat) * np.cos(Az[i]) * np.cos(El[i])) / np.cos(dec[i])))
        else:
            h.append(
                np.arccos((np.cos(Lat) * np.sin(El[i]) - np.sin(Lat) * np.cos(Az[i]) * np.cos(El[i])) / np.cos(dec[i])))

        alpha.append(LST[i] - h[i])

    # Observation point position vectors
    r_list = [];
    for i in range(0, len(LST)):
        Lx = np.cos(alpha[i]) * np.cos(dec[i])
        Ly = np.sin(alpha[i]) * np.cos(dec[i])
        Lz = np.sin(dec[i])

        # print('\n' + str([Lx,Ly,Lz]) + '\n')

        r_list.append([k * rho[i] for k in [Lx, Ly, Lz]])

    # print(r_list)
    rho_1 = r_list[0];
    rho_2 = r_list[1];
    rho_3 = r_list[2]

    # print(R_list)
    R_1 = R_list[0];
    R_2 = R_list[1];
    R_3 = R_list[2]

    r_1 = [i + j for i, j in zip(R_1, rho_1)]
    r_2 = [i + j for i, j in zip(R_2, rho_2)]
    r_3 = [i + j for i, j in zip(R_3, rho_3)]

    Gibbs_Method(r_1, r_2, r_3)

def Split_PC(r1,r2,delta_i_deg,phi_deg=0,v1=0,v2=0):
    '''Determines the optimal weight of the first and second plane change impulses for a split plane change maneuver.
       If the transfer is a Hohmann transfer, phi and nu are left as 0. If not, phi is the flight path angle in the
       transfer orbit at r2 and v1 and v2 are the velocities in the transfer orbit at perigee and r2, respectively'''
    #Import Libraries
    import numpy as np
    from numpy.linalg import norm
    from scipy import sin,cos,sqrt
    from scipy.optimize import fsolve

    #Constants
    delta_i = np.radians(delta_i_deg)
    phi = np.radians(phi_deg)

    if norm(r1) > 500:
        mu = 3.986e5
    else:
        mu = 1.0

    if phi == 0:
        a_t = (r1 + r2)/2
        V_ta = sqrt(2 * mu / r2 - mu / a_t)
        V_tp = sqrt(2 * mu / r1 - mu / a_t)
    else:
        V_tp = v1
        V_ta = v2

    V_c1 = sqrt(mu/r1)
    V_c2 = sqrt(mu/r2)


    def f(alpha1):
        dV1 = sqrt(V_tp**2 + V_c1**2 - 2*V_tp*V_c1*cos(alpha1))
        dV2 = sqrt(V_c2**2 + (V_ta*cos(phi))**2 - 2*V_c2*V_ta*cos(phi)*cos(delta_i-alpha1))

        func = V_tp*V_c1*sin(alpha1)/dV1 - V_c2*V_ta*cos(phi)*sin(delta_i-alpha1)/dV2

        return func

    alpha1 = fsolve(f,0)
    alpha1 = alpha1[0]

    dV1 = sqrt(V_tp ** 2 + V_c1 ** 2 - 2 * V_tp * V_c1 * cos(alpha1))
    dV2 = sqrt(V_c2 ** 2 + (V_ta * cos(phi)) ** 2 - 2 * V_c2 * V_ta * cos(phi) * cos(delta_i - alpha1))
    dV_T = dV1+dV2

    return alpha1*180/np.pi,dV1,dV2,dV_T

def ToF(a,e,nu1_deg,nu0_deg=0):
    if abs(a) > 50:
        mu = 3.986e5
    else:
        mu = 1.0

    nu1 = np.radians(nu1_deg)
    nu0 = np.radians(nu0_deg)

    if e <1:
        E1 = v_to_E(nu1,e)
        E0 = v_to_E(nu0,e)
        ToF = np.sqrt(a**3/mu)*( (E1 - e*np.sin(E1)) - (E0 - e*np.sin(E0)) )
    elif e ==1:
        p = a * (1 - e ** 2)
        D1 = np.sqrt(p)*np.tan(nu1/2)
        D0 = np.sqrt(p)*np.tan(nu0/2)
        ToF = 1/(2*np.sqrt(mu)) * ( (p*D1 + (1/3)*D1**3) - (p*D0 + (1/3)*D0**3) )
    elif e > 1:
        F1 = v_to_F(nu1,e)
        F0 = v_to_F(nu0,e)
        ToF = np.sqrt((-a)**3/mu) * ((e*np.sinh(F1) - F1) - (e*np.sinh(F0) - F0))

    return ToF












