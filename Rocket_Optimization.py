def R_opt(I_sps, eps, V_star):

    # Import libraries
    import numpy as np
    import sympy as sp
    from sympy import nsolve, Symbol, Eq

    import matplotlib.pyplot as plt

    ###################################
    #              Setup              #
    ###################################
    g = 9.81/1000
    I_sp = np.array(I_sps)
    ep = np.array(eps)


    ###################################
    #          Calculations           #
    ###################################

    V_e = np.zeros(3)
    #import pdb; pdb.set_trace()

    for i in range(0,len(I_sp)):
        V_e[i] = I_sp[i]*g

    lam = Symbol('lam')

    Equation = 0 #Symbol('Equation')
    for i in range(0, len(I_sp)):
        Equation = Equation + V_e[i]*sp.log(ep[i] - ep[i]/(1 + lam*V_e[i]))

    print(Equation)

    func = Eq(0,Equation + V_star)
    sol = nsolve(func,lam)

    #print(sol)





Isp=[290, 350, 455]
eps = [0.10, 0.12, 0.18]
V_star = 12

R_opt(Isp, eps, V_star)













