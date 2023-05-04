# This is an automatically generated file converting from the apm format

import casadi as cs

def hs011():
    # The optimal objective is (if given in):
    f_opt = -8.498464223
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)
    
    x0[0]= 4.9
    x0[1]= 0.1
    
    ubg[0] = 0
    g[0] = x[0]**2  - x[1]
    obj = (x[0] - 5)**2 + x[1]**2 -25

    #a = 7.5*cs.sqrt(6) + cs.sqrt(338.5)
    #x_opt = cs.DM([(a - 1/a)/cs.sqrt(6), (a**2 - 2 + a**(-2))/6])
    x_opt = [1.23477, 1.52466]
    
    f = obj
    print('hs011: Acc. to H&S optimal solution is [(a - 1/a)/cs.sqrt(6), (a**2 - 2 + a**(-2))/6]')
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
