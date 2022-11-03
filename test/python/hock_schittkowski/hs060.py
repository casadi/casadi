# This is an automatically generated file converting from the apm format

import casadi as cs

def hs060():
    # The optimal objective is (if given in):
    f_opt = 0.03256820025
    x_opt = cs.DM([1.10486, 1.19667, 1.53526])
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)
    
    x0[0:3]= 2
    lbx[0:3]= -10
    ubx[0:3]= 10
    
    lbg[0] = 4 + 3*cs.sqrt(2)
    ubg[0] = 4 + 3*cs.sqrt(2)
    g[0] = x[0]*(1 + x[1]**2) + x[2]**4
    obj = (x[0] - 1)**2 + (x[0] - x[1])**2 + (x[1] - x[2])**4
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
