# This is an automatically generated file converting from the apm format

import casadi as cs

def hs061():
    # The optimal objective is (if given in):
    f_opt = - 143.6461422
    x_opt = cs.DM([5.32677, -2.119, 3.21046])
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    lbg[0] =  7
    ubg[0] =  7
    g[0] = 3*x[0] - 2*x[1]**2
    lbg[1] = 11
    ubg[1] = 11
    g[1] = 4*x[0] - x[2]**2
    obj = 4*x[0]**2 + 2*x[1]**2 + 2*x[2]**2 - 33*x[0] + 16*x[1] - 24*x[2]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
