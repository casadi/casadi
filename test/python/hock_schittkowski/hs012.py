# This is an automatically generated file converting from the apm format

import casadi as cs

def hs012():
    # The optimal objective is (if given in):
    f_opt = -30
    x_opt = cs.DM([2, 3])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)
    
    x0[0]= 0
    x0[1]= 0
    
    ubg[0] = 25
    g[0] = 4*x[0]**2 + x[1]**2
    obj = (x[0]**2)/2 + x[1]**2 - x[0]*x[1] - 7*x[0] - 7*x[1]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
