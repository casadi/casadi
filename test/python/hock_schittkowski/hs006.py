# This is an automatically generated file converting from the apm format

import casadi as cs

def hs006():
    # The optimal objective is (if given in):
    f_opt = 0
    x_opt = cs.DM([1, 1])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)
    
    x0[0]= -1.2
    x0[1]= 1
    
    lbg[0] =  0
    ubg[0] =  0
    g[0] = 10*(x[1] - x[0]**2)
    obj = (1-x[0])**2
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
