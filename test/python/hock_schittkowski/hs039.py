# This is an automatically generated file converting from the apm format

import casadi as cs

def hs039():
    # The optimal objective is (if given in): 
    x_opt = cs.DM([1, 1, 0, 0])
    f_opt = -1
    x = cs.MX.sym('x', 4)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(4, 1)
    lbx = -cs.inf*cs.DM.ones(4, 1)
    ubx = cs.inf*cs.DM.ones(4, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    x0[0:4]= 2
    
    lbg[0] =  0
    ubg[0] =  0
    g[0] = x[1] - x[0]**3 - x[2]**2
    lbg[1] =  0
    ubg[1] =  0
    g[1] = x[0]**2 - x[1] - x[3]**2
    obj = -x[0]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
