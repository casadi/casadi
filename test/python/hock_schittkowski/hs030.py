# This is an automatically generated file converting from the apm format

import casadi as cs

def hs030():
    # The optimal objective is (if given in):
    x_opt = cs.DM([1, 0, 0])
    f_opt = 1
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)
    
    x0[0]= 1
    x0[1]= 1
    x0[2]= 1
    
    lbg[0] = 1
    g[0] = x[0]**2 + x[1]**2
    lbx[0] = 1
    ubx[0] =  10
    lbx[1] = -10
    ubx[1] =  10
    lbx[2] = -10
    ubx[2] =  10
    obj = x[0]**2 + x[1]**2 + x[2]**2
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
