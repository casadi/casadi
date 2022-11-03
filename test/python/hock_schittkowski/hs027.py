# This is an automatically generated file converting from the apm format

import casadi as cs

def hs027():
    # The optimal objective is (if given in):
    x_opt = cs.DM([-1,1,0])
    f_opt = 0.04
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)
    
    x0[0]= 2
    x0[1]= 2
    x0[2]= 2
    
    lbg[0] =  -1
    ubg[0] =  -1
    g[0] = x[0] + x[2]**2
    obj = (x[0] - 1)**2/100 + (x[1] - x[0]**2)**2
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
