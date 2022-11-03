# This is an automatically generated file converting from the apm format

import casadi as cs

def hs037():
    # The optimal objective is (if given in):
    x_opt = cs.DM([24, 12, 12])
    f_opt = -3456
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    x0[0:3]= 10
    lbx[0:3]=0
    ubx[0:3]=42
    
    ubg[0] = 72
    g[0] = x[0] + 2*x[1] + 2*x[2]
    lbg[1] = 0
    g[1] = x[0] + 2*x[1] + 2*x[2]
    obj = -x[0]*x[1]*x[2]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
