# This is an automatically generated file converting from the apm format

import casadi as cs

def hs023():
    # The optimal objective is (if given in):
    f_opt = 2
    x_opt = cs.DM([1,1])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    g = cs.MX.zeros(5, 1)
    lbg = -cs.inf*cs.DM.ones(5, 1)
    ubg = cs.inf*cs.DM.ones(5, 1)
    
    x0[0]= 3
    lbx[0]=-50
    ubx[0]=50
    x0[1]= 1
    lbx[1]=-50
    ubx[1]=50
    
    lbg[0] = 1
    g[0] = x[0] + x[1]
    lbg[1] = 1
    g[1] = x[0]**2 + x[1]**2
    lbg[2] = 9
    g[2] = 9*x[0]**2 + x[1]**2
    lbg[3] = 0
    g[3] = x[0]**2 - x[1]
    lbg[4] = 0
    g[4] = x[1]**2 - x[0]
    obj = x[0]**2 + x[1]**2
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
