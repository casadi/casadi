# This is an automatically generated file converting from the apm format

import casadi as cs

def hs052():
    # The optimal objective is (if given in):
    x_opt = cs.DM([-33/349, 11/349, 180/349, -158/349, 11/349])
    f_opt = 1859/349
    x = cs.MX.sym('x', 5)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(5, 1)
    lbx = -cs.inf*cs.DM.ones(5, 1)
    ubx = cs.inf*cs.DM.ones(5, 1)
    g = cs.MX.zeros(3, 1)
    lbg = -cs.inf*cs.DM.ones(3, 1)
    ubg = cs.inf*cs.DM.ones(3, 1)
    
    x0[0:5]= 2
    
    lbg[0] =  0
    ubg[0] =  0
    g[0] = x[0] + 3*x[1]
    lbg[1] =  0
    ubg[1] =  0
    g[1] = x[2] +   x[3] - 2*x[4]
    lbg[2] =  0
    ubg[2] =  0
    g[2] = x[1] -   x[4]
    obj = (4*x[0]-x[1])**2 + (x[1]+x[2]-2)**2 + (x[3]-1)**2 + (x[4]-1)**2
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
