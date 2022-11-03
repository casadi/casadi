# This is an automatically generated file converting from the apm format

import casadi as cs

def hs033():
    # The optimal objective is (if given in):
    x_opt = cs.DM([0, cs.sqrt(2), cs.sqrt(2)])
    f_opt = cs.sqrt(2) - 6
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    x0[0]= 0
    lbx[0]=0
    x0[1]= 0
    lbx[1]=0
    x0[2]= 3
    lbx[2]=0
    
    ubg[0] = 0
    g[0] = x[0]**2 + x[1]**2  - x[2]**2
    lbg[1] = 4
    g[1] = x[0]**2 + x[1]**2 + x[2]**2
    ubx[2] =  5
    obj = (x[0] - 1)*(x[0] - 2)*(x[0] - 3) + x[2]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
