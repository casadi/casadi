# This is an automatically generated file converting from the apm format

import casadi as cs

def hs035():
    # The optimal objective is (if given in):
    x_opt = cs.DM([4/3, 7/9, 4/9])
    f_opt = 1/9
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)
    
    x0[0]= 0.5
    lbx[0]=0
    x0[1]= 0.5
    lbx[1]=0
    x0[2]= 0.5
    lbx[2]=0
    
    ubg[0] = 3
    g[0] = x[0] + x[1] + 2*x[2]
    obj = 9 - 8*x[0] - 6*x[1] - 4*x[2] + 2*x[0]**2 + 2*x[1]**2 + x[2]**2 + 2*x[0]*x[1] + 2*x[0]*x[2]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
