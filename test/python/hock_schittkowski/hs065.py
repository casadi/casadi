# This is an automatically generated file converting from the apm format

import casadi as cs

def hs065():
    # The optimal objective is (if given in):
    f_opt = 0.9535288567
    x_opt = cs.DM([3.65046, 3.65046, 4.62042])
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)
    
    x0[0]= -5
    x0[1]=  5
    x0[2]=  0
    
    ubg[0] = 48
    g[0] = x[0]**2 + x[1]**2 + x[2]**2
    lbx[0] = -4.5
    ubx[0] =  4.5
    lbx[1] = -4.5
    ubx[1] =  4.5
    lbx[2] = -5
    ubx[2] =    5
    obj = (x[0] - x[1])**2 + (x[0] + x[1] - 10)**2/9 + (x[2] - 5)**2
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
