# This is an automatically generated file converting from the apm format

import casadi as cs

def hs043():
    # The optimal objective is (if given in):
    x_opt = cs.DM([0, 1, 2, -1])
    f_opt = -44
    x = cs.MX.sym('x', 4)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(4, 1)
    lbx = -cs.inf*cs.DM.ones(4, 1)
    ubx = cs.inf*cs.DM.ones(4, 1)
    g = cs.MX.zeros(3, 1)
    lbg = -cs.inf*cs.DM.ones(3, 1)
    ubg = cs.inf*cs.DM.ones(3, 1)
    
    x0[0:4]= 0
    
    ubg[0] = 8
    g[0] = x[0]**2 +x[1]**2 +x[2]**2 +x[3]**2 +x[0] -x[1] +x[2] -x[3]
    ubg[1] = 10
    g[1] = x[0]**2 +2*x[1]**2 +x[2]**2 +2*x[3]**2 -x[0] -x[3]
    ubg[2] = 5
    g[2] = 2*x[0]**2 +x[1]**2 +x[2]**2 +2*x[0] -x[1] -x[3]
    obj = x[0]**2 + x[1]**2 + 2*x[2]**2 + x[3]**2 - 5*x[0] - 5*x[1] - 21*x[2] + 7*x[3]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
