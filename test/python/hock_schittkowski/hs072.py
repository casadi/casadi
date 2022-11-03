# This is an automatically generated file converting from the apm format

import casadi as cs

def hs072():
    # The optimal objective is (if given in):
    f_opt = 727.67937
    x_opt = cs.DM([193.4071, 179.5475, 185.0186, 168.7062])
    x = cs.MX.sym('x', 4)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(4, 1)
    lbx = -cs.inf*cs.DM.ones(4, 1)
    ubx = cs.inf*cs.DM.ones(4, 1)
    a = cs.DM.zeros(2, 4)
    b = cs.DM.zeros(2, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    a[0, 0] = 4
    a[0, 1] = 2.25
    a[0, 2] = 1
    a[0, 3] = 0.25
    a[1, 0] = 0.16
    a[1, 1] = 0.36
    a[1, 2] = 0.64
    a[1, 3] = 0.64
    b[0] = 0.0401
    b[1] = 0.010085
    
    x0[0:4]= 1
    lbx[0:4] = 0.001
    ubx[0] = 4e5
    ubx[1] = 3e5
    ubx[2] = 2e5
    ubx[3] = 1e5
    
    ubg[0] = b[0]
    g[0] = a[0, 0]/x[0] + a[0, 1]/x[1] + a[0, 2]/x[2] + a[0, 3]/x[3]
    ubg[1] = b[1]
    g[1] = a[1, 0]/x[0] + a[1, 1]/x[1] + a[1, 2]/x[2] + a[1, 3]/x[3]
    obj = 1 + x[0] + x[1] + x[2] + x[3]
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
