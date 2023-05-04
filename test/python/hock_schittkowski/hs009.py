# This is an automatically generated file converting from the apm format

import casadi as cs

def hs009():
    # The optimal objective is (if given in):
    f_opt = -0.5
    x_opt = cs.DM([-471, -628])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)
    
    x0[0]= 0
    x0[1]= 0
    
    lbg[0] =  0
    ubg[0] =  0
    g[0] = 4*x[0] - 3*x[1]
    obj = cs.sin(cs.pi * x[0]/12) * cs.cos(cs.pi * x[1]/16)
    
    
    f = obj
    print('hs009: There are infinitely many solutions (12k-3, 16k-4), k=0 +-1,...')
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
