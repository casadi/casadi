# This is an automatically generated file converting from the apm format

import casadi as cs

def hs019():
    # The optimal objective is (if given in):
    f_opt = -6961.81381
    x_opt = cs.DM([14.095, 0.84296079])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    x0[0]= 20.1
    x0[1]= 5.84
    
    lbg[0] = 100
    g[0] = (x[0] - 5)**2 + (x[1] - 5)**2
    ubg[1] = 82.81
    g[1] = (x[1] - 5)**2 + (x[0] - 6)**2
    lbx[0] = 13
    ubx[0] =  100
    lbx[1] = 0
    ubx[1] =  100
    obj = (x[0] - 10)**3 + (x[1] - 20)**3
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
