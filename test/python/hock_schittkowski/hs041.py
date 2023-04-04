# This is an automatically generated file converting from the apm format

import casadi as cs

def hs041():
    # The optimal objective is (if given in):
    x_opt = cs.DM([2/3, 1/3, 1/3, 2])
    f_opt = 52/27
    x = cs.MX.sym('x', 4)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(4, 1)
    lbx = -cs.inf*cs.DM.ones(4, 1)
    ubx = cs.inf*cs.DM.ones(4, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)
    
    x0[0:4]= 2
    lbx[0:4]= 0
    
    lbg[0] =  0
    ubg[0] =  0
    g[0] = x[0] + 2*x[1] + 2*x[2] - x[3]
    ubx[0] =  1
    ubx[1] =  1
    ubx[2] =  1
    ubx[3] =  2
    obj = 2-x[0]*x[1]*x[2]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
