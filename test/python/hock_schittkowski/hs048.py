# This is an automatically generated file converting from the apm format

import casadi as cs

def hs048():
    # The optimal objective is (if given in):
    x_opt = cs.DM([1,1,1,1,1])
    f_opt = 0
    x = cs.MX.sym('x', 5)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(5, 1)
    lbx = -cs.inf*cs.DM.ones(5, 1)
    ubx = cs.inf*cs.DM.ones(5, 1)
    p = cs.MX.zeros(5, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    x0[0]= 3
    x0[1]= 5
    x0[2]= -3
    x0[3]= 2
    x0[4]= -2
    
    lbg[0] =  5
    ubg[0] =  5
    g[0] = cs.sum1(x)
    lbg[1] =  -3
    ubg[1] =  -3
    g[1] = x[2] - 2*(x[3]+x[4])
    obj = (x[0]-1)**2 + (x[1]-x[2])**2 + (x[3]-x[4])**2
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
