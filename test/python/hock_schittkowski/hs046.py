# This is an automatically generated file converting from the apm format

import casadi as cs

def hs046():
    # The optimal objective is (if given in):
    x_opt = cs.DM([1,1,1,1,1])
    f_opt = 0
    x = cs.MX.sym('x', 5)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(5, 1)
    lbx = -cs.inf*cs.DM.ones(5, 1)
    ubx = cs.inf*cs.DM.ones(5, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    x0[0]= cs.sqrt(2)/2
    x0[1]= 1.75
    x0[2]= 0.5
    x0[3]= 2
    x0[4]= 2
    
    lbg[0] =  1
    ubg[0] =  1
    g[0] = x[0]**2*x[3] + cs.sin(x[3] - x[4])
    lbg[1] =  2
    ubg[1] =  2
    g[1] = x[1] + x[2]**4*x[3]**2
    obj = (x[0]-x[1])**2 + (x[2]-1)**2 + (x[3]-1)**4 + (x[4]-1)**6
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
