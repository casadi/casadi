# This is an automatically generated file converting from the apm format

import casadi as cs

def hs042():
    # The optimal objective is (if given in):
    x_opt = cs.DM([2, 2, 0.6*cs.sqrt(2), 0.8*cs.sqrt(2)])
    f_opt = 28 - 10*cs.sqrt(2)
    x = cs.MX.sym('x', 4)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(4, 1)
    lbx = -cs.inf*cs.DM.ones(4, 1)
    ubx = cs.inf*cs.DM.ones(4, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)
    
    x0[0:4]= 1
    lbx[0:4]=0
    
    lbx[0] =  2
    ubx[0] =  2
    lbg[0] =  2
    ubg[0] =  2
    g[0] = x[2]**2 + x[3]**2
    obj = (x[0]-1)**2 + (x[1]-2)**2 + (x[2]-3)**2 + (x[3]-4)**2
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
