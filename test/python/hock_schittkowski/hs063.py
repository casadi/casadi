# This is an automatically generated file converting from the apm format

import casadi as cs

def hs063():
    # The optimal objective is (if given in):
    f_opt = 961.7151721
    x_opt = cs.DM([3.51212, 0.216988, 3.55217])
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    x0[0:3]= 2
    lbx[0:3]=0
    
    lbg[0] =  56
    ubg[0] =  56
    g[0] = 8*x[0] + 14*x[1] + 7*x[2]
    lbg[1] =  25
    ubg[1] =  25
    g[1] = x[0]**2 +  x[1]**2 + x[2]**2
    obj = 1000 - x[0]**2 - 2*x[1]**2 - x[2]**2 - x[0]*x[1] - x[0]*x[2]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
