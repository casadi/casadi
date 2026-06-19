# This is an automatically generated file converting from the apm format

import casadi as cs
import numpy as np

def hs022():
    # The optimal objective is (if given in):
    f_opt = 1
    x_opt = cs.DM([1,1])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -np.inf*cs.DM.ones(2, 1)
    ubx = np.inf*cs.DM.ones(2, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -np.inf*cs.DM.ones(2, 1)
    ubg = np.inf*cs.DM.ones(2, 1)
    
    x0[0]= 2
    x0[1]= 2
    
    ubg[0] = 2
    g[0] = x[0] + x[1]
    lbg[1] = 0
    g[1] = -x[0]**2 + x[1]
    obj = (x[0] - 2)**2 + (x[1] - 1)**2
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
