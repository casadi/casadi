# This is an automatically generated file converting from the apm format

import casadi as cs
import numpy as np

def hs007():
    # The optimal objective is (if given in):
    f_opt = -cs.sqrt(3)
    x_opt = cs.DM([0, cs.sqrt(3)])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -np.inf*cs.DM.ones(2, 1)
    ubx = np.inf*cs.DM.ones(2, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -np.inf*cs.DM.ones(1, 1)
    ubg = np.inf*cs.DM.ones(1, 1)
    
    x0[0]= 2
    x0[1]= 2
    
    lbg[0] =  4
    ubg[0] =  4
    g[0] = (1+x[0]**2)**2 + x[1]**2
    obj = cs.log(1+x[0]**2) - x[1]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
