# This is an automatically generated file converting from the apm format

import casadi as cs
import numpy as np

def hs018():
    # The optimal objective is (if given in):
    f_opt = 5
    x_opt = cs.DM([cs.sqrt(250), cs.sqrt(2.5)])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -np.inf*cs.DM.ones(2, 1)
    ubx = np.inf*cs.DM.ones(2, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -np.inf*cs.DM.ones(2, 1)
    ubg = np.inf*cs.DM.ones(2, 1)
    
    x0[0]= 2
    lbx[0]=2
    ubx[0]=50
    x0[1]= 2
    lbx[1]=0
    ubx[1]=50
    
    lbg[0] = 25
    g[0] = x[0]*x[1]
    lbg[1] = 25
    g[1] = x[0]**2 + x[1]**2
    obj = x[0]**2/100 + x[1]**2
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
