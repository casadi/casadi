# This is an automatically generated file converting from the apm format

import casadi as cs

def hs004():
    # The optimal objective is (if given in):
    f_opt = 8/3
    x_opt = cs.DM([1, 0])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    g = []
    lbg = []
    ubg = []
    
    x0[0]= 1.125
    lbx[0]= 1
    x0[1]= 0.125
    lbx[1]= 0
    
    obj = (x[0]+1)**3/3 + x[1]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
