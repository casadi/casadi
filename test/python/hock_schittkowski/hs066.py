# This is an automatically generated file converting from the apm format

import casadi as cs

def hs066():
    # The optimal objective is (if given in):
    f_opt = 0.5181632741
    x_opt = cs.DM([0.184127, 1.20217, 3.32732])
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    x0[0]= 0
    x0[1]= 1.05
    x0[2]= 2.9
    
    lbg[0] = 0
    g[0] = x[1] - cs.exp(x[0])
    lbg[1] = 0
    g[1] = x[2] - cs.exp(x[1])
    lbx[0] = 0
    ubx[0] =  100
    lbx[1] = 0
    ubx[1] =  100
    lbx[2] = 0
    ubx[2] =   10
    obj = 0.2*x[2] - 0.8*x[0]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
