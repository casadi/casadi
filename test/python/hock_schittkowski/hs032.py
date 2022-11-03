# This is an automatically generated file converting from the apm format

import casadi as cs

def hs032():
    # The optimal objective is (if given in):
    x_opt = cs.DM([0, 0, 1])
    f_opt = 1
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    x0[0]= 0.1
    lbx[0]=0
    x0[1]= 0.7
    lbx[1]=0
    x0[2]= 0.2
    lbx[2]=0
    
    lbg[0] = 3
    g[0] = 6*x[1] + 4*x[2] - x[0]**3
    lbg[1] =  1
    ubg[1] =  1
    g[1] = x[0] + x[1] + x[2]
    obj = (x[0] + 3*x[1] + x[2])**2 + 4*(x[0] - x[1])**2
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
