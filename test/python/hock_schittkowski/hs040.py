# This is an automatically generated file converting from the apm format

import casadi as cs

def hs040():
    # The optimal objective is (if given in):
    x_opt = cs.DM([2**(-1/3), 2**(1/2), 0.529732, 0.840896])
    f_opt = -0.25
    x = cs.MX.sym('x', 4)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(4, 1)
    lbx = -cs.inf*cs.DM.ones(4, 1)
    ubx = cs.inf*cs.DM.ones(4, 1)
    g = cs.MX.zeros(3, 1)
    lbg = -cs.inf*cs.DM.ones(3, 1)
    ubg = cs.inf*cs.DM.ones(3, 1)
    
    x0[0:4]= 0.8
    
    lbg[0] =  1
    ubg[0] =  1
    g[0] = x[0]**3 + x[1]**2
    lbg[1] =  0
    ubg[1] =  0
    g[1] = x[0]**2*x[3] - x[2]
    lbg[2] =  0
    ubg[2] =  0
    g[2] = x[3]**2 - x[1]
    obj = -x[0]*x[1]*x[2]*x[3]
    
    
    f = obj
    print('hs040: first two components of optimal point are from H&S, the other two from numerical sim.')
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
