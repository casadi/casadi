# This is an automatically generated file converting from the apm format

import casadi as cs

def hs029():
    # The optimal objective is (if given in):
    x_opt = cs.DM([4, 2*cs.sqrt(2), 2])
    f_opt = -16*cs.sqrt(2)
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)
    
    x0[0]= 1
    x0[1]= 1
    x0[2]= 1
    
    ubg[0] = 48
    g[0] = x[0]**2 + 2*x[1]**2 + 4*x[2]**2
    obj = -x[0]*x[1]*x[2]
    
    
    f = obj
    print('hs029: Possible optimal points could be (2, -2sqrt(2), -2), (-2, 2sqrt(2), -2), (-4, -2sqrt(2), 2)')
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
