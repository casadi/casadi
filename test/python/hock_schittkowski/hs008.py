# This is an automatically generated file converting from the apm format

import casadi as cs

def hs008():
    # The optimal objective is (if given in):
    f_opt = -1
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    x0[0]= 2
    x0[1]= 1
    
    lbg[0] =  25
    ubg[0] =  25
    g[0] = x[0]**2 + x[1]**2
    lbg[1] =  9
    ubg[1] =  9
    g[1] = x[0]*x[1]
    obj = -1
    
    a = cs.sqrt((25 + cs.sqrt(301))/2)
    b = cs.sqrt((25 - cs.sqrt(301))/2)
    x_opt = cs.DM([a, 9/a])
    
    f = obj
    print('hs008: There are more feasible points: (-a, -9/a), (b, 9/b), (-b, -9/b)')
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
