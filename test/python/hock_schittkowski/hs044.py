# This is an automatically generated file converting from the apm format

import casadi as cs

def hs044():
    # The optimal objective is (if given in):
    x_opt = cs.DM([0, 3, 0, 4])
    f_opt = -15
    x = cs.MX.sym('x', 4)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(4, 1)
    lbx = -cs.inf*cs.DM.ones(4, 1)
    ubx = cs.inf*cs.DM.ones(4, 1)
    g = cs.MX.zeros(6, 1)
    lbg = -cs.inf*cs.DM.ones(6, 1)
    ubg = cs.inf*cs.DM.ones(6, 1)
    
    x0[0:4]= 0
    lbx[0:4]=0
    
    ubg[0] = 8
    g[0] = x[0] + 2*x[1]
    ubg[1] = 12
    g[1] = 4*x[0] + x[1]
    ubg[2] = 12
    g[2] = 3*x[0] + 4*x[1]
    ubg[3] = 8
    g[3] = 2*x[2] + x[3]
    ubg[4] = 8
    g[4] = x[2] + 2*x[3]
    ubg[5] = 5
    g[5] = x[2] + x[3]
    obj = x[0] - x[1] - x[2] - x[0]*x[2] + x[0]*x[3] + x[1]*x[2] - x[1]*x[3]
    
    
    f = obj
    print('hs044: Other local minimum is f_opt=-13 at (3, 0, 4, 0).')
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
