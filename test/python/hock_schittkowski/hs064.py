# This is an automatically generated file converting from the apm format

import casadi as cs

def hs064():
    # The optimal objective is (if given in):
    f_opt = 6299.842428
    x_opt = cs.DM([108.7347175, 85.12613942, 204.3247078])
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)
    
    x0[0:3]= 1
    lbx[0:3]= 1e-5
    
    ubg[0] = 1
    g[0] = 4/x[0] + 32/x[1] + 120/x[2]
    obj = 5*x[0] + 50000/x[0] + 20*x[1] + 72000/x[1] + 10*x[2] + 144000/x[2]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
