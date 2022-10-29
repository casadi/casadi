# This is an automatically generated file converting from the apm format

import casadi as cs

def hs024():
    # The optimal objective is (if given in):
    f_opt = -1
    x_opt = cs.DM([3, cs.sqrt(3)])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    g = cs.MX.zeros(3, 1)
    lbg = -cs.inf*cs.DM.ones(3, 1)
    ubg = cs.inf*cs.DM.ones(3, 1)
    
    x0[0]= 1
    lbx[0]=0
    x0[1]= 1/2
    lbx[1]=0
    
    lbg[0] = 0
    g[0] = x[0]/cs.sqrt(3) - x[1]
    lbg[1] = 0
    g[1] = x[0] + cs.sqrt(3)*x[1]
    lbg[2] = -6
    g[2] = -x[0] - cs.sqrt(3)*x[1]
    obj = ((x[0] - 3)**2 - 9) * x[1]**3 / (27*cs.sqrt(3))
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
