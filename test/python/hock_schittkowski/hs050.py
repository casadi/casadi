# This is an automatically generated file converting from the apm format

import casadi as cs

def hs050():
    # The optimal objective is (if given in):
    x_opt = cs.DM([1,1,1,1,1])
    f_opt = 0
    x = cs.MX.sym('x', 5)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(5, 1)
    lbx = -cs.inf*cs.DM.ones(5, 1)
    ubx = cs.inf*cs.DM.ones(5, 1)
    g = cs.MX.zeros(3, 1)
    lbg = -cs.inf*cs.DM.ones(3, 1)
    ubg = cs.inf*cs.DM.ones(3, 1)
    
    x0[0]=  35
    x0[1]= -31
    x0[2]=  11
    x0[3]=   5
    x0[4]=  -5
    
    lbg[0] =  6
    ubg[0] =  6
    g[0] = x[0] + 2*x[1] + 3*x[2]
    lbg[1] =  6
    ubg[1] =  6
    g[1] = x[1] + 2*x[2] + 3*x[3]
    lbg[2] =  6
    ubg[2] =  6
    g[2] = x[2] + 2*x[3] + 3*x[4]
    obj = (x[0]-x[1])**2 + (x[1]-x[2])**2 + (x[2]-x[3])**4 + (x[3]-x[4])**2
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
