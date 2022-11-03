# This is an automatically generated file converting from the apm format

import casadi as cs

def hs058():
    # The optimal objective is (if given in):
    f_opt = 3.19033354957
    x_opt = cs.DM([-0.786151, 0.618034])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    g = cs.MX.zeros(3, 1)
    lbg = -cs.inf*cs.DM.ones(3, 1)
    ubg = cs.inf*cs.DM.ones(3, 1)
    
    x0[0]= -2
    lbx[0]=-2
    ubx[0]=0.5
    x0[1]= 1
    
    lbg[0] = 0
    g[0] =  x[1]**2 - x[0]
    lbg[1] = 0
    g[1] =  x[0]**2 - x[1]
    lbg[2] = 1
    g[2] =  x[0]**2 + x[1]**2
    obj = 100*(x[1]-x[0]**2)**2+(1-x[0])**2
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
