# This is an automatically generated file converting from the apm format

import casadi as cs

def hs014():
    # The optimal objective is (if given in):
    f_opt = 9 - 2.875*cs.sqrt(7)
    x_opt = cs.DM([0.5*(cs.sqrt(7)-1), 0.25*(cs.sqrt(7)+1)])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    x0[0]= 2
    lbx[0]=0
    x0[1]= 2
    lbx[1]=0
    
    ubg[0] = 1
    g[0] = x[0]**2/4 + x[1]**2
    lbg[1] =  -1
    ubg[1] =  -1
    g[1] = x[0] - 2*x[1]
    obj = (x[0] - 2)**2 + (x[1]-1)**2
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
