# This is an automatically generated file converting from the apm format

import casadi as cs

def hs073():
    # The optimal objective is (if given in):
    f_opt = 29.894378
    x_opt = cs.DM([0.635522, 2.06628e-10, 0.312702, 0.0517765])
    x = cs.MX.sym('x', 4)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(4, 1)
    lbx = -cs.inf*cs.DM.ones(4, 1)
    ubx = cs.inf*cs.DM.ones(4, 1)
    g = cs.MX.zeros(3, 1)
    lbg = -cs.inf*cs.DM.ones(3, 1)
    ubg = cs.inf*cs.DM.ones(3, 1)
    
    x0[0:4]= 1
    lbx[0:4]=0
    
    lbg[0] = 5
    g[0] = 2.3*x[0] + 5.6*x[1] + 11.1*x[2] + 1.3*x[3]
    lbg[1] = 21
    g[1] = 12*x[0] + 11.9*x[1] + 41.8*x[2] + 52.1*x[3]  - 1.645*cs.sqrt(0.28*x[0]**2 + 0.19*x[1]**2 + 20.5*x[2]**2 + 0.62*x[3]**2)
    lbg[2] =  1
    ubg[2] =  1
    g[2] = x[0]+x[1]+x[2]+x[3]
    obj = 24.55*x[0] + 26.75*x[1] + 39*x[2] + 40.50*x[3]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
