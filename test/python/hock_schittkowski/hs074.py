# This is an automatically generated file converting from the apm format

import casadi as cs

def hs074():
    # The optimal objective is (if given in):
    f_opt = 5126.4981
    x_opt = cs.DM([679.945, 1026.07, 0.118876, -0.396234])
    x = cs.MX.sym('x', 4)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(4, 1)
    lbx = -cs.inf*cs.DM.ones(4, 1)
    ubx = cs.inf*cs.DM.ones(4, 1)
    l = cs.DM.zeros(4, 1)
    u = cs.DM.zeros(4, 1)
    g = cs.MX.zeros(4, 1)
    lbg = -cs.inf*cs.DM.ones(4, 1)
    ubg = cs.inf*cs.DM.ones(4, 1)
    
    a = 0.55
    l[0] = 0
    l[1] = 0
    l[2] = -a
    l[3] = -a
    u[0] = 1200
    u[1] = 1200
    u[2] = a
    u[3] = a
    
    x0[0:4]= 0
    lbx[0:4]=l[0:4]
    ubx[0:4]=u[0:4]
    
    g[0] =  x[3] - x[2]
    lbg[0] = -a
    ubg[0] =  a
    lbg[1] = 894.8
    ubg[1] = 894.8
    g[1] = x[0]  - 1000*cs.sin(-x[2] - 0.25) - 1000*cs.sin(-x[3] - 0.25)
    lbg[2] = 894.8
    ubg[2] = 894.8
    g[2] = x[1]  - 1000*cs.sin(x[2] - 0.25) - 1000*cs.sin(x[2]-x[3] - 0.25)
    lbg[3] =  -1294.8
    ubg[3] =  -1294.8
    g[3] = 1000*cs.sin(x[3] - 0.25) + 1000*cs.sin(x[3] - x[2] - 0.25)
    obj = 3*x[0] + 1.0e-6*x[0]**3 + 2*x[1] + 2.0e-6*x[1]**3/3
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
