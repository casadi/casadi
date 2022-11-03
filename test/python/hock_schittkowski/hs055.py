# This is an automatically generated file converting from the apm format

import casadi as cs

def hs055():
    # The optimal objective is (if given in):
    x_opt = cs.DM([0, 4/3, 5/3, 1, 2/3, 1/3])
    f_opt = 19/3
    x = cs.MX.sym('x', 6)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(6, 1)
    lbx = -cs.inf*cs.DM.ones(6, 1)
    ubx = cs.inf*cs.DM.ones(6, 1)
    u = cs.DM.zeros(6, 1)
    g = cs.MX.zeros(6, 1)
    lbg = -cs.inf*cs.DM.ones(6, 1)
    ubg = cs.inf*cs.DM.ones(6, 1)
    
    x0[0]= 1
    lbx[0]=0
    ubx[0]=1
    x0[1]= 2
    lbx[1]=0
    x0[2:5]= 0
    lbx[2:5]=0
    lbx[2] = 1 # revised not original bound
    ubx[3] =  1
    x0[5]= 2

    lbg[0] =  6
    ubg[0] =  6
    g[0] = x[0] + 2*x[1] + 5*x[4]
    lbg[1] =  3
    ubg[1] =  3
    g[1] = x[0] + x[1] + x[2]
    lbg[2] =  2
    ubg[2] =  2
    g[2] = x[3] + x[4] + x[5]
    lbg[3] =  1
    ubg[3] =  1
    g[3] = x[0] + x[3]
    lbg[4] =  2
    ubg[4] =  2
    g[4] = x[1] + x[4]
    lbg[5] =  2
    ubg[5] =  2
    g[5] = x[2] + x[5]
    f = x[0] + 2*x[1] + 4*x[4] + cs.exp(x[0]*x[3])
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
