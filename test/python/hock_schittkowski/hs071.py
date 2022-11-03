# This is an automatically generated file converting from the apm format

import casadi as cs

def hs071():
    # The optimal objective is (if given in):
    x_opt = cs.DM([1, 4.743, 3.82115, 1.37941])
    f_opt = 17.0140173
    x = cs.MX.sym('x', 4)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(4, 1)
    lbx = -cs.inf*cs.DM.ones(4, 1)
    ubx = cs.inf*cs.DM.ones(4, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    lbx[0]=1
    ubx[0]=5
    x0[1]= 5
    lbx[1]=1
    ubx[1]=5
    x0[2]= 5
    lbx[2]=1
    ubx[2]=5
    x0[3]= 1
    lbx[3]=1
    ubx[3]=5
    
    lbg[0] = 25
    g[0] = x[0]*x[1]*x[2]*x[3]
    lbg[1] =  40
    ubg[1] =  40
    g[1] = x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2
    obj = x[0]*x[3]*(x[0]+x[1]+x[2]) + x[2]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
