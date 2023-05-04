# This is an automatically generated file converting from the apm format

import casadi as cs

def hs062():
    # The optimal objective is (if given in):
    x_opt = cs.DM([0.617813, 0.328202, 0.0539851])
    f_opt = -26272.51448
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)
    
    x0[0]= 0.7
    lbx[0]=0
    ubx[0]=1
    x0[1]= 0.2
    lbx[1]=0
    ubx[1]=1
    x0[2]= 0.1
    lbx[2]=0
    ubx[2]=1
    
    lbg[0] =  1
    ubg[0] =  1
    g[0] = x[0] + x[1] + x[2]
    obj = -32.174*(255*cs.log((x[0]+x[1]+x[2]+0.03)/(0.09*x[0] + x[1] + x[2] + 0.03)) +280*cs.log((x[1]+x[2]+0.03)/(0.07*x[1] + x[2] + 0.03))+290*cs.log((x[2]+0.03)/(0.13*x[2] + 0.03)))
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
