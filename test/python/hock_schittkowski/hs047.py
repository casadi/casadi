# This is an automatically generated file converting from the apm format

import casadi as cs

def hs047():
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
    
    x0[0]= 2
    x0[1]= cs.sqrt(2)
    x0[2]= -1
    x0[3]= 2-cs.sqrt(2)
    x0[4]= 1/2
    
    lbg[0] =  3
    ubg[0] =  3
    g[0] = x[0] + x[1]**2 + x[2]**3
    lbg[1] =  1
    ubg[1] =  1
    g[1] = x[1] - x[2]**2 + x[3]
    lbg[2] =  1
    ubg[2] =  1
    g[2] = x[0]*x[4]
    obj = (x[0]-x[1])**2 + (x[1]-x[2])**3 + (x[2]-x[3])**4 + (x[3]-x[4])**4
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
