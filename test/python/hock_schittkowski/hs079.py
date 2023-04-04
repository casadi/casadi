# This is an automatically generated file converting from the apm format

import casadi as cs

def hs079():
    # The optimal objective is (if given in):
    f_opt = 0.0787768209
    x_opt = cs.DM([1.19113, 1.3626, 1.47282, 1.63502, 1.67908])
    x = cs.MX.sym('x', 5)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(5, 1)
    lbx = -cs.inf*cs.DM.ones(5, 1)
    ubx = cs.inf*cs.DM.ones(5, 1)
    g = cs.MX.zeros(3, 1)
    lbg = -cs.inf*cs.DM.ones(3, 1)
    ubg = cs.inf*cs.DM.ones(3, 1)
    
    x0[0:5]= 2
    
    lbg[0] = 2 + 3*cs.sqrt(2)#0
    ubg[0] = 2 + 3*cs.sqrt(2)#0
    g[0] = x[0] + x[1]**2 + x[2]**3  #- ( 2 + 3*cs.sqrt(2))
    lbg[1] = -2 + 2*cs.sqrt(2)#0
    ubg[1] = -2 + 2*cs.sqrt(2)#0
    g[1] = x[1] - x[2]**2 + x[3]   #- ( -2 + 2*cs.sqrt(2))
    lbg[2] =  2
    ubg[2] =  2
    g[2] = x[0]*x[4]
    obj = (x[0]-1)**2 + (x[0]-x[1])**2 + (x[1]-x[2])**2 + (x[2]-x[3])**4 + (x[3]-x[4])**4
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
