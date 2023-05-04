# This is an automatically generated file converting from the apm format

import casadi as cs

def hs005():
    # The optimal objective is (if given in):
    x_opt = cs.DM([-cs.pi/3 + 0.5, -cs.pi/3 - 0.5])
    f_opt = -(cs.sqrt(3)/2 + cs.pi/3)
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    g = []
    lbg = []
    ubg = []
    
    x0[0]= 0
    lbx[0]= -1.5
    ubx[0]= 4
    x0[1]= 0
    lbx[1]= -3
    ubx[1]=3
    
    obj = cs.sin(x[0]+x[1]) + (x[0]-x[1])**2 - 1.5*x[0] + 2.5*x[1] + 1
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
