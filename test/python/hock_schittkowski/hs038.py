# This is an automatically generated file converting from the apm format

import casadi as cs

def hs038():
    # The optimal objective is (if given in):
    x_opt = cs.DM([1, 1, 1, 1])
    f_opt = 0
    x = cs.MX.sym('x', 4)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(4, 1)
    lbx = -cs.inf*cs.DM.ones(4, 1)
    ubx = cs.inf*cs.DM.ones(4, 1)
    g = []
    lbg = []
    ubg = []
    
    x0[0]= -3
    lbx[0]=-10
    ubx[0]=10
    x0[1]= -1
    lbx[1]=-10
    ubx[1]=10
    x0[2]= -3
    lbx[2]=-10
    ubx[2]=10
    x0[3]= -1
    lbx[3]=-10
    ubx[3]=10
    
    obj = 100*(x[1]-x[0]**2)**2 + (1-x[0])**2 + 90*(x[3]-x[2]**2)**2 + (1-x[2])**2 + 10.1*( (x[1]-1)**2 + (x[3]-1)**2 ) + 19.8*(x[1]-1)*(x[3]-1)
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
