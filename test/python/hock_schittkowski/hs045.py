# This is an automatically generated file converting from the apm format

import casadi as cs

def hs045():
    # The optimal objective is (if given in):
    x_opt = cs.DM([1,2,3,4,5])
    f_opt = 1
    x = cs.MX.sym('x', 5)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(5, 1)
    lbx = -cs.inf*cs.DM.ones(5, 1)
    ubx = cs.inf*cs.DM.ones(5, 1)
    g = []
    lbg = []
    ubg = []
    
    x0[0]= 0
    lbx[0]=0
    ubx[0]=1
    x0[1]= 0
    lbx[1]=0
    ubx[1]=2
    x0[2]= 0
    lbx[2]=0
    ubx[2]=3
    x0[3]= 0
    lbx[3]=0
    ubx[3]=4
    x0[4]= 0
    lbx[4]=0
    ubx[4]=5
    
    obj = 2 - x[0]*x[1]*x[2]*x[3]*x[4]/120
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
