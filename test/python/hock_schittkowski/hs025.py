# This is an automatically generated file converting from the apm format

import casadi as cs

def hs025():
    # The optimal objective is (if given in):
    x_opt = cs.DM([50, 25, 1.5])
    f_opt = 0
    x = cs.MX.sym('x', 3)
    obj = cs.MX.zeros(99, 1)
    x0 = cs.DM.zeros(3, 1)
    lbx = -cs.inf*cs.DM.ones(3, 1)
    ubx = cs.inf*cs.DM.ones(3, 1)
    s = cs.MX.zeros(99, 1)
    u = cs.MX.zeros(99, 1)
    g = []
    lbg = []
    ubg = []
    
    x0[0]= 100
    x0[1]= 12.5
    x0[2]= 3
    
    s = cs.DM(range(1,100))
    u[0:99] = 25 + (-50*cs.log(s[0:99]/100))**(2/3)
    
    lbx[0] = 1/10
    ubx[0] =  100
    lbx[1] = 0
    ubx[1] =  25.6
    lbx[2] = 0
    ubx[2] =  5
    obj[0:99] = (-s[0:99]/100 + cs.exp(-(u[0:99] - x[1])**x[2]/x[0]))**2
    
    
    f = cs.sum1(obj)
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
