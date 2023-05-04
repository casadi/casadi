# This is an automatically generated file converting from the apm format

import casadi as cs

def hs002():
    # The optimal objective is (if given in):
    f_opt = 0.0504261879 #That is the global solution, but for this starting point it is a different local solution
    f_opt_local = 4.9412294
    x_opt_local = cs.DM([-1.22103, 1.5])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    g = []
    lbg = []
    ubg = []
    
    x0[0]= -2
    x0[1]= 1
    lbx[1]= 1.5
    
    obj = 100*(x[1] - x[0]**2)**2 + (1-x[0])**2
    
    
    f = obj
    print('hs002: Problem should converge to local solution, but H&S give global solution: 0.0504261879')
    print('hs002: global x_opt: [1.2243707, 1.5]')
    return (x_opt_local, f_opt_local, x, f, g, lbg, ubg, lbx, ubx, x0)
