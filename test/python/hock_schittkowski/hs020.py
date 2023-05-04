# This is an automatically generated file converting from the apm format

import casadi as cs

def hs020():
    # The optimal objective is (if given in):
    #f_opt = 81.5 - 25*cs.sqrt(3)
    f_opt_local = 40.198727306534536
    x_opt_local = cs.DM([-0.5, 0.5*cs.sqrt(3)])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    g = cs.MX.zeros(3, 1)
    lbg = -cs.inf*cs.DM.ones(3, 1)
    ubg = cs.inf*cs.DM.ones(3, 1)
    
    x0[0]= -2
    x0[1]= 1
    
    lbg[0] = 0
    g[0] = x[0] + x[1]**2
    lbg[1] = 0
    g[1] = x[0]**2 + x[1]
    lbg[2] = 1
    g[2] = x[0]**2 + x[1]**2
    lbx[0] = -1/2
    ubx[0] =  1/2
    obj = 100*(x[1] - x[0]**2)**2 + (1-x[0])**2
    
    
    f = obj
    
    print('hs020: Problem starts closer to local minimum.')
    print('hs020: H&S give global solution (0.5, 0.5*sqrt(3)) with opt f = 81.5 -25*sqrt(3)')
    return (x_opt_local, f_opt_local, x, f, g, lbg, ubg, lbx, ubx, x0)
