# This is an automatically generated file converting from the apm format

import casadi as cs

def hs056():
    # The optimal objective is (if given in):
    f_opt = -3.456
    x_opt = cs.DM([2.4, 1.2, 1.2, ])
    x = cs.MX.sym('x', 7)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(7, 1)
    lbx = -cs.inf*cs.DM.ones(7, 1)
    ubx = cs.inf*cs.DM.ones(7, 1)
    xi = cs.DM.zeros(7, 1)
    g = cs.MX.zeros(4, 1)
    lbg = -cs.inf*cs.DM.ones(4, 1)
    ubg = cs.inf*cs.DM.ones(4, 1)
    
    a = cs.asin(cs.sqrt(1/4.2))
    b = cs.asin(cs.sqrt(5/7.2))
    c = cs.asin(cs.sqrt(4/7))
    d = cs.asin(cs.sqrt(2/7))
    xi[0:3] = 1
    xi[3:6] = a
    xi[6] = b
    
    x0[0:7]=  xi[0:7]
    lbx[0:7]= 0
    
    lbg[0] =  0
    ubg[0] =  0
    g[0] = x[0] - 4.2*cs.sin(x[3])**2
    lbg[1] =  0
    ubg[1] =  0
    g[1] = x[1] - 4.2*cs.sin(x[4])**2
    lbg[2] =  0
    ubg[2] =  0
    g[2] = x[2] - 4.2*cs.sin(x[5])**2
    lbg[3] =  0
    ubg[3] =  0
    g[3] = x[0] + 2*x[1] + 2*x[2] - 7.2*cs.sin(x[6])**2
    obj = -x[0]*x[1]*x[2]
    
    x_opt = cs.DM([2.4, 1.2, 1.2, 0.857072, 2.57765, 2.57765, 14.1372])
    f = obj
    
    print('hs056: there are infinitely many solutions for the last 4 components of x.')
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
