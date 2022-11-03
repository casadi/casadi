# This is an automatically generated file converting from the apm format

import casadi as cs

def hs080():
    # The optimal objective is (if given in):
    f_opt = 0.0539498478
    x_opt = cs.DM([-1.71714, 1.59571, 1.82725, -0.763643, -0.763643])
    x = cs.MX.sym('x', 5)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(5, 1)
    lbx = -cs.inf*cs.DM.ones(5, 1)
    ubx = cs.inf*cs.DM.ones(5, 1)
    l = cs.DM.zeros(5, 1)
    u = cs.DM.zeros(5, 1)
    xi = cs.DM.zeros(5, 1)
    g = cs.MX.zeros(3, 1)
    lbg = -cs.inf*cs.DM.ones(3, 1)
    ubg = cs.inf*cs.DM.ones(3, 1)
    
    l[0] = -2.3
    l[1] = -2.3
    l[2] = -3.2
    l[3] = -3.2
    l[4] = -3.2
    u[0] = 2.3
    u[1] = 2.3
    u[2] = 3.2
    u[3] = 3.2
    u[4] = 3.2
    xi[0] = -2
    xi[1] =  2
    xi[2] =  2
    xi[3] = -1
    xi[4] = -1
    
    x0[0:5]= xi[0:5]
    lbx[0:5]=l[0:5]
    ubx[0:5]=u[0:5]
    
    lbg[0] =  10
    ubg[0] =  10
    g[0] = x[0]**2+x[1]**2+x[2]**2+x[3]**2+x[4]**2
    lbg[1] =  0
    ubg[1] =  0
    g[1] = x[1]*x[2] - 5*x[3]*x[4]
    lbg[2] =  -1
    ubg[2] =  -1
    g[2] = x[0]**3 + x[1]**3
    obj = cs.exp( x[0]*x[1]*x[2]*x[3]*x[4] )
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
