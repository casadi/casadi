# This is an automatically generated file converting from the apm format

import casadi as cs

def hs083():
    # The optimal objective is (if given in):
    f_opt = -30665.53867
    x_opt = cs.DM([78, 33, 29.9953, 45, 36.7758])
    x = cs.MX.sym('x', 5)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(5, 1)
    lbx = -cs.inf*cs.DM.ones(5, 1)
    ubx = cs.inf*cs.DM.ones(5, 1)
    a = cs.DM.zeros(12, 1)
    l = cs.DM.zeros(5, 1)
    u = cs.DM.zeros(5, 1)
    xi = cs.DM.zeros(5, 1)
    g = cs.MX.zeros(3, 1)
    lbg = -cs.inf*cs.DM.ones(3, 1)
    ubg = cs.inf*cs.DM.ones(3, 1)
    
    a[0] = 85.334407
    a[1] = 0.0056858
    a[2] = 0.0006262
    a[3] = 0.0022053
    a[4] = 80.51249
    a[5] = 0.0071317
    a[6] = 0.0029955
    a[7] = 0.0021813
    a[8] = 9.300961
    a[9] = 0.0047026
    a[10] = 0.0012547
    a[11] = 0.0019085
    l[0] = 78
    l[1] = 33
    l[2] = 27
    l[3] = 27
    l[4] = 27
    u[0] = 102
    u[1] = 45
    u[2] = 45
    u[3] = 45
    u[4] = 45
    xi[0] = 78
    xi[1] = 33
    xi[2] = 27
    xi[3] = 27
    xi[4] = 27
    
    x0[0:5]= xi[0:5]
    lbx[0:5]=l[0:5]
    ubx[0:5]=u[0:5]
    
    g[0] =  a[0] + a[1]*x[1]*x[4] + a[2]*x[0]*x[3] - a[3]*x[2]*x[4]
    lbg[0] = 0
    ubg[0] =  92
    g[1] =  a[4] + a[5]*x[1]*x[4] + a[6]*x[0]*x[1] + a[7]*x[2]**2 
    lbg[1] = 90
    ubg[1] =  110
    g[2] =  a[8] + a[9]*x[2]*x[4] + a[10]*x[0]*x[2] + a[11]*x[2]*x[3] 
    lbg[2] = 20
    ubg[2] =  25
    obj = 5.3578547*x[2]**2 + 0.8356891*x[0]*x[4] + 37.293239*x[0] - 40792.141
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
