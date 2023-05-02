# This is an automatically generated file converting from the apm format

import casadi as cs

def hs059():
    # The optimal objective is (if given in):
    f_opt = -7.804226324
    x_opt = cs.DM([13.5501, 51.66])
    x = cs.MX.sym('x', 2)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)
    u = cs.DM.zeros(2, 1)
    g = cs.MX.zeros(3, 1)
    lbg = -cs.inf*cs.DM.ones(3, 1)
    ubg = cs.inf*cs.DM.ones(3, 1)
    
    u[0] = 75
    u[1] = 65
    
    x0[0]= 90
    lbx[0]=0
    ubx[0]=u[0]
    x0[1]= 10
    lbx[1]=0
    ubx[1]=u[1]
    
    lbg[0] = 700
    g[0] = x[0]*x[1]
    lbg[1] = 0
    g[1] = x[1] - x[0]**2/125
    lbg[2] = 0
    g[2] = (x[1] - 50)**2 - 5*(x[0] - 55)
    obj = -75.196 + 3.8112*x[0] -0.12694*x[0]**2 + 0.0020567*x[0]**3 - 1.0345e-5*x[0]**4 + 6.8306*x[1] - 0.030234*x[0]*x[1] + 1.28134e-3*x[1]*x[0]**2 + 2.266e-7*x[0]**4*x[1] - 0.25645*x[1]**2 + 0.0034604*x[1]**3 - 1.3514e-5*x[1]**4 + 28.106/(x[1] + 1) + 5.2375e-6*x[0]**2*x[1]**2 + 6.3e-8*x[0]**3*x[1]**2 - 7e-10*x[0]**3*x[1]**3 - 3.405e-4*x[0]*x[1]**2 + 1.6638e-6*x[0]*x[1]**3 + 2.8673*cs.exp(0.0005*x[0]*x[1]) - 3.5256e-5*x[0]**3*x[1]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
