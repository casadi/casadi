# This is an automatically generated file converting from the apm format

import casadi as cs

def hs069():
    # The optimal objective is (if given in):
    f_opt = -956.71288
    x_opt = cs.DM([0.0293714, 1.19025, 0.233947, 0.791668])
    x = cs.MX.sym('x', 4)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(4, 1)
    lbx = -cs.inf*cs.DM.ones(4, 1)
    ubx = cs.inf*cs.DM.ones(4, 1)
    l = cs.DM.zeros(4, 1)
    u = cs.DM.zeros(4, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)
    
    a = 0.1
    b = 1000
    d = 1
    n = 4
    l[0] = 0.0001
    l[1] = 0
    l[2] = 0
    l[3] = 0
    u[0] = 100
    u[1] = 100
    u[2] = 2
    u[3] = 2
    
    x0[0:4]= 1
    lbx[0:4]=l[0:4]
    ubx[0:4]=u[0:4]
    
    argn = -1*x[1] - d*cs.sqrt(n)
    arg0 = -1*x[1]
    argp = -1*x[1] + d*cs.sqrt(n)
    phin = (1/2)+(1/2)*cs.erf(argn/cs.sqrt(2))
    phi0 = (1/2)+(1/2)*cs.erf(arg0/cs.sqrt(2))
    phip = (1/2)+(1/2)*cs.erf(argp/cs.sqrt(2))

    lbg[0] =  0
    ubg[0] =  0
    g[0] = x[2] - 2*phi0
    lbg[1] =  0
    ubg[1] =  0
    g[1] = x[3] - phip - phin
    obj = ( a*n - (b*(cs.exp(x[0])-1) - x[2])*x[3]/(cs.exp(x[0]) - 1 + x[3]) )/x[0]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
