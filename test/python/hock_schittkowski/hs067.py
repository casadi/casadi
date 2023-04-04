# This is an automatically generated file converting from the apm format

import casadi as cs

def hs067():
    # The optimal objective is (if given in):
    f_opt = -1162.02698006
    x_opt = cs.DM([1728.37, 16000, 98.1362, 3056.04, 2000, 90.6185, 94.1896, 10.4144, 2.61574, 149.569])
    x = cs.MX.sym('x', 10)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(10, 1)
    lbx = -cs.inf*cs.DM.ones(10, 1)
    ubx = cs.inf*cs.DM.ones(10, 1)
    g = cs.MX.zeros(7, 1)
    lbg = -cs.inf*cs.DM.ones(7, 1)
    ubg = cs.inf*cs.DM.ones(7, 1)
    
    x0[0]= 1745
    lbx[0]=1e-5
    ubx[0]= 2000
    x0[1]= 12000
    lbx[1]=1e-5
    ubx[1]= 16000
    x0[2]= 110
    lbx[2]=1e-5
    ubx[2]= 120
    
    lbx[3] = 0
    ubx[3] =  5000
    lbx[4] = 0
    ubx[4] =  2000
    lbx[5] = 85
    ubx[5] =  93
    lbx[6] = 90
    ubx[6] =  95
    lbx[7] = 3
    ubx[7] =  12
    lbx[8] = 0.01
    ubx[8] =  4
    lbx[9] = 145
    ubx[9] =  162
    lbg[0] = 0
    ubg[0] = 0
    g[0] = x[4]  - ( 1.22*x[3] - x[0])
    lbg[1] = 0
    ubg[1] = 0
    g[1] = x[7]  - ( (x[1]+x[4])/x[0])
    lbg[2] = 0
    ubg[2] = 0
    g[2] = x[3]  - ( 0.01*x[0]*(112 + 13.167*x[7] - 0.6667*x[7]**2))
    lbg[3] = 0
    ubg[3] = 0
    g[3] = x[6]  - ( 86.35 + 1.098*x[7] - 0.038*x[7]**2 + 0.325*(x[5]-89))
    lbg[4] = 0
    ubg[4] = 0
    g[4] = x[9]  - ( 3*x[6] - 133)
    lbg[5] = 0
    ubg[5] = 0
    g[5] = x[8]  - ( 35.82 - 0.222*x[9])
    lbg[6] = 0
    ubg[6] = 0
    g[6] = x[5]  - ( 98000*x[2]/(x[3]*x[8] + 1000*x[2]))
    obj = -(0.063*x[3]*x[6] - 5.04*x[0] - 3.36*x[4] - 0.035*x[1] - 10*x[2])
    
    
    f = obj
    print('hs067: Currently this is an alternative problem formulation.')
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
