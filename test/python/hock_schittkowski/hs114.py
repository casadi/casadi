# This is an automatically generated file converting from the apm format
import casadi as cs


def hs114():
    # The optimal objective is (if given in):
    f_opt = -1768.80696
    x_opt = cs.DM([1698.09, 15818.6, 54.1027, 3031.23, 2000, 90.1154, 95, 10.4933, 1.56164, 153.535])
    x = cs.MX.sym('x', 10)
    x0 = cs.DM.zeros(10, 1)
    lbx = -cs.inf*cs.DM.ones(10, 1)
    ubx = cs.inf*cs.DM.ones(10, 1)
    x0 = cs.DM.zeros(10, 1)
    lb = cs.DM.zeros(10, 1)
    ub = cs.DM.zeros(10, 1)
    g = cs.MX.zeros(11, 1)
    lbg = -cs.inf*cs.DM.ones(11, 1)
    ubg = cs.inf*cs.DM.ones(11, 1)

    a = .99
    b = .9
    x0[0] = 1745
    x0[1] = 12000
    x0[2] = 110
    x0[3] = 3048
    x0[4] = 1974
    x0[5] = 89.2
    x0[6] = 92.8
    x0[7] = 8
    x0[8] = 3.6
    x0[9] = 145
    lb[0] = .00001
    lb[1] = .00001
    lb[2] = .00001
    lb[3] = .00001
    lb[4] = .00001
    lb[5] = 85
    lb[6] = 90
    lb[7] = 3
    lb[8] = 1.2
    lb[9] = 145
    ub[0] = 2000
    ub[1] = 16000
    ub[2] = 120
    ub[3] = 5000
    ub[4] = 2000
    ub[5] = 93
    ub[6] = 95
    ub[7] = 12
    ub[8] = 4
    ub[9] = 162

    x0[0:10] = x0[0:10]
    lbx[0:10] = lb[0:10]
    ubx[0:10] = ub[0:10]

    G1 = 35.82 - .222*x[9] - b*x[8]
    G2 = -133 + 3*x[6] - a*x[9]
    G5 = 1.12*x[0] + .13167*x[0]*x[7] - .00667*x[0]*x[7]**2 - a*x[3]
    G6 = 57.425 + 1.098*x[7] - .038*x[7]**2 + .325*x[5] - a*x[6]

    lbg[0] = 0
    g[0] = G1
    lbg[1] = 0
    g[1] = G2
    lbg[2] = 0
    g[2] = -G1 + x[8]*(1/b - b)
    lbg[3] = 0
    g[3] = -G2 + (1/a - a)*x[9]
    lbg[4] = 0
    g[4] = G5
    lbg[5] = 0
    g[5] = G6
    lbg[6] = 0
    g[6] = -G5 + (1/a - a)*x[3]
    lbg[7] = 0
    g[7] = -G6 + (1/a - a)*x[6]
    lbg[8] = 0
    ubg[8] = 0
    g[8] = 1.22*x[3] - x[0] - x[4]
    lbg[9] = 0
    ubg[9] = 0
    g[9] = 98000*x[2]/(x[3]*x[8] + 1000*x[2]) - x[5]
    lbg[10] = 0
    ubg[10] = 0
    g[10] = (x[1] + x[4])/x[0] - x[7]
    f = 5.04*x[0] + .035*x[1] + 10*x[2] + 3.36*x[4] - .063*x[3]*x[6]

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
