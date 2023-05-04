# This is an automatically generated file converting from the apm format
import casadi as cs


def hs111():
    # The optimal objective is (if given in):
    f_opt = -47.76109026
    x_opt = cs.DM([-3.20231, -1.91237, -0.244427, -6.56118, -0.723098, -7.27423, -3.59724, -4.02032, -3.28838, -2.33437])
    x = cs.MX.sym('x', 10)
    obj = cs.MX.zeros(10, 1)
    x0 = cs.DM.zeros(10, 1)
    lbx = -cs.inf*cs.DM.ones(10, 1)
    ubx = cs.inf*cs.DM.ones(10, 1)
    c = cs.DM.zeros(10, 1)
    sum = cs.MX.zeros(11, 1)
    g = cs.MX.zeros(3, 1)
    lbg = -cs.inf*cs.DM.ones(3, 1)
    ubg = cs.inf*cs.DM.ones(3, 1)

    c[0] = -6.089
    c[1] = -17.164
    c[2] = -34.054
    c[3] = -5.914
    c[4] = -24.721
    c[5] = -14.986
    c[6] = -24.100
    c[7] = -10.708
    c[8] = -26.662
    c[9] = -22.179

    x0[0:10] = -2.3
    lbx[0:10] = -100
    ubx[0:10] = 100

    sum[0:10] = cs.exp(x[0:10])
    sum[10] = cs.sum1(sum[0:10])

    lbg[0] = 2
    ubg[0] = 2
    g[0] = cs.exp(x[0]) + 2*cs.exp(x[1]) + 2*cs.exp(x[2]) + cs.exp(x[5]) + cs.exp(x[9])
    lbg[1] = 1
    ubg[1] = 1
    g[1] = cs.exp(x[3]) + 2*cs.exp(x[4]) + cs.exp(x[5]) + cs.exp(x[6])
    lbg[2] = 1
    ubg[2] = 1
    g[2] = cs.exp(x[2]) + cs.exp(x[6]) + cs.exp(x[7]) + 2*cs.exp(x[8]) + cs.exp(x[9])
    obj[0:10] = cs.exp(x[0:10])*(c[0:10] + x[0:10] - cs.log(sum[10]))

    f = cs.sum1(obj)

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
