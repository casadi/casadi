# This is an automatically generated file converting from the apm format
import casadi as cs


def hs112():
    # The optimal objective is (if given in):
    f_opt = -47.76109026
    x_opt = cs.DM([0.0406681, 0.14773, 0.783153, 0.00141422, 0.485247, 0.000693175, 0.0273993, 0.0179473, 0.0373144, 0.0968713])
    x = cs.MX.sym('x', 10)
    obj = cs.MX.zeros(10, 1)
    x0 = cs.DM.zeros(10, 1)
    lbx = -cs.inf*cs.DM.ones(10, 1)
    ubx = cs.inf*cs.DM.ones(10, 1)
    c = cs.DM.zeros(10, 1)
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

    x0[0:10] = 0.1
    lbx[0:10] = 1.0e-6

    sum = cs.sum1(x)

    lbg[0] = 2
    ubg[0] = 2
    g[0] = x[0] + 2*x[1] + 2*x[2] + x[5] + x[9]
    lbg[1] = 1
    ubg[1] = 1
    g[1] = x[3] + 2*x[4] + x[5] + x[6]
    lbg[2] = 1
    ubg[2] = 1
    g[2] = x[2] + x[6] + x[7] + 2*x[8] + x[9]
    obj[0:10] = x[0:10] * (c[0:10] + cs.log(x[0:10]/sum))

    f = cs.sum1(obj)

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
