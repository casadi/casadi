# This is an automatically generated file converting from the apm format
import casadi as cs


def hs106():
    # The optimal objective is (if given in):
    f_opt = 7049.24776020595
    x_opt = cs.DM([579.307, 1359.97, 5109.97, 182.018, 295.601, 217.982, 286.417, 395.601])
    x = cs.MX.sym('x', 8)
    x0 = cs.DM.zeros(8, 1)
    lbx = -cs.inf*cs.DM.ones(8, 1)
    ubx = cs.inf*cs.DM.ones(8, 1)
    g = cs.MX.zeros(6, 1)
    lbg = -cs.inf*cs.DM.ones(6, 1)
    ubg = cs.inf*cs.DM.ones(6, 1)

    a = 0.0025
    b = 0.01
    c = 833.3325
    d = 100
    e = 83333.33
    f = 1250
    k = 1250000
    h = 2500

    x0[0] = 5000
    x0[1] = 5000
    x0[2] = 5000
    x0[3] = 200
    x0[4] = 350
    x0[5] = 150
    x0[6] = 225
    x0[7] = 425

    lbg[0] = -1
    g[0] = - a * (x[3] + x[5])
    lbg[1] = -1
    g[1] = - a * (x[4] + x[6] - x[3])
    lbg[2] = -1
    g[2] = - b * (x[7] - x[4])
    lbg[3] = -e
    g[3] = x[0] * x[5] - c * x[3] - d * x[0]
    lbg[4] = 0
    g[4] = x[1] * x[6] - f * x[4] - x[1] * x[3] + f * x[3]
    lbg[5] = k
    g[5] = x[2] * x[7] - x[2] * x[4] + h * x[4]
    lbx[0] = 100
    ubx[0] = 10000
    lbx[1:3] = 1000
    ubx[1:3] = 10000
    lbx[3:8] = 10
    ubx[3:8] = 1000
    f = x[0] + x[1] + x[2]

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
