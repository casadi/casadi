# This is an automatically generated file converting from the apm format
import casadi as cs


def hs116():
    # The optimal objective is (if given in):
    f_opt = 97.588409
    x_opt = cs.DM([0.803773, 0.899986, 0.970982, 0.1, 0.190813, 0.460655, 574.077, 74.0775, 500.016, 0.1, 20.233, 77.3477, 0.00672893])
    x = cs.MX.sym('x', 13)
    x0 = cs.DM.zeros(13, 1)
    lbx = -cs.inf*cs.DM.ones(13, 1)
    ubx = cs.inf*cs.DM.ones(13, 1)
    g = cs.MX.zeros(15, 1)
    lbg = -cs.inf*cs.DM.ones(15, 1)
    ubg = cs.inf*cs.DM.ones(15, 1)

    a = 0.002
    b = 1.262626
    c = 1.231059
    d = 0.03475
    e = 0.975
    f = 0.00975

    x0[0] = 0.5
    lbx[0] = 0
    x0[1] = 0.8
    lbx[1] = 0
    x0[2] = 0.9
    lbx[2] = 0
    x0[3] = 0.1
    lbx[3] = 0
    x0[4] = 0.14
    lbx[4] = 0
    x0[5] = 0.5
    lbx[5] = 0
    x0[6] = 489
    lbx[6] = 0
    x0[7] = 80
    lbx[7] = 0
    x0[8] = 650
    lbx[8] = 0
    x0[9] = 450
    lbx[9] = 0
    x0[10] = 150
    lbx[10] = 0
    x0[11] = 150
    lbx[11] = 0
    x0[12] = 150
    lbx[12] = 0

    lbg[0] = 0
    g[0] = x[2] - x[1]
    lbg[1] = 0
    g[1] = x[1] - x[0]
    lbg[2] = 0
    g[2] = 1 - a * x[6] + a * x[7]
    lbg[3] = 50
    g[3] = x[10] + x[11] + x[12]
    lbg[4] = 0
    g[4] = x[12] - b * x[9] + c * x[2] * x[9]
    lbg[5] = 0
    g[5] = x[4] - d * x[1] - e * x[1] * x[4] + f * x[1]**2
    lbg[6] = 0
    g[6] = x[5] - d * x[2] - e * x[2] * x[5] + f * x[2]**2
    lbg[7] = 0
    g[7] = x[3] - d * x[0] - e * x[0] * x[3] + f * x[0]**2
    lbg[8] = 0
    g[8] = x[11] - b * x[8] + c * x[1] * x[8]
    lbg[9] = 0
    g[9] = x[10] - b * x[7] + c * x[0] * x[7]
    lbg[10] = 0
    g[10] = x[4] * x[6] - x[0] * x[7] - x[3] * x[6] + x[3] * x[7]
    lbg[11] = 0
    g[11] = 1 - a * (x[1] * x[8] + x[4] * x[7] - x[0] * x[7] - x[5] * x[8]) - x[4] - x[5]
    lbg[12] = 0
    g[12] = x[1] * x[8] - x[2] * x[9] - x[5] * x[8] - 500 * x[1] + 500 * x[5] + x[1] * x[9]
    lbg[13] = 0
    g[13] = x[1] - 0.9 - a * (x[1] * x[9] - x[2] * x[9])
    ubg[14] = 250
    g[14] = x[10] + x[11] + x[12]
    lbx[0] = 0.1
    ubx[0] = 1
    lbx[1] = 0.1
    ubx[1] = 1
    lbx[2] = 0.1
    ubx[2] = 1
    lbx[3] = 0.0001
    ubx[3] = 0.1
    lbx[4] = 0.1
    ubx[4] = 0.9
    lbx[5] = 0.1
    ubx[5] = 0.9
    lbx[6] = 0.1
    ubx[6] = 1000
    lbx[7] = 0.1
    ubx[7] = 1000
    lbx[8] = 500
    ubx[8] = 1000
    lbx[9] = 0.1
    ubx[9] = 500
    lbx[10] = 1
    ubx[10] = 150
    lbx[11] = 0.0001
    ubx[11] = 150
    lbx[12] = 0.0001
    ubx[12] = 150
    f = x[10] + x[11] + x[12]

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
