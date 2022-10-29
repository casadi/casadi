# This is an automatically generated file converting from the apm format
import casadi as cs


def hs118():
    # The optimal objective is (if given in):
    f_opt = 664.8204500
    x_opt = cs.DM([8, 49, 3, 1, 56, 0, 0.999999, 63, 6, 3, 70, 12, 5, 77, 18])
    x = cs.MX.sym('x', 15)
    obj = cs.MX.zeros(5, 1)
    x0 = cs.DM.zeros(15, 1)
    lbx = -cs.inf*cs.DM.ones(15, 1)
    ubx = cs.inf*cs.DM.ones(15, 1)
    g = cs.MX.zeros(17, 1)
    lbg = -cs.inf*cs.DM.ones(17, 1)
    ubg = cs.inf*cs.DM.ones(17, 1)

    x0[0] = 20
    lbx[0] = 8
    ubx[0] = 21
    x0[1] = 55
    lbx[1] = 43
    ubx[1] = 57
    x0[2] = 15
    lbx[2] = 3
    ubx[2] = 16
    x0[3] = 20
    lbx[3] = 0
    ubx[3] = 90
    x0[4] = 60
    lbx[4] = 0
    ubx[4] = 120
    x0[5] = 20
    lbx[5] = 0
    ubx[5] = 60
    x0[6] = 20
    lbx[6] = 0
    ubx[6] = 90
    x0[7] = 60
    lbx[7] = 0
    ubx[7] = 120
    x0[8] = 20
    lbx[8] = 0
    ubx[8] = 60
    x0[9] = 20
    lbx[9] = 0
    ubx[9] = 90
    x0[10] = 60
    lbx[10] = 0
    ubx[10] = 120
    x0[11] = 20
    lbx[11] = 0
    ubx[11] = 60
    x0[12] = 20
    lbx[12] = 0
    ubx[12] = 90
    x0[13] = 60
    lbx[13] = 0
    ubx[13] = 120
    x0[14] = 20
    lbx[14] = 0
    ubx[14] = 60

    g[0] = x[3] - x[0] + 7
    lbg[0] = 0
    ubg[0] = 13
    g[1] = x[6] - x[3] + 7
    lbg[1] = 0
    ubg[1] = 13
    g[2] = x[9] - x[6] + 7
    lbg[2] = 0
    ubg[2] = 13
    g[3] = x[12] - x[9] + 7
    lbg[3] = 0
    ubg[3] = 13
    g[4] = x[4] - x[1] + 7
    lbg[4] = 0
    ubg[4] = 14
    g[5] = x[7] - x[4] + 7
    lbg[5] = 0
    ubg[5] = 14
    g[6] = x[10] - x[7] + 7
    lbg[6] = 0
    ubg[6] = 14
    g[7] = x[13] - x[10] + 7
    lbg[7] = 0
    ubg[7] = 14
    g[8] = x[5] - x[2] + 7
    lbg[8] = 0
    ubg[8] = 13
    g[9] = x[8] - x[5] + 7
    lbg[9] = 0
    ubg[9] = 13
    g[10] = x[11] - x[8] + 7
    lbg[10] = 0
    ubg[10] = 13
    g[11] = x[14] - x[11] + 7
    lbg[11] = 0
    ubg[11] = 13
    lbg[12] = 60
    g[12] = x[0] + x[1] + x[2]
    lbg[13] = 50
    g[13] = x[3] + x[4] + x[5]
    lbg[14] = 70
    g[14] = x[6] + x[7] + x[8]
    lbg[15] = 85
    g[15] = x[9] + x[10] + x[11]
    lbg[16] = 100
    g[16] = x[12] + x[13] + x[14]
    obj[0] = (2.3*x[0] + 0.0001*x[0]**2 + 1.7*x[1] + 0.0001*x[1]**2 + 2.2*x[2] + 0.00015*x[2]**2)
    obj[1] = (2.3*x[3] + 0.0001*x[3]**2 + 1.7*x[4] + 0.0001*x[4]**2 + 2.2*x[5] + 0.00015*x[5]**2)
    obj[2] = (2.3*x[6] + 0.0001*x[6]**2 + 1.7*x[7] + 0.0001*x[7]**2 + 2.2*x[8] + 0.00015*x[8]**2)
    obj[3] = (2.3*x[9] + 0.0001*x[9]**2 + 1.7*x[10] + 0.0001*x[10]**2 + 2.2*x[11] + 0.00015*x[11]**2)
    obj[4] = (2.3*x[12] + 0.0001*x[12]**2 + 1.7*x[13] + 0.0001*x[13]**2 + 2.2*x[14] + 0.00015*x[14]**2)

    f = cs.sum1(obj)

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
