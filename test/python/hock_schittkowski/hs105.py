# This is an automatically generated file converting from the apm format
import casadi as cs


def hs105():
    # The optimal objective is (if given in):
    f_opt = 1136.36
    x_opt = cs.DM([0.393749, 0.449, 130.984, 164.653, 221.471, 12.2415, 17.7368, 17.0334])
    x = cs.MX.sym('x', 8)
    obj = cs.MX.zeros(235, 1)
    x0 = cs.DM.zeros(8, 1)
    lbx = -cs.inf*cs.DM.ones(8, 1)
    ubx = cs.inf*cs.DM.ones(8, 1)
    y = cs.DM.zeros(235, 1)
    a = cs.MX.zeros(235, 1)
    b = cs.MX.zeros(235, 1)
    c = cs.MX.zeros(235, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)

    y[0] = 95
    y[1] = 105
    y[2:6] = 110
    y[6:10] = 115
    y[10:25] = 120
    y[25:40] = 125
    y[40:55] = 130
    y[55:68] = 135
    y[68:89] = 140
    y[89:101] = 145
    y[101:118] = 150
    y[118:122] = 155
    y[122:142] = 160
    y[142:150] = 165
    y[150:167] = 170
    y[167:175] = 175
    y[175:181] = 180
    y[181:187] = 185
    y[187:194] = 190
    y[194:198] = 195
    y[198:201] = 200
    y[201:204] = 205
    y[204:212] = 210
    y[212] = 215
    y[213:219] = 220
    y[219:224] = 230
    y[224] = 235
    y[225:232] = 240
    y[232] = 245
    y[233:235] = 250

    x0[0] = .1
    x0[1] = .2
    x0[2] = 100
    x0[3] = 125
    x0[4] = 175
    x0[5] = 11.2
    x0[6] = 13.2
    x0[7] = 15.8

    a[0:235] = x[0] / x[5] * cs.exp(-(y[0:235] - x[2])**2 / (2 * x[5]**2))
    b[0:235] = x[1] / x[6] * cs.exp(-(y[0:235] - x[3])**2 / (2 * x[6]**2))
    c[0:235] = (1 - x[1] - x[0]) / x[7] * cs.exp(-(y[0:235] - x[4])**2 / (2 * x[7]**2))

    lbg[0] = -1
    g[0] = - x[0] - x[1]
    lbx[0] = .001
    ubx[0] = .499
    lbx[1] = .001
    ubx[1] = .449
    lbx[2] = 100
    ubx[2] = 180
    lbx[3] = 130
    ubx[3] = 210
    lbx[4] = 170
    ubx[4] = 240
    lbx[5] = 5
    ubx[5] = 25
    lbx[6] = 5
    ubx[6] = 25
    lbx[7] = 5
    ubx[7] = 25
    obj[0:235] = - (cs.log((a[0:235] + b[0:235] + c[0:235]) / cs.sqrt(2 * cs.pi)))

    f = cs.sum1(obj)

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
