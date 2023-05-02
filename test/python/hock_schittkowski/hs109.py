# This is an automatically generated file converting from the apm format
import casadi as cs


def hs109():
    # The optimal objective is (if given in):
    x_opt = cs.DM([675.025, 1134.02, 0.133485, -0.37119, 252, 252, 201.466, 426.619, 368.488])
    f_opt = 5362.06928
    x = cs.MX.sym('x', 9)
    x0 = cs.DM.zeros(9, 1)
    lbx = -cs.inf*cs.DM.ones(9, 1)
    ubx = cs.inf*cs.DM.ones(9, 1)
    l = cs.DM.zeros(9, 1)
    u = cs.DM.zeros(9, 1)
    g = cs.MX.zeros(10, 1)
    lbg = -cs.inf*cs.DM.ones(10, 1)
    ubg = cs.inf*cs.DM.ones(10, 1)

    a = 50.176
    b1 = .25
    b = cs.sin(b1)
    c = cs.cos(b1)
    l[0:2] = 0
    l[2:4] = -0.55
    l[4:7] = 196
    l[7:9] = -400
    u[0:2] = cs.inf
    u[2:4] = 0.55
    u[4:7] = 252
    u[7:9] = 800

    x0[0:9] = 0
    lbx[0:9] = l[0:9]
    ubx[0:9] = u[0:9]

    lbg[0] = 0
    g[0] = x[3] - x[2] + .55
    lbg[1] = 0
    g[1] = x[2] - x[3] + .55
    lbg[2] = 0
    g[2] = 2250000 - x[0]**2 - x[7]**2
    lbg[3] = 0
    g[3] = 2250000 - x[1]**2 - x[8]**2

    lbg[4] = 0
    ubg[4] = 0
    g[4] = x[4] * x[5] * cs.sin(-x[2] - .25) + x[4] * x[6] * cs.sin(-x[3] - .25) + 2 * b * x[4]**2 - a * x[0] + 400 * a
    lbg[5] = 0
    ubg[5] = 0
    g[5] = x[4] * x[5] * cs.sin(x[2] - .25) + x[5] * x[6] * cs.sin(x[2] - x[3] - .25) + 2 * b * x[5]**2 - a * x[1] + 400 * a
    lbg[6] = 0
    ubg[6] = 0
    g[6] = x[4] * x[6] * cs.sin(x[3] - .25) + x[5] * x[6] * cs.sin(x[3] - x[2] - .25) + 2 * b * x[6]**2 + 881.779 * a
    lbg[7] = 0
    ubg[7] = 0
    g[7] = a * x[7] + x[4] * x[5] * cs.cos(-x[2] - .25) + x[4] * x[6] * cs.cos(-x[3] - .25) - 200 * a - 2 * c * x[4]**2 + .7533e-3 * a * x[4]**2

    lbg[8] = 0
    ubg[8] = 0
    g[8] = a * x[8] + x[4] * x[5] * cs.cos(x[2] - .25) + x[5] * x[6] * cs.cos(x[2] - x[3] - .25) - 2 * c * x[5]**2 + .7533e-3 * a * x[5]**2 - 200 * a
    lbg[9] = 0
    ubg[9] = 0
    g[9] = x[4] * x[6] * cs.cos(x[3] - .25) + x[5] * x[6] * cs.cos(x[3] - x[2] - .25) - 2 * c * x[6]**2 - 22.938 * a + .7533e-3 * a * x[6] **2

    f = 3 * x[0] + 1e-6 * x[0]**3 + 2 * x[1] + .522074e-6 * x[1]**3

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
