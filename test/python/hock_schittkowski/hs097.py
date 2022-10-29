# This is an automatically generated file converting from the apm format
import casadi as cs


def hs097():
    # The optimal objective is (if given in):
    f_opt = 3.1358091
    x_opt = cs.DM([0.268564, 0, 0, 0, 0.028, 0.0134])
    x = cs.MX.sym('x', 6)
    x0 = cs.DM.zeros(6, 1)
    lbx = -cs.inf*cs.DM.ones(6, 1)
    ubx = cs.inf*cs.DM.ones(6, 1)
    u = cs.DM.zeros(6, 1)
    g = cs.MX.zeros(4, 1)
    lbg = -cs.inf*cs.DM.ones(4, 1)
    ubg = cs.inf*cs.DM.ones(4, 1)

    u[0] = 0.31
    u[1] = 0.046
    u[2] = 0.068
    u[3] = 0.042
    u[4] = 0.028
    u[5] = 0.0134

    x0[0:6] = 0
    lbx[0:6] = 0
    ubx[0:6] = u[0:6]

    lbg[0] = 32.97
    g[0] = 17.1*x[0] + 38.2*x[1] + 204.2*x[2] + 212.3*x[3] + 623.4*x[4] + 1495.5*x[5] - 169*x[0]*x[2] - 3580*x[2]*x[4] - 3810*x[3]*x[4] - 18500*x[3]*x[5] - 24300*x[4]*x[5]
    lbg[1] = 25.12
    g[1] = 17.9*x[0] + 36.8*x[1] + 113.9*x[2] + 169.7*x[3] + 337.8*x[4] + 1385.2*x[5] - 139*x[0]*x[2] - 2450*x[3]*x[4] - 16600*x[3]*x[5] - 17200*x[4]*x[5]
    lbg[2] = -29.08
    g[2] = -273*x[1] - 70*x[3] - 819*x[4] + 26000*x[3]*x[4]
    lbg[3] = -78.02
    g[3] = 159.9*x[0] - 311*x[1] + 587*x[3] + 391*x[4] + 2198*x[5] - 14000*x[0]*x[5]
    f = 4.3*x[0] + 31.8*x[1] + 63.3*x[2] + 15.8*x[3] + 68.5*x[4] + 4.7*x[5]

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
