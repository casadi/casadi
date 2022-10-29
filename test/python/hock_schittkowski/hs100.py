# This is an automatically generated file converting from the apm format
import casadi as cs


def hs100():
    # The optimal objective is (if given in):
    f_opt = 680.6300573
    x_opt = cs.DM([2.3305, 1.95137, -0.477541, 4.36573, -0.624487, 1.03813, 1.59423])
    x = cs.MX.sym('x', 7)
    x0 = cs.DM.zeros(7, 1)
    lbx = -cs.inf*cs.DM.ones(7, 1)
    ubx = cs.inf*cs.DM.ones(7, 1)
    g = cs.MX.zeros(4, 1)
    lbg = -cs.inf*cs.DM.ones(4, 1)
    ubg = cs.inf*cs.DM.ones(4, 1)

    x0[0] = 1
    x0[1] = 2
    x0[2] = 0
    x0[3] = 4
    x0[4] = 0
    x0[5] = 1
    x0[6] = 1

    ubg[0] = 127
    g[0] = 2*x[0]**2 + 3*x[1]**4 + x[2] + 4*x[3]**2 + 5*x[4]
    ubg[1] = 282
    g[1] = 7*x[0] + 3*x[1] + 10*x[2]**2 + x[3] - x[4]
    ubg[2] = 196
    g[2] = 23*x[0] + x[1]**2 + 6*x[5]**2 - 8*x[6]
    lbg[3] = 0
    g[3] = -4*x[0]**2 - x[1]**2 + 3*x[0]*x[1] - 2*x[2]**2 - 5*x[5] + 11*x[6]
    f = (x[0]-10)**2 + 5*(x[1]-12)**2 + x[2]**4 + 3*(x[3]-11)**2 + 10*x[4]**6 + 7*x[5]**2 + x[6]**4 - 4*x[5]*x[6] - 10*x[5] - 8*x[6]

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
