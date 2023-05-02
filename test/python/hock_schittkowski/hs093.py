# This is an automatically generated file converting from the apm format
import casadi as cs


def hs093():
    # The optimal objective is (if given in):
    f_opt = 135.075961
    x_opt = cs.DM([5.33267, 4.65674, 10.433, 12.0823, 0.752607, 0.878651])
    x = cs.MX.sym('x', 6)
    x0 = cs.DM.zeros(6, 1)
    lbx = -cs.inf*cs.DM.ones(6, 1)
    ubx = cs.inf*cs.DM.ones(6, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)

    x0[0] = 5.54
    x0[1] = 4.4
    x0[2] = 12.02
    x0[3] = 11.82
    x0[4] = 0.702
    x0[5] = 0.852

    lbg[0] = 2.07
    g[0] = 0.001 * x[0]*x[1]*x[2]*x[3]*x[4]*x[5]
    ubg[1] = 1
    g[1] = 0.00062*x[0]*x[3]*x[4]**2*(x[0] + x[1] + x[2]) + 0.00058*x[1]*x[2]*x[5]**2*(x[0] + 1.57*x[1] + x[3])
    f = 0.0204*x[0]*x[3]*(x[0] + x[1] + x[2]) + 0.0187*x[1]*x[2]*(x[0] + 1.57*x[1] + x[3]) + 0.0607*x[0]*x[3]*x[4]**2*(x[0] + x[1] + x[2]) + 0.0437*x[1]*x[2]*x[5]**2*(x[0] + 1.57*x[1] + x[3])

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
