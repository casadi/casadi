# This is an automatically generated file converting from the apm format
import casadi as cs


def hs108():
    # The optimal objective is (if given in):
    f_opt = -.8660254038
    x_opt = cs.DM([-0.554603, -0.831898, 0.44323, -0.896408, -0.554697, -0.832052, 0.443331, -0.896141, 0.0002167])
    x = cs.MX.sym('x', 9)
    x0 = cs.DM.zeros(9, 1)
    lbx = -cs.inf*cs.DM.ones(9, 1)
    ubx = cs.inf*cs.DM.ones(9, 1)
    g = cs.MX.zeros(13, 1)
    lbg = -cs.inf*cs.DM.ones(13, 1)
    ubg = cs.inf*cs.DM.ones(13, 1)

    x0[0:9] = 1

    lbg[0] = -1
    g[0] = -x[2]**2-x[3]**2
    lbg[1] = -1
    g[1] = -x[4]**2-x[5]**2
    lbg[2] = -1
    g[2] = -x[8]**2
    lbg[3] = -1
    g[3] = -x[0]**2-(x[1]-x[8])**2
    lbg[4] = -1
    g[4] = -(x[0]-x[4])**2-(x[1]-x[5])**2
    lbg[5] = -1
    g[5] = -(x[0]-x[6])**2-(x[1]-x[7])**2
    lbg[6] = -1
    g[6] = -(x[2]-x[6])**2-(x[3]-x[7])**2
    lbg[7] = -1
    g[7] = -(x[2]-x[4])**2-(x[3]-x[5])**2
    lbg[8] = -1
    g[8] = -x[6]**2-(x[7]-x[8])**2
    lbg[9] = 0
    g[9] = x[0]*x[3]-x[1]*x[2]
    lbg[10] = 0
    g[10] = x[2]*x[8]
    lbg[11] = 0
    g[11] = -x[4]*x[8]
    lbg[12] = 0
    g[12] = x[4]*x[7]-x[5]*x[6]
    lbx[8] = 0
    f = -.5*(x[0]*x[3]-x[1]*x[2]+x[2]*x[8]-x[4]*x[8]+x[4]*x[7]-x[5]*x[6]) + (3*x[0] - 2*x[1])**2 + (3*x[4] - 2*x[5])**2

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
