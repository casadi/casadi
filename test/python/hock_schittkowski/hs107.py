# This is an automatically generated file converting from the apm format
import casadi as cs


def hs107():
    # The optimal objective is (if given in):
    f_opt = 5055.011803
    x_opt = cs.DM([0.667013, 1.02238, 0.228287, 0.184822, 1.0909, 1.0909, 1.06904, 0.106611, -0.338788])
    x = cs.MX.sym('x', 9)
    x0 = cs.DM.zeros(9, 1)
    lbx = -cs.inf*cs.DM.ones(9, 1)
    ubx = cs.inf*cs.DM.ones(9, 1)
    g = cs.MX.zeros(6, 1)
    lbg = -cs.inf*cs.DM.ones(6, 1)
    ubg = cs.inf*cs.DM.ones(6, 1)

    c = (48.4/50.176)*cs.sin(.25)
    d = (48.4/50.176)*cs.cos(.25)

    x0[0] = .8
    x0[1] = .8
    x0[2] = .2
    x0[3] = .2
    x0[4] = 1.0454
    x0[5] = 1.0454
    x0[6] = 0
    x0[7] = 0
    x0[8] = 0

    y1 = cs.sin(x[7])
    y2 = cs.cos(x[7])
    y3 = cs.sin(x[8])
    y4 = cs.cos(x[8])
    y5 = cs.sin(x[7]-x[8])
    y6 = cs.cos(x[7]-x[8])

    lbg[0] = -0.4
    ubg[0] = -0.4
    g[0] = -x[0]+2*c*x[4]**2-x[4]*x[5]*(d*y1+c*y2)-x[4]*x[6]*(d*y3+c*y4)
    lbg[1] = -0.4
    ubg[1] = -0.4
    g[1] = -x[1]+2*c*x[5]**2+x[4]*x[5]*(d*y1-c*y2)+x[5]*x[6]*(d*y5-c*y6)
    lbg[2] = -.8
    ubg[2] = -.8
    g[2] = 2*c*x[6]**2+x[4]*x[6]*(d*y3-c*y4)-x[5]*x[6]*(d*y5+c*y6)
    lbg[3] = -0.2
    ubg[3] = -0.2
    g[3] = -x[2]+2*d*x[4]**2+x[4]*x[5]*(c*y1-d*y2)+x[4]*x[6]*(c*y3-d*y4)
    lbg[4] = -0.2
    ubg[4] = -0.2
    g[4] = -x[3]+2*d*x[5]**2-x[4]*x[5]*(c*y1+d*y2)-x[5]*x[6]*(c*y5+d*y6)
    lbg[5] = .337
    ubg[5] = .337
    g[5] = 2*d*x[6]**2-x[4]*x[6]*(c*y3+d*y4)+x[5]*x[6]*(c*y5-d*y6)
    lbx[0] = 0
    lbx[1] = 0
    lbx[4] = .90909
    lbx[5] = .90909
    lbx[6] = .90909
    ubx[4] = 1.0909
    ubx[5] = 1.0909
    ubx[6] = 1.0909
    f = 3000*x[0]+1000*x[0]**3+2000*x[1]+666.667*x[1]**3

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
