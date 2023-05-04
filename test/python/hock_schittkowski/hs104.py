# This is an automatically generated file converting from the apm format
import casadi as cs


def hs104():
    # The optimal objective is (if given in):
    f_opt = 3.9511634396
    x_opt = cs.DM([6.46511, 2.23271, 0.667397, 0.595756, 5.93268, 5.52723, 1.01332, 0.400668])
    x = cs.MX.sym('x', 8)
    x0 = cs.DM.zeros(8, 1)
    lbx = -cs.inf*cs.DM.ones(8, 1)
    ubx = cs.inf*cs.DM.ones(8, 1)
    g = cs.MX.zeros(5, 1)
    lbg = -cs.inf*cs.DM.ones(5, 1)
    ubg = cs.inf*cs.DM.ones(5, 1)

    x0[0] = 6
    lbx[0] = .1
    ubx[0] = 10
    x0[1] = 3
    lbx[1] = .1
    ubx[1] = 10
    x0[2] = 0.4
    lbx[2] = .1
    ubx[2] = 10
    x0[3] = 0.2
    lbx[3] = .1
    ubx[3] = 10
    x0[4] = 6
    lbx[4] = .1
    ubx[4] = 10
    x0[5] = 6
    lbx[5] = .1
    ubx[5] = 10
    x0[6] = 1
    lbx[6] = .1
    ubx[6] = 10
    x0[7] = 0.5
    lbx[7] = .1
    ubx[7] = 10

    lbg[0] = -1
    g[0] = -.0588*x[4]*x[6]-.1*x[0]
    lbg[1] = -1
    g[1] = -.0588*x[5]*x[7]-.1*x[0]-.1*x[1]
    lbg[2] = -1
    g[2] = -4*x[2]/x[4]-2/(x[2]**.71*x[4])-.0588*x[6]/x[2]**1.3
    lbg[3] = -1
    g[3] = -4*x[3]/x[5]-2/(x[3]**.71*x[5])-.0588*x[7]/x[3]**1.3
    lbg[4] = .1
    g[4] = .4*x[0]**.67*x[6]**(-.67)+.4*x[1]**.67*x[7]**(-.67)+10-x[0]-x[1]
    ubg[4] = 4.2
    f = .4*x[0]**.67*x[6]**(-.67)+.4*x[1]**.67*x[7]**(-.67)+10-x[0]-x[1]

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
