# This is an automatically generated file converting from the apm format
import casadi as cs


def hs070():
    # The optimal objective is (if given in):
    f_opt = 0.007498464
    x_opt = cs.DM([12.27695, 4.631788, 0.3128625, 2.029290])
    x = cs.MX.sym('x', 4)
    obj = cs.MX.zeros(19, 1)
    x0 = cs.DM.zeros(4, 1)
    lbx = -cs.inf*cs.DM.ones(4, 1)
    ubx = cs.inf*cs.DM.ones(4, 1)
    u = cs.DM.zeros(4, 1)
    c = cs.DM.zeros(19, 1)
    y_obs = cs.DM.zeros(19, 1)
    y_cal = cs.MX.zeros(19, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)

    u[0] = 100
    u[1] = 100
    u[2] = 1
    u[3] = 100
    c[0] = 0.1
    c[1] = 1
    c[2] = 2
    c[3] = 3
    c[4] = 4
    c[5] = 5
    c[6] = 6
    c[7] = 7
    c[8] = 8
    c[9] = 9
    c[10] = 10
    c[11] = 11
    c[12] = 12
    c[13] = 13
    c[14] = 14
    c[15] = 15
    c[16] = 16
    c[17] = 17
    c[18] = 18
    y_obs[0] = 0.00189
    y_obs[1] = 0.1038
    y_obs[2] = 0.268
    y_obs[3] = 0.506
    y_obs[4] = 0.577
    y_obs[5] = 0.604
    y_obs[6] = 0.725
    y_obs[7] = 0.898
    y_obs[8] = 0.947
    y_obs[9] = 0.845
    y_obs[10] = 0.702
    y_obs[11] = 0.528
    y_obs[12] = 0.385
    y_obs[13] = 0.257
    y_obs[14] = 0.159
    y_obs[15] = 0.0869
    y_obs[16] = 0.0453
    y_obs[17] = 0.01509
    y_obs[18] = 0.00189

    x0[0] = 2
    lbx[0] = 0.00001
    ubx[0] = u[0]
    x0[1] = 4
    lbx[1] = 0.00001
    ubx[1] = u[1]
    x0[2] = 0.04
    lbx[2] = 0.00001
    ubx[2] = u[2]
    x0[3] = 2
    lbx[3] = 0.00001
    ubx[3] = u[3]

    b = x[2] + (1-x[2])*x[3]
    d = 7.658
    e = -1
    y_cal[0:19] = (1 + 1/(12*x[1]))**e * x[2]*b**x[1] * (x[1]/6.2832)**(1/2) *\
        (c[0:19]/d)**(x[1] - 1) * cs.exp(x[1] - b*c[0:19]*x[1]/7.658) + \
            (1 + 1/(12*x[0]))**e * (1 - x[2]) * (b/x[3])**x[0] * \
                (x[0]/6.2832)**(1/2) * (c[0:19]/7.658)**(x[0] - 1) * \
                    cs.exp(x[0] - b*c[0:19]*x[0]/(7.658*x[3]))

    lbg[0] = 0
    g[0] = x[2] + (1-x[2])*x[3]
    obj[0:19] = (y_cal[0:19] - y_obs[0:19])**2

    f = cs.sum1(obj)

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
