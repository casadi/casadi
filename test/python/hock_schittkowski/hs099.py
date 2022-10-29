# This is an automatically generated file converting from the apm format
import casadi as cs


def hs099():
    # The optimal objective is (if given in):
    f_opt = -0.831079892e+9
    x_opt = cs.DM([0.542468, 0.529021, 0.508449, 0.480269, 0.451236, 0.409183, 0.352788])
    x = cs.MX.sym('x', 7)
    x0 = cs.DM.zeros(7, 1)
    lbx = -cs.inf*cs.DM.ones(7, 1)
    ubx = cs.inf*cs.DM.ones(7, 1)
    a = cs.DM.zeros(8, 1)
    t = cs.DM.zeros(8, 1)
    r_tmp = cs.MX.zeros(8, 1)
    q_tmp = cs.MX.zeros(8, 1)
    s_tmp = cs.MX.zeros(8, 1)
    g = cs.MX.zeros(2, 1)
    lbg = -cs.inf*cs.DM.ones(2, 1)
    ubg = cs.inf*cs.DM.ones(2, 1)

    a[0] = 0
    a[1] = 50
    a[2] = 50
    a[3] = 75
    a[4] = 75
    a[5] = 75
    a[6] = 100
    a[7] = 100
    b = 32
    t[0] = 0
    t[1] = 25
    t[2] = 50
    t[3] = 100
    t[4] = 150
    t[5] = 200
    t[6] = 290
    t[7] = 380

    x0[0:7] = 0.5
    lbx[0:7] = 0
    ubx[0:7] = 1.58

    r_tmp[0] = 0
    for i in range(1, 8):
        r_tmp[i] = a[i]*(t[i]-t[i-1]) * cs.cos(x[i-1]) + r_tmp[i-1]

    s_tmp[0] = 0
    for i in range(1, 8):
        s_tmp[i] = (t[i]-t[i-1])*(a[i]*cs.sin(x[i-1]) - b) + s_tmp[i-1]
    q_tmp[0] = 0
    for i in range(1, 8):
        q_tmp[i] = 0.5*(t[i]-t[i-1])**2*(a[i]*cs.sin(x[i-1]) - b) + (t[i]-t[i-1])*s_tmp[i-1] + q_tmp[i-1]

    lbg[0] = 1.0e+5
    ubg[0] = 1.0e+5
    g[0] = q_tmp[7]
    lbg[1] = 1.0e+3
    ubg[1] = 1.0e+3
    g[1] = s_tmp[7]

    f = -r_tmp[7]**2

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
