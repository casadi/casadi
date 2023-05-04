# This is an automatically generated file converting from the apm format
import casadi as cs


def hs110():
    # The optimal objective is (if given in):
    f_opt = -45.77846971
    x_opt = cs.DM([9.35027, 9.35027, 9.35027, 9.35027, 9.35027, 9.35027, 9.35027, 9.35027, 9.35027, 9.35027])
    x = cs.MX.sym('x', 10)
    obj = cs.MX.zeros(10, 1)
    x0 = cs.DM.zeros(10, 1)
    lbx = -cs.inf*cs.DM.ones(10, 1)
    ubx = cs.inf*cs.DM.ones(10, 1)
    g = []
    lbg = []
    ubg = []

    x0[0:10] = 9
    lbx[0:10] = 2.001
    ubx[0:10] = 9.999

    prod = 1
    for i in range(10):
        prod = prod*x[i]**0.2

    obj[0:10] = cs.log(x[0:10]-2)**2 + cs.log(10-x[0:10])**2

    f = cs.sum1(obj) - prod

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
