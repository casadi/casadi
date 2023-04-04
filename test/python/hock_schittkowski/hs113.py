# This is an automatically generated file converting from the apm format
import casadi as cs


def hs113():
    # The optimal objective is (if given in):
    f_opt = 24.3062091
    x_opt = cs.DM([2.172, 2.36368, 8.77393, 5.09598, 0.990655, 1.43057, 1.32164, 9.82873, 8.28009, 8.37593])
    x = cs.MX.sym('x', 10)
    x0 = cs.DM.zeros(10, 1)
    lbx = -cs.inf*cs.DM.ones(10, 1)
    ubx = cs.inf*cs.DM.ones(10, 1)
    g = cs.MX.zeros(8, 1)
    lbg = -cs.inf*cs.DM.ones(8, 1)
    ubg = cs.inf*cs.DM.ones(8, 1)

    x0[0] = 2
    x0[1] = 3
    x0[2] = 5
    x0[3] = 5
    x0[4] = 1
    x0[5] = 2
    x0[6] = 7
    x0[7] = 3
    x0[8] = 6
    x0[9] = 10

    lbg[0] = 0
    g[0] = 105-4*x[0]-5*x[1]+3*x[6]-9*x[7]
    lbg[1] = 0
    g[1] = -10*x[0]+8*x[1]+17*x[6]-2*x[7]
    lbg[2] = 0
    g[2] = 8*x[0]-2*x[1]-5*x[8]+2*x[9]+12
    lbg[3] = 0
    g[3] = -3*(x[0]-2)**2-4*(x[1]-3)**2-2*x[2]**2+7*x[3]+120
    lbg[4] = 0
    g[4] = -5*x[0]**2-8*x[1]-(x[2]-6)**2+2*x[3]+40
    lbg[5] = 0
    g[5] = -.5*(x[0]-8)**2-2*(x[1]-4)**2-3*x[4]**2+x[5]+30
    lbg[6] = 0
    g[6] = -x[0]**2-2*(x[1]-2)**2+2*x[0]*x[1]-14*x[4]+6*x[5]
    lbg[7] = 0
    g[7] = 3*x[0]-6*x[1]-12*(x[8]-8)**2+7*x[9]
    f = x[0]**2+x[1]**2+x[0]*x[1]-14*x[0]-16*x[1]+(x[2]-10)**2+4*(x[3]-5)**2 +(x[4]-3)**2+2*(x[5]-1)**2+5*x[6]**2+7*(x[7]-11)**2+2*(x[8]-10)**2 +(x[9]-7)**2+45

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
