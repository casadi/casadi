# This is an automatically generated file converting from the apm format

import casadi as cs

def hs054():
    # The optimal objective is (if given in):
    f_opt = -cs.exp(-27/280)
    x_opt = cs.DM([91600/7, 79/70, 2e6, 10, 1e-3, 1e8])
    x = cs.MX.sym('x', 6)
    lbx = cs.DM([0, -10, 0, 0, -1, 0])
    ubx = cs.DM([2e4, 10, 1e7, 20, 1, 2e8])
    y_1 = (x[0] - 10000)/8000
    y_2 = (x[1] - 1)/1
    y_3 = (x[2] - 2000000)/7000000
    y_4 = (x[3] - 10)/50
    y_5 = (x[4] - 1/1000)*20
    y_6 = (x[5] - 100000000)/500000000
    h1 = (y_1**2 + y_1*y_2*2/5 + y_2**2)*25/24
    h2 = y_3**2 + y_4**2 + y_5**2 + y_6**2
    f = -cs.exp(-(h1 + h2)/2)

    g = x[0] + 4e3*x[1] - 1.76e4
    lbg = cs.DM([0])
    ubg = cs.DM([0])

    x0 = cs.DM([6e3, 1.5, 4e6, 2, 3e-3, 5e7])
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
