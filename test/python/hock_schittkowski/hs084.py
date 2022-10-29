# This is an automatically generated file converting from the apm format

import casadi as cs

def hs084():
    # The optimal objective is (if given in):
    f_opt = -5280335.133
    x_opt = cs.DM([4.53743, 2.4, 60, 9.3, 7])
    x = cs.MX.sym('x', 5)
    obj = cs.MX.zeros(1, 1)
    x0 = cs.DM.zeros(5, 1)
    lbx = -cs.inf*cs.DM.ones(5, 1)
    ubx = cs.inf*cs.DM.ones(5, 1)
    a = cs.DM.zeros(21, 1)
    l = cs.DM.zeros(5, 1)
    u = cs.DM.zeros(5, 1)
    xi = cs.DM.zeros(5, 1)
    g = cs.MX.zeros(3, 1)
    lbg = -cs.inf*cs.DM.ones(3, 1)
    ubg = cs.inf*cs.DM.ones(3, 1)
    
    a[0] =    -24345
    a[1] =  -8720288.849
    a[2] =    150512.5253
    a[3] =      -156.6950325
    a[4] =    476470.3222
    a[5] =    729482.8271
    a[6] =   -145421.402
    a[7] =      2931.1506
    a[8] =       -40.427932
    a[9] =     5106.192
    a[10] =    15711.36
    a[11] =  -155011.1084
    a[12] =     4360.53352
    a[13] =       12.9492344
    a[14] =    10236.884
    a[15] =    13176.786
    a[16] =  -326669.5104
    a[17] =     7390.68412
    a[18] =      -27.8986976
    a[19] =    16643.076
    a[20] =    30988.146
    l[0] = 0
    l[1] = 1.2
    l[2] = 20
    l[3] = 9
    l[4] = 6.5
    u[0] = 1000
    u[1] = 2.4
    u[2] = 60
    u[3] = 9.3
    u[4] = 7
    xi[0] = 2.52
    xi[1] = 2
    xi[2] = 37.5
    xi[3] = 9.25
    xi[4] = 6.8
    
    x0[0:5]= xi[0:5]
    lbx[0:5]=l[0:5]
    ubx[0:5]=u[0:5]
    
    g[0] =  a[6]*x[0] + a[7]*x[0]*x[1] + a[8]*x[0]*x[2] + a[9]*x[0]*x[3] + a[10]*x[0]*x[4]
    lbg[0] = 0
    ubg[0] =  294000
    g[1] =  a[11]*x[0] + a[12]*x[0]*x[1] + a[13]*x[0]*x[2] + a[14]*x[0]*x[3] + a[15]*x[0]*x[4]
    lbg[1] = 0
    ubg[1] =  294000
    g[2] =  a[16]*x[0] + a[17]*x[0]*x[1] + a[18]*x[0]*x[2] + a[19]*x[0]*x[3] + a[20]*x[0]*x[4]
    lbg[2] = 0
    ubg[2] =  277200
    obj = -a[0] - a[1]*x[0] - a[2]*x[0]*x[1] - a[3]*x[0]*x[2] - a[4]*x[0]*x[3] - a[5]*x[0]*x[4]
    
    
    f = obj
    
    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
