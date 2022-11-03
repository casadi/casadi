# This is an automatically generated file converting from the apm format
import casadi as cs


def hs088():
    # The optimal objective is (if given in):
    f_opt = 1.36265681
    x_opt = cs.DM([1.07431, 0.456611])
    x = cs.MX.sym('x', 2)
    x0 = cs.DM.zeros(2, 1)
    lbx = -cs.inf*cs.DM.ones(2, 1)
    ubx = cs.inf*cs.DM.ones(2, 1)

    mu = cs.DM.zeros(30, 1)
    mu12 = cs.DM.zeros(30, 1)
    tg12 = cs.DM.zeros(30, 1)
    mu24 = cs.DM.zeros(30, 1)
    tg24 = cs.DM.zeros(30, 1)
    mu48 = cs.DM.zeros(30, 1)
    snmu = cs.DM.zeros(30, 1)
    csmu = cs.DM.zeros(30, 1)
    snmuxcsmu = cs.DM.zeros(30, 1)
    aux2 = cs.DM.zeros(30, 1)
    aux4 = cs.DM.zeros(30, 1)
    A = cs.DM.zeros(30, 1)

    emx_1 = cs.MX.zeros(30, 1)
    emx_2 = cs.MX.zeros(30, 1)
    rhoaux_2 = cs.MX.zeros(30, 1)
    aux1 = cs.MX.zeros(30, 1)
    aux3 = cs.MX.zeros(30, 1)
    hsum = cs.MX.zeros(30, 1)

    rho = cs.MX.zeros(30, 1)
    g = cs.MX.zeros(1, 1)
    lbg = -cs.inf*cs.DM.ones(1, 1)
    ubg = cs.inf*cs.DM.ones(1, 1)

    n = 2
    mu12[0] = 0.86033358901938
    mu12[1] = 3.42561845948172
    mu12[2] = 6.43729817917194
    mu12[3] = 9.52933440536196
    mu12[4] = 12.6452872238566
    mu12[5] = 15.7712848748158
    mu12[6] = 18.90240995686
    mu12[7] = 22.0364967279385
    mu12[8] = 25.1724463266466
    mu12[9] = 28.309642854452
    mu12[10] = 31.4477146375462
    mu12[11] = 34.5864242152889
    mu12[12] = 37.7256128277765
    mu12[13] = 40.865170330488
    mu12[14] = 44.0050179208308
    mu12[15] = 47.145097736761
    mu12[16] = 50.2853663377736
    mu12[17] = 53.4257904773946
    mu12[18] = 56.5663442798215
    mu12[19] = 59.7070073053354
    mu12[20] = 62.8477631944544
    mu12[21] = 65.9885986984903
    mu12[22] = 69.1295029738952
    mu12[23] = 72.2704670603089
    mu12[24] = 75.4114834888481
    mu12[25] = 78.5525459842429
    mu12[26] = 81.6936492356016
    mu12[27] = 84.8347887180422
    mu12[28] = 87.9759605524932
    mu12[29] = 91.1171613944647

    tg12[0:30] = cs.tan(mu12[0:30])
    mu24[0:30] = (1 + mu12[0:30]**2*(1 + tg12[0:30]**2)) / (tg12[0:30] + mu12[0:30]*(1 + tg12[0:30]**2))
    tg24[0:30] = cs.tan(mu24[0:30])
    mu48[0:30] = (1 + mu24[0:30]**2*(1 + tg24[0:30]**2)) / (tg24[0:30] + mu24[0:30]*(1 + tg24[0:30]**2))
    mu[0:30] = mu48[0:30]
    isign = cs.DM.ones(30)
    isign[1::2] = -1
    snmu[0:30] = isign[0:30] * cs.sqrt(1/(1 + mu[0:30]**2))
    csmu[0:30] = isign[0:30] * cs.sqrt(mu[0:30]**2/(1 + mu[0:30]**2))
    snmuxcsmu[0:30] = mu[0:30]/(1 + mu[0:30]**2)
    aux2[0:30] = snmuxcsmu[0:30]/(2*mu[0:30]) + 1/2
    aux4[0:30] = (-2)*snmu[0:30]/mu[0:30] + 2*csmu[0:30]
    A[0:30] = 2*snmu[0:30]/(mu[0:30] + snmuxcsmu[0:30])

    # Intermediates
    emx_1[0:30] = cs.exp(-mu[0:30]**2*x[0]**2)
    emx_2[0:30] = cs.exp(-mu[0:30]**2*x[1]**2)
    rhoaux_2[0:30] = 1*emx_1[0:30]*emx_2[0:30] - 2*emx_2[0:30] + 1
    rho[0:30] = (-1)*rhoaux_2[0:30]/mu[0:30]**2
    aux1[0:30] = A[0:30]**2*rhoaux_2[0:30]**2
    aux3[0:30] = A[0:30]*rho[0:30]
    hsum[0:30] = aux1[0:30]*aux2[0:30] + aux3[0:30]*aux4[0:30]

    x0[0] = 0.5
    x0[1] = -0.5

    ubg[0] = 0.0001
    g[0] = 2/15 + cs.sum1(hsum)
    f = x[0]**2 + x[1]**2

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
