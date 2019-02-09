// C-REPLACE "fabs" "casadi_fabs"
// C-REPLACE "sign" "casadi_sign"

// NOLINT(legal/copyright)
// SYMBOL "lsqr_sym_ortho"
template<typename T1>
void casadi_lsqr_sym_ortho(T1 a, T1 b, T1* cs, T1* sn, T1* rho) {
    T1 tau;
    if (b == 0) {
        *cs = sign(a);
        *sn = 0;
        *rho = fabs(a);
    } else if (a==0) {
        *cs = 0;
        *sn = sign(b);
        *rho = fabs(b);
    } else if (fabs(b)>fabs(a)) {
        tau = a/b;
        *sn = sign(b)/sqrt(1+tau*tau);
        *cs = (*sn)*tau;
        *rho = b/(*sn);
    } else {
        tau = b/a;
        *cs = sign(a)/sqrt(1+tau*tau);
        *sn = (*cs)*tau;
        *rho = a/(*cs);
    }
}

// NOLINT(legal/copyright)
// SYMBOL "lsqr_single_solve"
// Ref: scipy
template<typename T1>
int casadi_lsqr_single_solve(const T1* A, T1* x, casadi_int tr, const casadi_int* sp, T1* w) {
    casadi_int m, n, i;
    T1 damp, atol, btol, conlim, ctol, anorm, acond, dampsq, ddnorm, res2, xnorm, xxnorm, z;
    T1 cs2, sn2, alpha, beta, rhobar, phibar, bnorm, rnorm, arnorm, rhobar1, cs1, sn1, psi;
    T1 cs, sn, rho, theta, phi, tau, t1, t2, n2dk, delta, gambar, rhs, zbar, gamma, res1;
    T1 r1sq, r1norm, test1, test2, test3, rtol;
    casadi_int iter_lim, itn, istop;
    T1 *u, *v, *xx, *ww, *dk;

    m = sp[0];
    n = sp[1];

    damp = 0;
    atol = 1e-15;
    btol = 1e-15;
    conlim = 1e8;
    iter_lim = 10000;

    itn = 0;
    istop = 0;

    ctol = 0;
    if (conlim > 0) ctol = 1/conlim;
    anorm = 0;
    acond = 0;
    dampsq = damp*damp;
    ddnorm = 0;
    res2 = 0;
    xnorm = 0;
    xxnorm = 0;
    z = 0;
    cs2 = -1;
    sn2 = 0;

    u = w;  w+= m; casadi_copy(x, m, u);
    v = w;  w+= n; casadi_clear(v, n);
    xx = w; w+= n; casadi_clear(xx, n);
    ww = w; w+= n; casadi_clear(v, n);
    dk = w; w+= n;

    alpha = 0;
    beta = casadi_norm_2(m, u);

    if (beta>0) {
      for (i=0;i<m;++i) u[i]*=1/beta;
      casadi_mv(A, sp, u, v, !tr);
      alpha = casadi_norm_2(n, v);
    }

    if (alpha>0) {
      for (i=0;i<n;++i) v[i]*=1/alpha;
      casadi_copy(v, n, ww);
    }

    rhobar = alpha;
    phibar = beta;
    bnorm = beta;
    rnorm = beta;
    arnorm = alpha * beta;

    while (itn<iter_lim) {
      itn++;
      for (i=0;i<m;++i) u[i]*=-alpha;
      casadi_mv(A, sp, v, u, tr);
      beta = casadi_norm_2(m, u);

      if (beta>0) {
        for (i=0;i<m;++i) u[i]*=1/beta;
        anorm = sqrt(anorm*anorm + alpha*alpha+beta*beta+damp*damp);
        for (i=0;i<n;++i) v[i]*=-beta;
        casadi_mv(A, sp, u, v, !tr);
        alpha = casadi_norm_2(n, v);
        if (alpha>0) for (i=0;i<n;++i) v[i]*=1/alpha;
      }

      rhobar1 = sqrt(rhobar*rhobar+damp*damp);

      cs1 = rhobar / rhobar1;
      sn1 = damp / rhobar1;
      psi = sn1 * phibar;
      phibar *= cs1;

      casadi_lsqr_sym_ortho(rhobar1, beta, &cs, &sn, &rho);

      theta = sn * alpha;
      rhobar = -cs * alpha;
      phi = cs * phibar;
      phibar *= sn;
      tau = sn * phi;

      t1 = phi / rho;
      t2 = -theta / rho;

      for (i=0;i<n;++i) dk[i]=ww[i]/rho;

      for (i=0; i<n; ++i) xx[i] += t1*ww[i];
      for (i=0; i<n; ++i) ww[i] = v[i] + t2*ww[i];

      n2dk = casadi_norm_2(n, dk);
      ddnorm += n2dk*n2dk;

      delta = sn2 * rho;
      gambar = -cs2 * rho;
      rhs = phi - delta * z;
      zbar = rhs / gambar;
      xnorm = sqrt(xxnorm + zbar*zbar);
      gamma = sqrt(gambar*gambar + theta*theta);
      cs2 = gambar / gamma;
      sn2 = theta / gamma;
      z = rhs / gamma;
      xxnorm += z*z;

      acond = anorm * sqrt(ddnorm);
      res1 = phibar*phibar;
      res2 += psi*psi;
      rnorm = sqrt(res1+res2);
      arnorm = alpha*fabs(tau);

      r1sq = rnorm*rnorm - dampsq * xxnorm;
      r1norm = sqrt(fabs(r1sq));
      if (r1sq < 0) r1norm = -r1norm;

      test1 = rnorm / bnorm;
      test2 = arnorm / (anorm * rnorm);
      test3 = 1 / acond;
      t1 = test1 / (1 + anorm * xnorm / bnorm);
      rtol = btol + atol * anorm * xnorm / bnorm;

      if (itn >= iter_lim) istop = 7;
      if (1 + test3 <= 1) istop = 6;
      if (1 + test2 <= 1) istop = 5;
      if (1 + t1 <= 1) istop = 4;

      if (test3 <= ctol) istop = 3;
      if (test2 <= atol) istop = 2;
      if (test1 <= rtol) istop = 1;

      if (istop != 0) break;

    }
    casadi_copy(xx, m, x);
    return 0;
}

// NOLINT(legal/copyright)
// SYMBOL "lsqr_solve"
template<typename T1>
int casadi_lsqr_solve(const T1* A, T1* x, casadi_int nrhs, casadi_int tr,
        const casadi_int* sp, T1* w) {
    casadi_int i;
    for (i=0; i<nrhs;++i) {
      if (casadi_lsqr_single_solve(A, x+i*sp[1], tr, sp, w)) return 1;
    }
    return 0;
}
