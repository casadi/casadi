c     ********************************************************************
      subroutine writewsc ()
      implicit double precision (a-h,o-z), integer (i-n)
      common /bqpdc/      irh1,na,na1,nb,nb1,ka1,kb1,kc1,irg1,lu1,lv,
     * lv1,ll1
      common /bqpd_count/ n_bqpd_calls, n_bqpd_print
      common /cpname/     char_l, pname
      common /densec/     nsd,ns1d,ntd,nt1d,nud,nu1d,mx1d,lcd,lc1d,lid,
     * li1d
      common /epsc/       bqpd_eps, tol, emin
      common /factorc/    m1,m2,mp,mq,lastr,irow
      common /hessc/      phl, phr, phc
      common /iprintc/    iprint
      common /kkll_maxc_/ kk_max, ll_max
      common /minorc/     c
      common /mxm1c/      mxm1
      common /noutc/      nout
      common /QPtolc/     QPtol
      common /refactorc/  nup,nfreq
      common /repc/       sgnf,nrep,npiv,nres
      common /scalec/     scale_mode, phe
      common /sparsec/    nss,ns1s,nts,nt1s,nus,nu1s,nxs,nx1s,nps,np1s,
     * nprofs,lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls,ls1,lt,lt1
      common /vstepc/     vstep
      common /wsc/        kk, ll, kkk, lll, maxwk, maxiwk


      print *, 'writing wsc'
      print *, irh1,na,na1,nb,nb1,ka1,kb1,kc1,irg1,lu1,lv,lv1,ll1
      print *, n_bqpd_calls, n_bqpd_print
      print *, char_l, pname
      print *, nsd,ns1d,ntd,nt1d,nud,nu1d,mx1d,lcd,lc1d,lid,li1d
      print *, 'writing epsc'
      print *, bqpd_eps, tol, emin
      print *, m1,m2,mp,mq,lastr,irow
      print *, phl, phr, phc
      print *, iprint
      print *, kk_max, ll_max
      print *, c
      print *, mxm1
      print *, nout
      print *, QPtol
      print *, nup,nfreq
      print *, sgnf,nrep,npiv,nres
      print *, scale_mode, phe
      print *, nss,ns1s,nts,nt1s,nus,nu1s,nxs,nx1s,nps,np1s,nprofs,lc
      print *, lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls,ls1,lt,lt1
      print *, vstep
      print *, kk, ll, kkk, lll, maxwk, maxiwk
      print *, 'stop writing wsc'

      return
      end



      subroutine defaultcommon()
      implicit double precision (a-h,o-z), integer (i-n)
      common/bqpdc/irh1,na,na1,nb,nb1,ka1,kb1,kc1,irg1,lu1,lv,lv1,ll1
      common/epsc/eps,tol,emin
      common/vstepc/vstep
      common/repc/sgnf,nrep,npiv,nres
      common/refactorc/nup,nfreq
      common/wsc/kk,ll,kkk,lll,mxws,mxlws
      common/alphac/alpha,rp,pj,qqj,qqj1
      common /sparsec/    nss,ns1s,nts,nt1s,nus,nu1s,nxs,nx1s,nps,np1s,
     * nprofs,lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls,ls1,lt,lt1
      common /factorc/    m1,m2,mp,mq,lastr,irow
      common /hessc/      phl, phr, phc
      common /minorc/     c
      common /scalec/     scale_mode, phe
      irh1  = 0
      na    = 0
      na1   = 0
      nb    = 0
      nb1   = 0
      ka1   = 0
      kb1   = 0
      kc1   = 0
      irg1  = 0
      lu1   = 0
      lv    = 0
      lv1   = 0
      ll1   = 0
      eps   = 1111.D-19
      tol   = 1.D-12 
      emin  = 0.D0
      vstep = 0.D0
      sgnf  = 1.D-4
      nrep  = 2
      npiv  = 3
      nres  = 2
      nup   = 0
      nfreq = 500 

      kk    = 0
      ll    = 0
      kkk   = 0
      lll   = 0
      mxws  = 0
      mxlws = 0

      alpha = 0
      rp    = 0
      pj    = 0
      qqj   = 0
      qqj1  = 0

      nss   = 0

      ns1s  = 0
      nts   = 0
      nt1s  = 0
      nus   = 0
      nu1s  = 0
      nxs   = 0
      nx1s  = 0
      nps   = 0
      np1s  = 0
      nprofs= 0
      lc    = 0
      lc1   = 0
      li    = 0
      li1   = 0
      lm    = 0
      lm1   = 0
      lp    = 0
      lp1   = 0
      lq    = 0
      lq1   = 0
      lr    = 0
      lr1   = 0
      ls    = 0
      ls1   = 0
      lt    = 0
      lt1   = 0

      m1    = 0
      m2    = 0
      mp    = 0
      mq    = 0
      lastr = 0
      irow  = 0

      phl   = 0.D0
      phr   = 0.D0
      phc   = 0.D0
      
      c     = 0.D0

      scale_mode = 0.D0
      phe        = 0.D0
c    *   /6.D-8, 1.D-6, 0.D0, 1.D-1,  2,    3,    3,   100,   0,  0/
c     double precision eps, tol, emin, sgnf
c     integer nrep, npriv, nres, nup, nfreq, kk, ll
c  single length tolerances
c     data  eps,   tol, emin, sgnf, nrep, npiv, nres, nfreq, kk, ll
c    *   /6.D-8, 1.D-6, 0.D0, 1.D-1,  2,    3,    3,   100,   0,  0/
c  double length tolerances
c     emin =  8.D-12
c     data  eps,    tol,   emin, sgnf, nrep, npiv, nres, nfreq, kk, ll 
c    * /1111.D-19, 1.D-12, 0.D0, 1.D-4,  2,    3,   2,   500,   0,  0/
      end
