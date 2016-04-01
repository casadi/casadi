Christen this file BqpdAux.f

c      subroutine gdotx (n, x, ws, lws, v)
c
cc     ==========================================================
cc     Computes G.x where G is Hessian and x is a vector for AMPL
cc     ==========================================================
c
c      implicit none
c
cc     ... declaration of passed parameters
c      integer          n, lws(0:*)
c      double precision x(n), v(n), ws(*)
c
cc     ... declaration of internal variables
c      integer i, j, pjp, row
c
cc     ========================  procedure body  =========================
c
cc     ... form v = W.d from sparse, upper triangular Hessian
c      pjp = lws(0)
c      if (pjp.gt.0) then
c         do i=1,n
c            do j=lws(pjp+i-1),lws(pjp+i)-1
c               row  = lws(j)
c               v(i) = v(i) + ws(j)*x(row)
c               if (row.ne.i) then
c                  v(row) = v(row) + ws(j)*x(i)
c               endif
c            enddo
c         enddo
c      endif
c
c      return
c      end
c
c     ********************************************************************

      subroutine setwsc (mxws0, mxlws0, kk0, ll0)
c     ==========================================================
c     Computes G.x where G is Hessian and x is a vector for AMPL
c     ==========================================================

      implicit none

c     ... declaration of passed parameters
      integer mxws0, mxlws0, kk0, ll0       

c     ... declaration of common blocks
      integer      kk, ll, kkk, lll, mxws, mxlws
      common /wsc/ kk, ll, kkk, lll, mxws, mxlws

c     ... storage map for hessian and scale_mode
      integer         scale_mode, phe
      common /scalec/ scale_mode, phe

c     ========================  procedure body  =========================

      mxws  = mxws0      
      mxlws = mxlws0
      kk    = kk0
      ll    = ll0
      scale_mode = 0
      phe = 1

c    common/wsc/kk,ll,kkk,lll,mxws,mxlws
c    write (*,*) 'kk = ', kk, ' ll = ', ll, ' kkk = ', kkk
c    write (*,*) 'lll = ', lll, ' mxws = ', mxws, ' mxlws = ', mxlws

      return
      end

c     ********************************************************************
      subroutine restorecommon ()
c     ==========================================================
c     Save common variables from the last solve.
c     ==========================================================
      integer char_l, char_lC
      character*10 pname, pnameC

      double precision bqpd_eps, tol, emin, sgnf
      double precision bqpd_epsC, tolC, eminC, sgnfC

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


      common /bqpdcC/      irh1C,naC,na1C,nbC,nb1C,ka1C,kb1C,kc1C,irg1C,
     *lu1C,lvC,lv1C,ll1C
      common /bqpd_countC/ n_bqpd_callsC, n_bqpd_printC
      common /cpnameC/     char_lC, pnameC
      common /densecC/     nsdC,ns1dC,ntdC,nt1dC,nudC,nu1dC,mx1dC,lcdC,
     * lc1dC,lidC,li1dC
      common /epscC/       bqpd_epsC, tolC, eminC
      common /factorcC/    m1C,m2C,mpC,mqC,lastrC,irowC
      common /hesscC/      phlC, phrC, phcC
      common /iprintcC/    iprintC
      common /kkll_maxc_C/ kk_maxC, ll_maxC
      common /minorcC/     cC
      common /mxm1cC/      mxm1C
      common /noutcC/      noutC
      common /QPtolcC/     QPtolC
      common /refactorcC/  nupC,nfreqC
      common /repcC/       sgnfC,nrepC,npivC,nresC
      common /scalecC/     scale_modeC, pheC
      common /sparsecC/    nssC,ns1sC,ntsC,nt1sC,nusC,nu1sC,nxsC,nx1sC,
     * npsC,np1sC,nprofsC,lcC,lc1C,liC,li1C,lmC,lm1C,lpC,lp1C,lqC,lq1C,
     * lrC,lr1C,lsC,ls1C,ltC,lt1C
      common /vstepcC/     vstepC
      common /wscC/        kkC, llC, kkkC, lllC, maxwkC, maxiwkC


      irh1         =  irh1C           
      na           =  naC             
      na1          =  na1C            
      nb           =  nbC             
      nb1          =  nb1C            
      ka1          =  ka1C            
      kb1          =  kb1C            
      kc1          =  kc1C            
      irg1         =  irg1C           
      lu1          =  lu1C            
      lv           =  lvC             
      lv1          =  lv1C            
      ll1          =  ll1C            
      n_bqpd_calls =  n_bqpd_callsC   
      n_bqpd_print =  n_bqpd_printC   
      char_l       =  char_lC         
      pname        =  pnameC          
      nsd          =  nsdC            
      ns1d         =  ns1dC           
      ntd          =  ntdC            
      nt1d         =  nt1dC           
      nud          =  nudC            
      nu1d         =  nu1dC           
      mx1d         =  mx1dC           
      lcd          =  lcdC            
      lc1d         =  lc1dC           
      lid          =  lidC            
      li1d         =  li1dC           
      bqpd_eps     =  bqpd_epsC       
      tol          =  tolC            
      emin         =  eminC           
      m1           =  m1C             
      m2           =  m2C             
      mp           =  mpC             
      mq           =  mqC             
      lastr        =  lastrC          
      irow         =  irowC           
      phl          =  phlC            
      phr          =  phrC            
      phc          =  phcC            
      iprint       =  iprintC         
      kk_max       =  kk_maxC         
      ll_max       =  ll_maxC         
      c            =  cC              
      mxm1         =  mxm1C           
      nout         =  noutC           
      QPtol        =  QPtolC          
      nup          =  nupC            
      nfreq        =  nfreqC          
      sgnf         =  sgnfC           
      nrep         =  nrepC           
      npiv         =  npivC           
      nres         =  nresC           
      scale_mode   =  scale_modeC     
      phe          =  pheC            
      nss          =  nssC            
      ns1s         =  ns1sC           
      nts          =  ntsC            
      nt1s         =  nt1sC           
      nus          =  nusC            
      nu1s         =  nu1sC           
      nxs          =  nxsC            
      nx1s         =  nx1sC           
      nps          =  npsC            
      np1s         =  np1sC           
      nprofs       =  nprofsC         
      lc           =  lcC             
      lc1          =  lc1C            
      li           =  liC             
      li1          =  li1C            
      lm           =  lmC             
      lm1          =  lm1C            
      lp           =  lpC             
      lp1          =  lp1C            
      lq           =  lqC             
      lq1          =  lq1C            
      lr           =  lrC             
      lr1          =  lr1C            
      ls           =  lsC             
      ls1          =  ls1C            
      lt           =  ltC             
      lt1          =  lt1C            
      vstep        =  vstepC          
      kk           =  kkC             
      ll           =  llC             
      kkk          =  kkkC            
      lll          =  lllC            
      maxwk        =  maxwkC          
      maxiwk       =  maxiwkC         

      return
      end

c     ********************************************************************
      subroutine savecommon()
c     ==========================================================
c     Restore common variables from the last checkpoint.
c     ==========================================================

      integer char_l, char_lC
      character*10 pname, pnameC

      double precision bqpd_eps, tol, emin
      double precision bqpd_epsC, tolC, eminC
      double precision sgnf, sgnfC

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


      common /bqpdcC/      irh1C,naC,na1C,nbC,nb1C,ka1C,kb1C,kc1C,irg1C,
     *lu1C,lvC,lv1C,ll1C
      common /bqpd_countC/ n_bqpd_callsC, n_bqpd_printC
      common /cpnameC/     char_lC, pnameC
      common /densecC/     nsdC,ns1dC,ntdC,nt1dC,nudC,nu1dC,mx1dC,lcdC,
     * lc1dC,lidC,li1dC
      common /epscC/       bqpd_epsC, tolC, eminC
      common /factorcC/    m1C,m2C,mpC,mqC,lastrC,irowC
      common /hesscC/      phlC, phrC, phcC
      common /iprintcC/    iprintC
      common /kkll_maxc_C/ kk_maxC, ll_maxC
      common /minorcC/     cC
      common /mxm1cC/      mxm1C
      common /noutcC/      noutC
      common /QPtolcC/     QPtolC
      common /refactorcC/  nupC,nfreqC
      common /repcC/       sgnfC,nrepC,npivC,nresC
      common /scalecC/     scale_modeC, pheC
      common /sparsecC/    nssC,ns1sC,ntsC,nt1sC,nusC,nu1sC,nxsC,nx1sC,
     * npsC,np1sC,nprofsC,lcC,lc1C,liC,li1C,lmC,lm1C,lpC,lp1C,lqC,lq1C,
     * lrC,lr1C,lsC,ls1C,ltC,lt1C
      common /vstepcC/     vstepC
      common /wscC/        kkC, llC, kkkC, lllC, maxwkC, maxiwkC


      irh1C         =   irh1
      naC           =   na
      na1C          =   na1
      nbC           =   nb
      nb1C          =   nb1
      ka1C          =   ka1
      kb1C          =   kb1
      kc1C          =   kc1
      irg1C         =   irg1
      lu1C          =   lu1
      lvC           =   lv
      lv1C          =   lv1
      ll1C          =   ll1
      n_bqpd_callsC =   n_bqpd_calls
      n_bqpd_printC =   n_bqpd_print
      char_lC       =   char_l
      pnameC        =   pname
      nsdC          =   nsd
      ns1dC         =   ns1d
      ntdC          =   ntd
      nt1dC         =   nt1d
      nudC          =   nud
      nu1dC         =   nu1d
      mx1dC         =   mx1d
      lcdC          =   lcd
      lc1dC         =   lc1d
      lidC          =   lid
      li1dC         =   li1d
      bqpd_epsC     =   bqpd_eps
      tolC          =   tol
      eminC         =   emin
      m1C           =   m1
      m2C           =   m2
      mpC           =   mp
      mqC           =   mq
      lastrC        =   lastr
      irowC         =   irow
      phlC          =   phl
      phrC          =   phr
      phcC          =   phc
      iprintC       =   iprint
      kk_maxC       =   kk_max
      ll_maxC       =   ll_max
      cC            =   c
      mxm1C         =   mxm1
      noutC         =   nout
      QPtolC        =   QPtol
      nupC          =   nup
      nfreqC        =   nfreq
      sgnfC         =   sgnf
      nrepC         =   nrep
      npivC         =   npiv
      nresC         =   nres
      scale_modeC   =   scale_mode
      pheC          =   phe
      nssC          =   nss
      ns1sC         =   ns1s
      ntsC          =   nts
      nt1sC         =   nt1s
      nusC          =   nus
      nu1sC         =   nu1s
      nxsC          =   nxs
      nx1sC         =   nx1s
      npsC          =   nps
      np1sC         =   np1s
      nprofsC       =   nprofs
      lcC           =   lc
      lc1C          =   lc1
      liC           =   li
      li1C          =   li1
      lmC           =   lm
      lm1C          =   lm1
      lpC           =   lp
      lp1C          =   lp1
      lqC           =   lq
      lq1C          =   lq1
      lrC           =   lr
      lr1C          =   lr1
      lsC           =   ls
      ls1C          =   ls1
      ltC           =   lt
      lt1C          =   lt1
      vstepC        =   vstep
      kkC           =   kk
      llC           =   ll
      kkkC          =   kkk
      lllC          =   lll
      maxwkC        =   maxwk
      maxiwkC       =   maxiwk

      return
      end
c     **********************************  e n d  ***************************

