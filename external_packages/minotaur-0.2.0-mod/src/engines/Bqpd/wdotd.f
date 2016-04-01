Christen this file wdotd.f

      subroutine wdotd (n, d, ws, lws, v)

c     ==========================================================
c     Computes W.d where W is Hessian and d is a vector for AMPL
c     Assumes v=0 on entry (OK, if called from gdotx, see QPsolve*.f)
c     ==========================================================

      implicit none

c     ... declaration of passed parameters
      integer          n, lws(0:*)
      double precision d(n), v(n), ws(*), vv, dd

c     ... declaration of internal variables
      integer i, j, pjp, row

c     ... storage map for hessian 
      integer        phl, phr, phc
      common /hessc/ phl, phr, phc

c     ========================  procedure body  =========================

c     ... form v = W.d from sparse, upper triangular Hessian
      phl = 1
      if (phl.gt.0) then
         pjp = lws(0)
c        write (*,*) 'pjp = ', pjp
         do i=1,n
            vv = v(i)
            dd = d(i)
            do j=lws(pjp+i-1),lws(pjp+i)-1
               row  = lws(j)
               vv = vv + ws(j)*d(row)
               if (row.ne.i) then
                  v(row) = v(row) + ws(j)*dd
               endif
            enddo
            v(i) = vv
         enddo
      endif

c      write(8,*) '  d = ',(d(i),i=1,n)
c      write(8,*) 'W.d = ',(v(i),i=1,n)

      return
      end

c     *******************************************************************

      subroutine ident_Hessian (n,lws,ws,a)

c     ========================================================================
c     Set initial Hessian = I and gradient = 0, to get l_2 closest feas. point 
c     ========================================================================

      implicit none

c     ... declaration of passed parameters
      integer          n, lws(0:*)
      double precision ws(*), a(*)

c     ... declaration of internal variables
      integer i, pjp

c     ... storage map for hessian 
      integer        phl, phr, phc
      common /hessc/ phl, phr, phc

c     ========================  procedure body  =========================

      phl    = 1
      pjp    = n+1
      lws(0) = pjp
      do i=1,n
         lws(pjp+i-1) = i
         lws(i)       = i
         ws(i)        = 1.D0
         a(i)         = 0.D0
      enddo
      lws(pjp+n) = n+1

      return
      end
