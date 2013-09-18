module SNOPTModule 
 
end module

subroutine wsnopt () bind ( C, name="snopt" )
      implicit none
      integer     maxF, maxn, nxname, nFname, lenA, lenG, lencw, leniw, lenrw
      parameter   ( maxF   = 1000, maxn  = 1500, lenA   = 5000, lenG  = 5000, nxname = 1, nFname = 1 )
      integer  iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG), xstate(maxn), Fstate(maxF)
      character*8    Prob, xnames(nxname), Fnames(nFname)
      double precision   ObjAdd, sInf,  A(lenA), Flow(maxF), Fupp(maxF), F(maxF), &
      xlow(maxn), xupp(maxn), x(maxn), Fmul(maxF), xmul(maxn)

!     SNOPT workspace

      parameter (lenrw = 1000000,  leniw =  500000, lencw =     500)
      double precision rw(lenrw)
      integer iw(leniw)
      integer DerOpt, Errors, neA, neG, ObjRow, INFO, iPrint, &
       iPrt, iSpecs, iSumm, iSum, itnlim, mincw, miniw, minrw, &
       nF, n, nInf, nOut, nS
      logical   byname
      character   lfile*20, cw(lencw)*8
      external    userf, userfg
!     ------------------------------------------------------------------
      integer            Cold,       Basis,      Warm
      parameter         (Cold   = 0, Basis  = 1, Warm  = 2)
!     ------------------------------------------------------------------
!     Specify some of the SNOPT files.
!     iSpecs  is the Specs file   (0 if none).
!     iPrint  is the Print file   (0 if none).
!     iSumm   is the Summary file (0 if none).
!     nOut    is an output file used here by sntest.

      iSpecs =  4  ! equivalenced to banana.spc
      iPrint =  9  ! equivalenced to banana.out
      iSumm  =  6  ! summary file goes to standard output...
      nOut   =  6  ! ... as do messages from this program.

      byname = .true.

      if ( byname ) then

!       Unix and DOS systems.  Open the Specs and print files.

         lfile = 't2bananaa.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 't2bananaa.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

!     ------------------------------------------------------------------
!     Set options to default values.
!     ------------------------------------------------------------------
      call snInit   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

!     Set up the problem to be solved.
!     No derivatives are provided for the first run.

      Errors = 0

      call t2dat0 ( Errors, maxF, maxn, Prob, nF, n, ObjAdd, &
       ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul, cw, &
       lencw, iw, leniw, rw, lenrw )
      if (Errors .ne. 0) go to 910

!     snOptA will compute the Jacobian by finite-differences.
!     Call  snJac  to define the coordinate arrays (iAfun,jAvar,A)
!     and (iGfun, jGvar).

      call snJac  ( INFO, nF, n, userf, iAfun, jAvar, lenA, neA, A, iGfun, &
       jGvar, lenG, neG, x, xlow, xupp, mincw, miniw, minrw, cw, lencw, &
       iw, leniw, rw, lenrw, cw, lencw, iw, leniw, rw, lenrw )
      if (INFO .ne. 102) go to 900

!     ------------------------------------------------------------------
!     Warn snOptA that userf does not provide the derivatives.
!     The parameters iPrt and iSum may refer to the Print and Summary
!     file respectively.  Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      DerOpt = 0
      iPrt   = 0
      iSum   = 0
      call snSeti ( 'Derivative option', DerOpt, iPrt, iSum, Errors, cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Go for it, using a Cold start (Start = 0).
!     ------------------------------------------------------------------
      call snOptA ( Cold, nF, n, nxname, nFname, &
          ObjAdd, ObjRow, Prob, userf, &
          iAfun, jAvar, lenA, neA, A, &
          iGfun, jGvar, lenG, neG, &
          xlow, xupp, xnames, Flow, Fupp, Fnames, &
          x, xstate, xmul, F, Fstate, Fmul, &
          INFO, mincw, miniw, minrw, &
          nS, nInf, sInf, &
          cw, lencw, iw, leniw, rw, lenrw, &
          cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'snOptA finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0) write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 30) go to 900

!     ------------------------------------------------------------------
!     Read a Specs file.
!     ------------------------------------------------------------------
      call snSpec ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 910
      end if

!     ------------------------------------------------------------------
!     Generate the same problem, but with derivatives provided.
!     ------------------------------------------------------------------
      Errors = 0
      call t2dat1( Errors, Prob, maxF, maxn, nF, n, &
          iAfun, jAvar, lenA, neA, A, &
          iGfun, jGvar, lenG, neG, &
          ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul, &
          cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .ne. 0) then
         go to 910
      end if

!     ------------------------------------------------------------------
!     Specify options that were not set in the Specs file.
!     ------------------------------------------------------------------
      DerOpt = 1
      call snSeti ( 'Derivative option',  DerOpt, iPrt, iSum, Errors, cw, lencw, iw, leniw, rw, lenrw )

      itnlim = 20000
      call snSeti( 'Iterations        ', itnlim, iPrt, iSum, Errors, cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Solve the problem again, using a Cold start (Start = 0).
!     ------------------------------------------------------------------
      call snOptA ( Cold, nF, n, nxname, nFname, &
          ObjAdd, ObjRow, Prob, userfg, &
          iAfun, jAvar, lenA, neA, A, &
          iGfun, jGvar, lenG, neG, &
          xlow, xupp, xnames, Flow, Fupp, Fnames, &
          x, xstate, xmul, F, Fstate, Fmul, &
          INFO, mincw, miniw, minrw, &
          nS, nInf, sInf, &
          cw, lencw, iw, leniw, rw, lenrw, &
          cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'snOptA finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0) write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 30) go to 900
      stop

!     ------------------------------------------------------------------
!     Error exit.
!     ------------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  900 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'
  910 stop

 4000 format(/  a, 2x, a  )

end subroutine wsnopt

      subroutine t2dat0 ( Errors, maxF, maxn, Prob, nF, n, &
         ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul, &
          cw, lencw, iw, leniw, rw, lenrw )

      implicit none
      integer Errors, maxF, maxn, nF, n, ObjRow, lencw, leniw, lenrw, xstate(maxn), iw(leniw)
      double precision ObjAdd, xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF), x(maxn), Fmul(maxF), rw(lenrw)
      character Prob*8, cw(lencw)*8

!     ==================================================================
!     t2dat0  generates data for the test problem t2banana
!
!     On entry,
!     maxF, maxn are upper limits on nF and n.
!
!     On exit,
!     Errors    is 0 if there is enough storage, 1 otherwise.
!     nF        is the number of problem functions
!               (objective and constraints, linear and nonlinear).
!     n         is the number of variables.
!     xlow      holds the lower bounds on x.
!     xupp      holds the upper bounds on x.
!     Flow      holds the lower bounds on F.
!     Fupp      holds the upper bounds on F.

!     xstate(1:n)  is a set of initial states for each x  (0,1,2,3,4,5).
!     Fstate(1:nF) is a set of initial states for each F  (0,1,2,3,4,5).
!     x (1:n)      is a set of initial values for x.
!     Fmul(1:nF)   is a set of initial values for the dual variables.
!
!     19 Jul 2000: First version of t4manne based on SNOPT 5.3 t4manne.
!     27 Oct 2002: Current version.
!     ==================================================================
      integer   j, Obj
!     ------------------------------------------------------------------
      double precision   bplus,            bminus
      parameter         (bplus  = 1.0d+20, bminus = -bplus)
      double precision   zero,          ten
      parameter         (zero = 0.0d+0, ten = 10.0d+0)
!     ------------------------------------------------------------------
!     Name the Problem.

      Prob = 't2Banana'

      nF   = 1
      n    = 2

!     Check if there is enough storage.

      Errors = 0
      if (nF     .gt. maxF) Errors = 1
      if (n      .gt. maxn) Errors = 1
      if (Errors .gt.   0 ) return

      Obj    = nF
      ObjRow = Obj
      ObjAdd = zero             ! no additive obj. constant

!     ------------------------------------------------------------------
!     Ranges for the problem functions.
!     ------------------------------------------------------------------
!     Variable ranges

      do j = 1, n
         xlow(j) = -ten
         xupp(j) =  ten
      end do

      Flow(1) = bminus
      Fupp(1) = bplus

!     ------------------------------------------------------------------
!     Initialize x and xstate.
!     Set the initial value and status of each variable.
!     For want of something better to do, make the variables x(1:n)
!     temporarily fixed at their current values.
!     The crash can set the rest.
!     ------------------------------------------------------------------
      x(1) = - 1.2d+0
      x(2) =   1.0d+0

      do j = 1, n
         xstate(j) = 0
      end do

      end ! subroutine t2dat0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine userf ( Status, n, x, needF, nF, F, needG, &
       lenG, G, cu, lencu, iu, leniu, ru, lenru )

      implicit none
      integer Status, needF, needG, nF, n, lenG, lencu, leniu, lenru, iu(leniu)
      double precision  F(nF), G(lenG), x(n), ru(lenru)
      character cu(lencu)*8

!     ==================================================================
!     This is userfg for problem t2banana.
!     ==================================================================
      integer nOut, Obj
      double precision x1, x2
!     ------------------------------------------------------------------
      double precision   one
      parameter         (one = 0.0d+0)
!     ------------------------------------------------------------------
      nOut   = 6

!     First entry.  Print something on the standard output.

      if (Status .eq. 1) then
         if (nOut .gt. 0) write(nOut, 1000)
      end if

!     -------------
!     Normal entry.
!     -------------
!     Set the objective element.

      x1     =   x(1)
      x2     =   x(2)

      Obj    = 1
      F(Obj) = 100.0d+0*(x2 - x1**2)**2  + (one - x1)**2

!     ------------
!     Final entry.
!     ------------
      if (Status .ge. 2) then
         if (nOut .gt. 0) write(nOut, 2000)
      end if
      return

 1000 format(/ ' Starting  problem  t2banana.'/)
 2000 format(/ ' Finishing problem  t2banana.'/)

      end ! subroutine userf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine t2dat1  ( Errors, Prob, maxF, maxn, nF, n, &
          iAfun, jAvar, lenA, neA, A, &
          iGfun, jGvar, lenG, neG, &
          ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul, &
          cw, lencw, iw, leniw, rw, lenrw )

      implicit none
      integer Errors, maxF, maxn, neA, lenA, neG, lenG, nF, n, &
         ObjRow, lencw, leniw, lenrw, xstate(maxn), &
         iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG), &
         iw(leniw)
      double precision ObjAdd, &
         A(lenA), xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF), &
         x(maxn), Fmul(maxF), rw(lenrw)
      character Prob*8, cw(lencw)*8

!     ==================================================================
!     t2dat1  generates data for the test problem t2banana
!
!     On entry,
!     maxF, maxn are upper limits on nF and n.
!
!     The  nF by n  Jacobian is written as the sum of
!     two  nF by n  sparse matrices G and A, i.e.,  J = A + G,  where
!     A  and  G  contain contain the constant and nonlinear
!     elements of J respectively.
!
!     The nonzero pattern of G and A is specified by listing the
!     coordinates (i.e., row and column indices) of each nonzero.
!     Note that the coordinates specify the overall STRUCTURE of the
!     sparsity, not the Jacobian entries that happen to be zero at
!     the initial point.
!
!     The coordinates of the kth nonzero of  G  are defined
!     by iGfun(k) and jGvar(k)  (i.e., if i=iGfun(k) and j=jGvar(k),
!     then G(k) is the ijth element of G.)  Any known values of G(k)
!     must be assigned by the user in the routine userfg.
!
!     The coordinates of the kth nonzero of  A  are defined by
!     iAfun(k) and jAvar(k)  (i.e., if i=iAfun(k) and j=jAvar(k),
!     then A(k) is the ijth element of A.)  All values of A must be
!     assigned before the call to SNOPT.
!
!     The elements of A and G can be stored in any order, (e.g., by rows
!     or columns or mixed).
!
!     RESTRICTIONS:
!     1.  A nonzero entry of J must be specified as either an element
!         of A, or an element of G, but NOT BOTH (i.e.,  coordinates of
!         A  and  G  must not overlap.  Elements that are a sum of a
!         constant and varying part must be included in G and loaded
!         by userfg.
!
!     2.  If the computed value of an element of G happens to be zero
!         at a given point, it must still be loaded in userfg. (The
!         order of the coordinates is meaningful in SNOPT.)
!
!     On exit,
!     Errors    is 0 if there is enough storage, 1 otherwise.
!     nF        is the number of problem functions
!               (objective and constraints, linear and nonlinear).
!     n         is the number of variables.
!     neG       is the number of nonzeros in Jn.
!     neA       is the number of nonzeros in Jc.
!     xlow      holds the lower bounds on x.
!     xupp      holds the upper bounds on x.
!     Flow      holds the lower bounds on F.
!     Fupp      holds the upper bounds on F.

!     xstate(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
!     Fstate(1:m) is a set of initial states for each F  (0,1,2,3,4,5).
!     x (1:n)     is a set of initial values for x.
!     Fmul(1:m)   is a set of initial values for the dual variables.
!
!     19 Jul 2000: First version of t2banana based on SNOPT 5.3 t2banana.
!     31 Oct 2002: Current version.
!     ==================================================================
      integer     j, Obj
!     ------------------------------------------------------------------
      double precision   bplus,            bminus
      parameter         (bplus  = 1.0d+20, bminus = -bplus)
      double precision   zero,          five,          ten
      parameter         (zero = 0.0d+0, five = 5.0d+0, ten = 10.0d+0)
!     ------------------------------------------------------------------
!     Name the Problem.

      Prob = 't2Banana'

      nF   = 1
      n    = 2

!     Check if there is enough storage.

      Errors = 0
      if (nF     .gt. maxF) Errors = 1
      if (n      .gt. maxn) Errors = 1
      if (Errors .gt.   0 ) return

      Obj    = nF
      ObjRow = Obj

      neG    = 0
      neA    = 0

      neG    = neG + 1
      iGfun(neG) = Obj
      jGvar(neG) = 1
!     G(neG) = -400.0d+0*(x2 - x1**2)*x1 - 2.0d+0*(1.0d+0 - x1)

      neG    = neG + 1
      iGfun(neG) = Obj
      jGvar(neG) = 2
!     g(neG) =  200.0d+0*(x2 - x1**2)

      ObjAdd = zero             ! no additive obj. constant

!     ------------------------------------------------------------------
!     Ranges for the problem functions.
!     ------------------------------------------------------------------
!     Variable ranges

      do j = 1, n
         xlow(j) = -ten
         xupp(j) =  ten
      end do

      Flow(1) = bminus
      Fupp(1) = bplus

!     ------------------------------------------------------------------
!     Initialize x and xstate.
!     Set the initial value and status of each variable.
!     For want of something better to do, make the variables x(1:n)
!     temporarily fixed at their current values.
!     The crash can set the rest.
!     ------------------------------------------------------------------
      x(1) = - 1.2d+0
      x(2) =   1.0d+0

      do j = 1, n
         xstate(j) = 0
      end do

      end ! subroutine t2dat1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine userfg ( Status, n, x, &
          needF, nF, F, &
          needG, lenG, G, &
          cu, lencu, iu, leniu, ru, lenru )

      implicit none
      integer Status, needF, needG, nF, n, lenG, lencu, leniu, lenru, iu(leniu)
      double precision F(nF), G(lenG), x(n), ru(lenru)
      character cu(lencu)*8

!     ==================================================================
!     This is userfg for problem t2banana.
!     ==================================================================
      integer nOut, neG, Obj
      double precision x1, x2
!     ------------------------------------------------------------------
      double precision   one,          two
      parameter         (one = 0.0d+0, two = 2.0d+0)
!     ------------------------------------------------------------------
      nOut   = 6

!     First entry.  Print something on the standard output.

      if (Status .eq. 1) then
         if (nOut .gt. 0) write(nOut, 1000)
      end if

!     -------------
!     Normal entry.
!     -------------
!     Set the objective element.

      Obj    = 1

      x1     =   x(1)
      x2     =   x(2)

      if (needF .gt. 0) then
         F(Obj) = 100.0d+0*(x2 - x1**2)**2  + (one - x1)**2
      end if

      if (needG .gt. 0) then
         neG    = 1
!        iGfun(neG) = Obj
!        jGvar(neG) = 1
         G(neG) = -400.0d+0*(x2 - x1**2)*x1 - two*(one - x1)

         neG    = neG + 1
!        iGfun(neG) = Obj
!        jGvar(neG) = 2
         G(neG) =  200.0d+0*(x2 - x1**2)

      end if

!     ------------
!     Final entry.
!     ------------
      if (Status .ge. 2) then
         if (nOut .gt. 0) write(nOut, 2000)
      end if
      return

 1000 format(/ ' Starting  problem  t2banana.'/)
 2000 format(/ ' Finishing problem  t2banana.'/)

      end ! subroutine userfg
