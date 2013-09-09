!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Example problem  HS118
!   The Hessian is given in a user-defined subroutine, which computes H*x
!   for a given vector x.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine wsqic (m, n, nnzA, indA, locA, valA, bl, bu, hEtype, hs, x, pi, rc, nnzH, indH, locH, valH, Obj) bind ( C, name="sqic" )
  use snModulePrecision, only : ip, rp
  use SQIC,              only : qpProb                     

  implicit none

  integer                   :: INFO
  integer(ip)               :: Errors, iObj, iPrint, iSumm, ncObj, &
                               m, n, nInf, nnH, nnzH, nNames, nnzA, nS

  real(rp)                  :: Obj, ObjAdd, sInf

  real(rp),     allocatable :: cObj(:)
  real(rp):: bl(n+m), bu(n+m), x(n+m), valA(nnzA), valH(nnzH) ,pi(m), rc(n+m)
  integer(ip):: indA(nnzA), locA(n+1), indH(nnzH), locH(n+1), hEtype(n+m), hs(n+m)

  character(8), allocatable :: Names(:)
  character(8)              :: probName

  type(qpProb)              :: QP

  Errors   = 0

  ncObj    =  0
  nnH      = n
  nNames   =  1
  ObjAdd = 0.0

  iSumm    = 6;  iPrint = 9;
  open ( iPrint, file='hs118.out', status='unknown' )

  ! Allocate space for problem.
  allocate ( cObj(ncObj) )
  allocate ( Names(nNames) )

  probName = 'QP'
  call QP%load ( probName, m, n, nnH, m, ObjAdd, &
                 nnzA, indA, locA, valA, bl, bu, ncObj, cObj, &
                 nNames, Names, hEtype, hs, x, pi, rc, nnzH, indH, locH, valH )

  ! Initialize SQIC.
  call QP%begin ( iPrint, iSumm )

  ! Set options.
  call qp%set   ( 'Print level        1', iPrint, iSumm, Errors )
  call qp%set   ( 'Summary level      1', iPrint, iSumm, Errors )
  call QP%set   ( 'Print frequency    1', iPrint, iSumm, Errors )
  call QP%set   ( 'Summary frequency  1', iPrint, iSumm, Errors )


  ! Solve the QP.
  call QP%solve ( 'Cold', INFO, nS, nInf, sInf, Obj )
  
  ! Finish up and deallocate.
  call QP%end

  if ( allocated(cObj) )   deallocate ( cObj )
  if ( allocated(Names) )  deallocate ( Names )

  close       ( iPrint )

end subroutine wsqic
