module snopt_wrapper

use iso_c_binding

implicit none

interface
  subroutine userfun &
      ( mode, nnObj, nnCon, nnJac, nnL, neJac, x, fObj, gObj,&
        fCon, gCon, nState, cu, lencu, iu, leniu, ru, lenru ) bind(c)
  import :: c_int, c_double, c_char
  integer(c_int), intent(in) :: lencu, leniu, lenru, nnObj, nnCon, nnJac, nnL, neJac, nState, iu(leniu)
  integer(c_int), intent(inout) :: mode
  real(c_double), intent(inout) :: fObj, fCon(nnCon), gCon(neJac), gObj(nnObj), ru(lenru), x(nnL)
  character(c_char), intent(inout) ::  cu(lencu*8)
  end subroutine
  !subroutine snSeti(buffer, ivalue, iPrint, iSumm,Errors, cw, lencw, iw, leniw, rw, lenrw)
  !  integer, intent(in)  :: lencw, leniw, lenrw, ivalue
  !  integer, intent(inout) :: iPrint, iSumm, Errors
  !  integer, intent(inout) :: iw(leniw)
  !  character*(*) , intent(in):: buffer
  !  real, intent(inout) :: rw(lenrw)
  !end subroutine
  
  subroutine snSTOP &
     &   ( iAbort, info, HQNType, KTcond, MjrPrt, minimz,
     &     m, maxS, n, nb, nnCon0, nnCon, nnObj0, nnObj, nS,
     &     itn, nMajor, nMinor, nSwap,
     &     condHz, iObj, sclObj, ObjAdd, fMrt, PenNrm, step,
     &     prInf, duInf, vimax, virel, hs,
     &     ne, nlocJ, locJ, indJ, Jcol, negCon,
     &     Ascale, bl, bu, fCon, gCon, gObj,
     &     yCon, pi, rc, rg, x,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw ) bind(c)
     
     import :: c_int, c_double, c_char
     
     logical(c_bool), intent(in) ::  KTcond(2)
     
     integer(c_int), intent(in) :: HQNType, info(6), iObj, itn,
     &     lencu, lencw, leniu, leniw, lenru, lenrw,
     &     MjrPrt, minimz, m, maxS, n, nb, ne, negCon, nlocJ,
     &     nnCon0, nnCon, nnObj0, nnObj, nMajor, nMinor, nS, nSwap,
     &     hs(nb), locJ(nlocJ), indJ(ne), iu(leniu), iw(leniw)
     
     real(c_double), intent(in) :: condHz, sclObj, ObjAdd, fMrt, PenNrm, virel, vimax, step,
     &     prInf, duInf, Ascale(nb), bl(nb), bu(nb), fCon(nnCon0),
     &     gCon(negCon), gObj(nnObj0), Jcol(ne), pi(m),
     &     rc(nb), rg(maxS), yCon(nnCon0), x(nb), ru(lenru), rw(lenrw)
     
     integer(c_int), intent(inout) :: iAbort
     
     character(c_char), intent(inout) ::  cu(lencu*8)

  end subroutine

end interface

contains

subroutine snopt_init (iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw) bind (c)

    integer(c_int), intent(in) :: iPrint, iSumm, lencw, leniw, lenrw

    integer(c_int), intent(inout) :: iw(leniw)
    character(c_char), intent(inout)  :: cw(lencw*8)
    real(c_double), intent(inout) :: rw(lenrw)

    call snInit ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

end  subroutine

subroutine snopt_getc (buffer, lenbuffer, cvalue, Errors, cw, lencw, iw, leniw, rw, lenrw) bind (c)
    integer(c_int), intent(in) :: lencw, leniw, lenrw, lenbuffer
    integer(c_int), intent(inout) :: iw(leniw), Errors
    character(c_char), intent(inout)  :: cw(lencw*8), cvalue(8)
    character(c_char), intent(in) :: buffer(lenbuffer)
    real(c_double), intent(inout) :: rw(lenrw)

    character(lenbuffer), pointer :: str
    type(C_PTR) p
    p = transfer(LOC(buffer(1)),p)
    call C_F_POINTER(p,str)

    call snGetc(buffer, cvalue, Errors, cw, lencw, iw, leniw, rw, lenrw)
end subroutine

subroutine snopt_geti (buffer, lenbuffer, ivalue, Errors, cw, lencw, iw, leniw, rw, lenrw) bind (c)
    integer(c_int), intent(in) :: lencw, leniw, lenrw, lenbuffer
    integer(c_int), intent(inout) :: ivalue
    integer(c_int), intent(inout) :: iw(leniw), Errors
    character(c_char), intent(inout)  :: cw(lencw*8)
    character(c_char), intent(in) :: buffer(lenbuffer)
    real(c_double), intent(inout) :: rw(lenrw)

    character(lenbuffer), pointer :: str
    type(C_PTR) p
    p = transfer(LOC(buffer(1)),p)
    call C_F_POINTER(p,str)

    call snGeti(str, ivalue, Errors, cw, lencw, iw, leniw, rw, lenrw)
end subroutine

subroutine snopt_getr (buffer, lenbuffer, rvalue, Errors, cw, lencw, iw, leniw, rw, lenrw) bind (c)
    integer(c_int), intent(in) :: lencw, leniw, lenrw, lenbuffer
    real(c_double), intent(inout) :: rvalue
    integer(c_int), intent(inout) :: iw(leniw), Errors
    character(c_char), intent(inout)  :: cw(lencw*8)
    character(c_char), intent(in) :: buffer(lenbuffer)
    real(c_double), intent(inout) :: rw(lenrw)

    character(lenbuffer), pointer :: str
    type(C_PTR) p
    p = transfer(LOC(buffer(1)),p)
    call C_F_POINTER(p,str)

    call snGetr(str, rvalue, Errors, cw, lencw, iw, leniw, rw, lenrw)
end subroutine

subroutine snopt_set (buffer, lenbuffer,iPrint, iSumm, Errors, cw, lencw, iw, leniw, rw, lenrw) bind (c)
    integer(c_int), intent(in) :: lencw, leniw, lenrw, lenbuffer
    integer(c_int), intent(inout) :: iw(leniw),iPrint, iSumm, Errors
    character(c_char), intent(inout)  :: cw(lencw*8)
    character(c_char), intent(in) :: buffer(lenbuffer)
    real(c_double), intent(inout) :: rw(lenrw)

    character(lenbuffer), pointer :: str
    type(C_PTR) p
    p = transfer(LOC(buffer(1)),p)
    call C_F_POINTER(p,str)

    call snSet(str, iPrint, iSumm,Errors, cw, lencw, iw, leniw, rw, lenrw)
end subroutine

subroutine snopt_seti (buffer, lenbuffer, ivalue, iPrint, iSumm,Errors, cw, lencw, iw, leniw, rw, lenrw) bind (c)
    integer(c_int), intent(in) :: lencw, leniw, lenrw, lenbuffer, ivalue
    integer(c_int), intent(inout) :: iPrint, iSumm, Errors
    integer(c_int), intent(inout) :: iw(leniw)
    character(c_char), intent(inout)  :: cw(lencw*8)
    character(c_char), intent(in) :: buffer(lenbuffer)
    real(c_double), intent(inout) :: rw(lenrw)

    character(lenbuffer), pointer :: str
    type(C_PTR) p
    p = transfer(LOC(buffer(1)),p)
    call C_F_POINTER(p,str)

    call snSeti(str, ivalue, iPrint, iSumm,Errors, cw, lencw, iw, leniw, rw, lenrw)

end subroutine


subroutine snopt_setr (buffer, lenbuffer, rvalue, iPrint, iSumm,Errors, cw, lencw, iw, leniw, rw, lenrw) bind (c)
    integer(c_int), intent(in) :: lencw, leniw, lenrw, lenbuffer
    real(c_double), intent(inout) :: rvalue
    integer(c_int), intent(inout) :: iw(leniw), iPrint, iSumm, Errors
    character(c_char), intent(inout)  :: cw(lencw*8)
    character(c_char), intent(in) :: buffer(lenbuffer)
    real(c_double), intent(inout) :: rw(lenrw)

    character(lenbuffer), pointer :: str
    type(C_PTR) p
    p = transfer(LOC(buffer(1)),p)
    call C_F_POINTER(p,str)

    call snSetr(str, rvalue, iPrint, iSumm,Errors, cw, lencw, iw, leniw, rw, lenrw)
end subroutine

subroutine snopt_spec (iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw) bind (c)
    integer(c_int), intent(in) :: iSpecs, lencw, leniw, lenrw
    integer(c_int), intent(inout) :: INFO
    integer(c_int), intent(inout) :: iw(leniw)
    character(c_char), intent(inout)  :: cw(lencw*8)
    real(c_double), intent(inout) :: rw(lenrw)
    call snSpec(iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw)
end subroutine

subroutine snopt_memb (INFO, m,n,neA, negCon, nnCon, nnJac, nnObj, mincw, miniw,minrw, cw, lencw, iw, leniw, rw, lenrw) bind (c)
    integer(c_int), intent(in) :: m,n,neA, negCon, nnCon, nnJac, nnObj, lencw, leniw, lenrw
    integer(c_int), intent(inout) :: INFO
    integer(c_int), intent(inout) :: iw(leniw)
    character(c_char), intent(inout)  :: cw(lencw*8)
    real(c_double), intent(inout) :: rw(lenrw)

    integer(c_int), intent(inout) :: mincw, miniw, minrw

    call snMemB(INFO, m,n,neA, negCon, nnCon, nnJac, nnObj, mincw, miniw,minrw, cw, lencw, iw, leniw, rw, lenrw)
end subroutine


subroutine snopt_c &
(Start, lenstart, m, n, neA, nName, nnCon, nnObj, nnJac, iObj,ObjAdd, Prob,userfun, snStop, Acol, indA, locA, bl, bu, Names, hs, x,&
pi, rc, INFO, mincw, miniw, minrw, nS, nInf, sInf, Obj,&
cu, lencu, iu, leniu, ru, lenru,&
cw, lencw, iw, leniw, rw, lenrw) bind (c)
    type(c_funptr) userfun
    external snLog, snLog2, sqLog
    type(c_funptr) snSTOP
    integer(c_int), intent(in) :: lenstart, INFO, iObj, lencu, leniu, lenru, lencw, leniw, lenrw, &
      mincw, miniw, minrw, m, n, neA, nName, nS, nInf, nnCon, nnObj, nnJac
    integer(c_int), intent(in) :: indA(neA), hs(n+m), locA(n+1)
    integer(c_int), intent(inout) :: iu(leniu), iw(leniw)
    real(c_double), intent(in) :: Obj, ObjAdd, sInf, Acol(neA), bl(n+m), bu(n+m)
    real(c_double), intent(inout) :: pi(m), rc(n+m), ru(lenru), rw(lenrw), x(n+m)
    character(c_char), intent(inout) :: Prob(8), Names(nName*8), cu(lencu*8), cw(lencw*8)
    character(c_char), intent(in) :: Start(lenstart)

    call snKerC &
      (Start, m, n, neA, nName, nnCon, nnObj, nnJac, iObj,ObjAdd, Prob, userfun, snSTOP, Acol, indA,&
       locA, bl, bu, Names, hs, x,&
       pi, rc, INFO, mincw, miniw, minrw, nS, nInf, sInf, Obj, cu, lencu, iu,&
       leniu, ru, lenru, cw, lencw, iw, leniw, rw, lenrw)
end subroutine


end module
