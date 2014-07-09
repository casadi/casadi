/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/// \cond INTERNAL

#pragma once

typedef void (*UserFun)(
    int * mode, int* nnObj, int * nnCon, int *nJac, int *nnL, int * neJac,
    double *x, double *fObj, double *gObj, double * fCon, double* gCon, int* nState,
    char* cu, int* lencu, int* iu, int* leniu, double* ru, int *lenru);

typedef void (*snStop)(
    int* iAbort, int* info, int* HQNType, int* KTcond, int* MjrPrt, int* minimz,
    int* m, int* maxS, int* n, int* nb,
    int* nnCon0, int* nnCon, int* nnObj0, int* nnObj, int* nS,
    int* itn, int* nMajor, int* nMinor, int* nSwap,
    double * condHz, int* iObj, double * sclObj,  double *ObjAdd,
    double * fMrt,  double * PenNrm,  double * step,
    double *prInf,  double *duInf,  double *vimax,  double *virel, int* hs,
    int* ne, int* nlocJ, int* locJ, int* indJ, double* Jcol, int* negCon,
    double* Ascale, double* bl, double* bu, double* fCon, double* gCon, double* gObj,
    double* yCon, double* pi, double* rc, double* rg, double* x,
    double*  cu, int * lencu, int* iu, int* leniu, double* ru, int *lenru,
    char*   cw, int* lencw,  int* iw, int *leniw, double* rw, int* lenrw);

typedef void (*dummyFun)();

#ifdef __cplusplus
extern "C" {
#endif

/*
  extern void snopt_init(const int * iPrint, const int * iSumm,
    char * cw, const int * lencw, int *iw, const  int *leniw, double * rw, const int * lenrw);

  extern void snopt_c(
    const char * Start, const int * lenstart, const int * m, const int * n, const int * neA,
    const int * nName, const int *nnCon, const int *nnObj, const int *nnJac, const int *iObj,
    const double *ObjAdd, const char* Prob , UserFun userfun,

    const double* Acol, const int* indA, const int *locA, double* bl, double* u,

    char* Names,

    // Initial values
    int* hs, double* x, double* pi, double * rc,

    // Outputs
    int *INFO, int* mincw, int* miniw, int* minrw, int * nS, int* nInf, double* sInf, double* Obj,

    // Working spaces for usrfun
    char* cu, const int* lencu, int* iu, const int* leniu, double* ru, const int* lenru,
    // Working spaces for SNOPT
    char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);

  // cvalue is always length 8
  extern void snopt_getc(const char *buffer, const int* lenbuffer, char *cvalue, int* Errors,
    char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);
    //  extern void snopt_geti(const char *buffer, const int* lenbuffer, int *ivalue, int* Errors,
    //    char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);

  extern void snopt_getr(const char *buffer, const int* lenbuffer, double *rvalue, int* Errors,
    char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);

  extern void snopt_set(const char *buffer, const int* lenbuffer, int * iPrint, int * iSumm,
    int* Errors,
    char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);
    //extern void snopt_seti(const char *buffer, const int* lenbuffer, const int *ivalue,
    //  int * iPrint, int * iSumm, int* Errors,
    //  char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);

  extern void snopt_setr(const char *buffer, const int* lenbuffer, const double *rvalue,
    int * iPrint, int * iSumm, int* Errors,
    char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);


  extern void snopt_spec(const int *iSpecs, int* INFO,
    char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);

  extern void snopt_memb(int *INFO, const int* m, const int* n, const int* neA,
    const int*  negCon, const int* nnCon, const int* nnJac, const int*  nnObj,
    int* mincw, int* miniw, int* minrw,
    char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw);
*/

  /* direct calls to the Fortran library */
  void sninit_(const int * iPrint, const int * iSumm,
               char* cw, const int* lencw,
               int* iw, const int* leniw,
               double* rw, const int* lenrw,
               const long cw_len8);
  void snseti_(const char *buffer, const int *ivalue, int * iPrint, int * iSumm, int* Errors,
               char* cw, const int* lencw,
               int* iw, const int* leniw,
               double* rw, const int* lenrw,
               const long buffer_ftn_len, const long cw_len8);
  void snsetr_(const char *buffer, const double *ivalue, int * iPrint, int * iSumm, int* Errors,
               char* cw, const int* lencw,
               int* iw, const int* leniw,
               double* rw, const int* lenrw,
               const long buffer_ftn_len, const long cw_len8);
  void snset_(const char *buffer, int * iPrint, int * iSumm, int* Errors,
              char* cw, const int* lencw,
              int* iw, const int* leniw,
              double* rw, const int* lenrw,
              const long buffer_ftn_len, const long cw_len8);
  void snmemb_(int *INFO, const int* m, const int* n, const int* neA, const int*  negCon,
               const int* nnCon, const int* nnJac, const int*  nnObj,
               int* mincw, int* miniw, int* minrw,
               char* cw, const int* lencw,
               int* iw, const int* leniw,
               double* rw, const int* lenrw,
               const long cw_len8);
  void snoptc_(const char * Start, const int * m, const int * n, const int * neA,
               const int * nName, const int *nnCon, const int *nnObj, const int *nnJac,
               const int *iObj, const double *ObjAdd, const char* Prob , UserFun userfun,
               const double* Acol, const int* indA, const int *locA, double* bl, double* bu,
               char* Names,
               /* Initial values */
               int* hs, double* x, double* pi, double * rc,
               /* Outputs */
               int *info, int* mincw, int* miniw, int* minrw, int * nS,
               int* nInf, double* sInf, double* Obj,
               /* Working spaces for usrfun */
               char* cu, const int* lencu, int* iu, const int* leniu, double* ru, const int* lenru,
               /* Working spaces for SNOPT */
               char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw,
               /* Fortran char array hack */
               const long start_len8, const long prob_len8, const long names_len8,
               const long cu_len8, const long cw_len8);

  void snkerc_(const char * Start, const int * m, const int * n, const int * neA,
               const int * nName, const int *nnCon, const int *nnObj, const int *nnJac,
               const int *iObj, const double *ObjAdd, const char* Prob ,
               UserFun userfun, dummyFun snlog, dummyFun snlog2, dummyFun sqlog, snStop snstop,
               const double* Acol, const int* indA, const int *locA, double* bl, double* bu,
               char* Names,
               /* Initial values */
               int* hs, double* x, double* pi, double * rc,
               /* Outputs */
               int *info, int* mincw, int* miniw, int* minrw, int * nS,
               int* nInf, double* sInf, double* Obj,
               /* Working spaces for usrfun */
               char* cu, const int* lencu, int* iu, const int* leniu, double* ru, const int* lenru,
               /* Working spaces for SNOPT */
               char* cw, const int* lencw, int* iw, const int* leniw, double* rw, const int* lenrw,
               /* Fortran char array hack */
               const long start_len8, const long prob_len8, const long names_len8,
               const long cu_len8, const long cw_len8);

  // Dummy symbols that point to functions in the snopt library
  void snlog_();
  void snlog2_();
  void sqlog_();

#ifdef __cplusplus
}
#endif

/// \endcond
