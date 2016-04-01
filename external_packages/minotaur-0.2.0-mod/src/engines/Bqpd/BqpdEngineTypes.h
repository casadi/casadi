// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

// /**
// \file BqpdEngineTypes.h
// \brief Define the types (real, fint, and others) that are used by bqpd
// solver.
// \author Sven Leyffer, Argonne National Laboratory
//
// */

#include <stdint.h>

#include "Types.h"

#ifndef MINOTAURBQPDENGINETYPES_H
#define MINOTAURBQPDENGINETYPES_H

using namespace Minotaur;

#if 1 // on linux, with g77, this works.

typedef double real;
typedef int    fint;
#define gdotx gdotx_
#define setwsc setwsc_
#define writewsc writewsc_
#define restorecommon restorecommon_
#define savecommon savecommon_
#define defaultcommon defaultcommon_

#else  // if we are on some other architecture, use some other defines etc.

#endif // end of if 1.

/// \todo explain this
int * convertPtrToInt(uintptr_t u);

/// \todo explain this
uintptr_t convertIntToPtr(int *iarray);

extern "C" {

  /// copy storage map for bqpd
  void setwsc(fint *mxws0, fint *mxlws0, fint *kk0, fint *ll0);
  void writewsc();
  void restorecommon();
  void savecommon();
  void defaultcommon();

  /// compute v=G.x (Hessian vector product); G stored in ws&lws
  void gdotx(fint *n, real *x, real *ws, fint *lws, real *v);

  /// call bqpd (QP) solver
  void bqpd_(fint *n, fint *m, fint *k, fint *kmax, real *a, fint *la, real *x, 
             real *bl, real *bu, real *f, real *fmin, real *g, real *r, 
             real *w, real *e, fint *ls, real *alp, fint *lp, fint *mlp,
             fint *peq, real *ws, fint *lws, fint *mode, fint *ifail, 
             fint *info, fint *iprint, fint *nout);
}

#endif // ifndef MINOTAURBQPDENGINETYPES_H

// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
