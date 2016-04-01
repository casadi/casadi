// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

/**
 * \file FilterSQPEngineTypes.h
 * \brief Define the types (real, fint, and others) that are used by FilterSQP
 * solver.
 * \author Sven Leyffer, Argonne National Laboratory
 */

#ifndef MINOTAURFILTERSQPENGINETYPES_H
#define MINOTAURFILTERSQPENGINETYPES_H

#include <stdint.h>

#include "Types.h"

using namespace Minotaur;

#if 1 // on linux, with g77, this works.

typedef double real;
typedef int    fint;
#define objfun objfun_
#define confun confun_
#define gradient gradient_
#define objgrad objgrad_
#define hessian hessian_

#else  // if we are on some other architecture, use some other defines etc.

#endif // end of if 1.

int * convertPtrToInt(uintptr_t u);

uintptr_t convertIntToPtr(int *iarray);

extern "C" {
void objfun(real *x, fint *n, real *f, real *user, fint *iuser, 
    fint *errflag);

void confun(real *x, fint *n, fint *m, real *c, real *a, fint *la,
    real *user, fint *iuser, fint *errflag);

void gradient(fint *N, fint *M, fint *mxa, real *x, real *a, fint *la,
    fint *maxa, real *user, fint *iuser, fint *errflag);

void objgrad (fint *n, fint *m, fint *mxa, real *x, real *a, fint *la, fint
    *maxa, real *user, fint *iuser, fint *errflag);

void hessian(real *x, fint *N, fint *M, fint *phase, real *lam,
    real *ws, fint *lws, real *user, fint *iuser,
    fint *l_hess, fint *li_hess, fint *errflag);

void filtersqp_ ( fint *n, fint *m, fint *kmax, fint *maxa, fint *maxf, 
    fint *mlp, fint *mxwk, fint *mxiwk, fint *iprint, fint *nout, fint *ifail, 
    real *rho, real *x, real *c, real *f, real *fmin, real *bl, real *bu, 
    real *s, real *a, fint *la, real *ws, fint *lws, real *lam, char *cstype, 
    real *user, fint *iuser, fint *maxiter, fint *istat, real *rstat, 
    long cstype_len);
}

#endif // ifndef MINOTAURFILTERSQPENGINETYPES_H

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
