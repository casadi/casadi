// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

#include "MinotaurConfig.h"
#include "Chol.h"

using namespace Minotaur;

#ifdef F77_FUNC

extern "C"
{
  void F77_FUNC(dpotrf, DPOTRF )(char *uplo, int *n, double *a, int *lda, 
                                 int *info );
}


CholCalculator::CholCalculator()
  :n_(0),
   A_(0),
   abstol_(1e-6)
{
}

#endif

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
