//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file YEqMonomial.cpp
 * \brief Define class for storing auxiliary variables for monomial
 * expressions.
 * \author Ashutosh Mahajan, Argonne National Laboratory.
 */

#include <cmath>

#include "MinotaurConfig.h"

#include "PolynomialFunction.h"
#include "Variable.h"
#include "YEqMonomial.h"

using namespace Minotaur;

YEqMonomial::YEqMonomial(UInt n)
{
  n_ = n;
  for (UInt i=0; i<n_; ++i) {
    rand_.push_back((double) rand()/(RAND_MAX));
  }
}


double YEqMonomial::evalHash_(MonomialFunPtr mf)
{
  double hash = 0.0;
  for (VarIntMapConstIterator it=mf->termsBegin(); it!=mf->termsEnd(); ++it) {
    hash += ((double) it->second)/(1.0+it->first->getId());
  }
  return hash;
}


VariablePtr YEqMonomial::findY(MonomialFunPtr mf)
{
  double hash = evalHash_(mf);
  VarIntMapConstIterator it, it2;
  bool found;

  for (UInt i=0; i<y_.size(); ++i) {
    if (fabs(hash-hash_[i])<1e-12 
        && mf->getDegree()==mf_[i]->getDegree()
        && fabs(mf->getCoeff()-mf_[i]->getCoeff())<1e-12) {
      it = mf->termsBegin();
      it2 = mf_[i]->termsBegin();
      found = true;
      for (; it!=mf->termsEnd(); ++it, ++it2) {
        if (it->first != it2->first || it->second != it2->second) {
          found = false;
          break;
        }
      }
      if (found) {
        return y_[i];
      }
    }
  }
  return VariablePtr();
}


void YEqMonomial::insert(VariablePtr auxvar, MonomialFunPtr mf)
{
  hash_.push_back(evalHash_(mf));
  y_.push_back(auxvar);
  mf_.push_back(mf);
}


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
