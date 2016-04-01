//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file Transformer.cpp
 * \brief Define class for reformulating a problem suitable for handlers.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iostream>

#include "MinotaurConfig.h"

#include "Constraint.h"
#include "LinearFunction.h"
#include "Variable.h"
#include "YEqLFs.h"

// #define SPEW 1

using namespace Minotaur;


YEqLFs::YEqLFs(UInt n)
  :n_(n)
{
  for (UInt i=0; i<n; ++i) {
    rand_.push_back((double) rand()/(RAND_MAX)*10.0);
  }
}


double YEqLFs::evalHash_(LinearFunctionPtr lf)
{
  double hash = 0.0;
  for (VariableGroupConstIterator it=lf->termsBegin(); it!=lf->termsEnd();
       ++it) {
    hash += it->second*rand_[it->first->getId()%n_];
  }
  return hash;
}


VariablePtr YEqLFs::findY(LinearFunctionPtr lf, double k)
{
  bool found;
  VariableGroupConstIterator it, it2;
  double hash = evalHash_(lf);
  for (UInt i=0; i<y_.size(); ++i) {
    if (fabs(k-k_[i])<1e-12 && lf->getNumTerms()==lf_[i]->getNumTerms() && 
        fabs(hash-hash_[i])<1e-12) {
      found = true;
      it = lf->termsBegin();
      it2 = lf_[i]->termsBegin();
      for (; it!=lf->termsEnd(); ++it, ++it2) {
        if (it->first!=it2->first || fabs(it->second-it2->second)>1e-12) {
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


void YEqLFs::insert(VariablePtr auxvar, LinearFunctionPtr lf, double k)
{
  hash_.push_back(evalHash_(lf));
  lf_.push_back(lf);
  y_.push_back(auxvar);
  k_.push_back(k);
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
