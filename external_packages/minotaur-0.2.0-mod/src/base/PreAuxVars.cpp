//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file PreAuxVars.cpp
 * \brief Define the PreAuxVars class for saving additional variables
 * variables added during presolve.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include "PreAuxVars.h"
using namespace Minotaur;


PreAuxVars::PreAuxVars() 
{
}


PreAuxVars::~PreAuxVars() 
{
  vars_.clear();
} 


void PreAuxVars::insert(VariablePtr v)
{
  vars_.push_back(v);
}


void PreAuxVars::postsolveGetX(const DoubleVector &x, DoubleVector *newx) 
{
  UInt n = x.size()-vars_.size();
  UInt i, j;
  std::deque<VariablePtr>::iterator it = vars_.begin();
  assert(x.size()>=vars_.size());
  newx->resize(n);

  i = 0; j = 0;
  for (j=0; j<n; ++i) {
    if (it!=vars_.end() && i==(*it)->getIndex()) {
      ++it;
    } else {
      (*newx)[j] = x[i];
      ++j;
    }
  }
}


UInt PreAuxVars::getSize() 
{
  return vars_.size();
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
