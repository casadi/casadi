//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file PreDelVars.cpp
 * \brief Postsolver for substituted variables.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include "MinotaurConfig.h"
#include "PreDelVars.h"

using namespace Minotaur;

PreDelVars::PreDelVars()
{
  vars_.clear();
}


PreDelVars::~PreDelVars()
{
  vars_.clear();
}


void PreDelVars::insert(VariablePtr v)
{
  vars_.push_front(v);
}


void PreDelVars::postsolveGetX(const DoubleVector &x, DoubleVector *newx)
{
  UInt n = x.size()+vars_.size();
  bool *filled = new bool[n];
  VariablePtr v;
  DoubleVector::const_iterator dit;
  DoubleVector::iterator dit2;

  std::fill(filled, filled+n, false);
  newx->resize(n);
  for (VarQueueConstIter it=vars_.begin(); it!=vars_.end(); ++it) {
    v = (*it);
    (*newx)[v->getIndex()] = v->getLb();
    filled[v->getIndex()] = true;
  }


  dit=x.begin();
  dit2=newx->begin();
  for (bool *it=filled; it!=filled+n; ++it, ++dit2) {
    if ((*it) == false) {
      *dit2 = *dit;
      ++dit;
    }
  }
  delete [] filled;
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
