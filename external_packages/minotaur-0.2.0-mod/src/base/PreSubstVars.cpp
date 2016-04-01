//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file PreSubstVars.cpp
 * \brief Postsolver for substituted variables.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include "MinotaurConfig.h"
#include "PreSubstVars.h"

using namespace Minotaur;

PreSubstVars::PreSubstVars()
{
  vars_.clear();
}


PreSubstVars::~PreSubstVars()
{
  for (std::deque<PreSubstVarData *>::const_iterator 
      it=vars_.begin(); it!=vars_.end(); ++it) {
    delete (*it);
  }
  vars_.clear();
}


void PreSubstVars::insert(VariablePtr vout, VariablePtr vin, double rat)
{
  PreSubstVarData *data = new PreSubstVarData();
  data->vout = vout;
  data->vinInd = vin->getIndex();
  data->rat = rat;
  vars_.push_front (data);
}


void PreSubstVars::postsolveGetX(const DoubleVector &x, DoubleVector *newx)
{
  // always called after PreDelVars::postsolveGetX(), so don't worry about
  // copying all values in x to newx.

  for (std::deque<PreSubstVarData *>::const_iterator 
      it=vars_.begin(); it!=vars_.end(); ++it) {
    (*newx)[(*it)->vout->getIndex()] = x[(*it)->vinInd] * (*it)->rat;
  }

  // Is important to do it twice. e.g let original variables be [w, x, y, z].
  // Suppose we substuted z by y, y by x, x by w in that order. After the
  // above loop, only the value of x has been restored, not y and z. Note that
  // the difference between the above loop and the below is that the rhs in
  // assignment is different.
  for (std::deque<PreSubstVarData *>::const_iterator 
      it=vars_.begin(); it!=vars_.end(); ++it) {
    (*newx)[(*it)->vout->getIndex()] = (*newx)[(*it)->vinInd] * (*it)->rat;
  }
}


UInt PreSubstVars::getSize()
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
