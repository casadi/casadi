//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file LinMods.cpp
 * \brief Implement the Modification class LinMods, that is used to store
 * modifications to the linear parts and rhs of a constraint.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <iostream>

#include "MinotaurConfig.h"
#include "LinMods.h"

using namespace Minotaur;


LinMods::LinMods()
{
  bmods_.clear();
  bmods2_.clear();
  lmods_.clear();
}


LinMods::~LinMods()
{
  bmods_.clear();
  bmods2_.clear();
  lmods_.clear();
}


void LinMods::applyToProblem(ProblemPtr problem) 
{
  for (VarBoundModConstIter it = bmods_.begin(); it != bmods_.end(); ++it) {
    (*it)->applyToProblem(problem);
  }
  for (VarBoundMod2ConstIter it = bmods2_.begin(); it != bmods2_.end(); ++it) {
    (*it)->applyToProblem(problem);
  }
  for (LinConModConstIter it = lmods_.begin(); it != lmods_.end(); ++it) {
    (*it)->applyToProblem(problem);
  }
}


void LinMods::undoToProblem(ProblemPtr problem) 
{
  for (VarBoundModConstIter it = bmods_.begin(); it != bmods_.end(); ++it) {
    (*it)->undoToProblem(problem);
  }
  for (VarBoundMod2ConstIter it = bmods2_.begin(); it != bmods2_.end(); ++it) {
    (*it)->undoToProblem(problem);
  }
  for (LinConModConstIter it = lmods_.begin(); it != lmods_.end(); ++it) {
    (*it)->undoToProblem(problem);
  }
}


void LinMods::insert(VarBoundModPtr bmod)
{
  bmods_.push_back(bmod);
}


void LinMods::insert(VarBoundMod2Ptr bmod2)
{
  bmods2_.push_back(bmod2);
}


void LinMods::insert(LinConModPtr lmod)
{
  lmods_.push_back(lmod);
}


bool LinMods::isEmpty() const
{
  return (lmods_.empty() && bmods2_.empty() && bmods_.empty());
}


void LinMods::write(std::ostream &out) const
{
  out << "LinMods: " << std::endl;
  for (VarBoundModConstIter it = bmods_.begin(); it != bmods_.end(); ++it) {
    (*it)->write(out);
  }
  for (VarBoundMod2ConstIter it = bmods2_.begin(); it != bmods2_.end(); ++it) {
    (*it)->write(out);
  }
  for (LinConModConstIter it = lmods_.begin(); it != lmods_.end(); ++it) {
    (*it)->write(out);
  }
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
