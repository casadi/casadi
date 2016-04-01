// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file Branch.cpp
 * \author Ashutosh Mahajan, Argonne National Laboratory.
 * \brief Define the class Branch for describing branches in branch-and-bound.
 */

#include <cmath>
#include <iostream>

#include "MinotaurConfig.h"
#include "BrCand.h"
#include "Branch.h"
#include "Modification.h"
#include "Operations.h"
#include "Variable.h"

using namespace Minotaur;

const std::string Branch::me_ = "Branch: "; 

Branch::Branch()
: pMods_(0),
  rMods_(0),
  activity_(INFINITY),
  brCand_(BrCandPtr()) // NULL
{

}


Branch::~Branch()
{
  pMods_.clear();
  rMods_.clear();
}


void Branch::addPMod(ModificationPtr mod) 
{
  pMods_.push_back(mod);
}


void Branch::addRMod(ModificationPtr mod) 
{
  rMods_.push_back(mod);
}


double Branch::getActivity() const
{
  return activity_;
}


void Branch::setActivity(double value) 
{
  activity_ = value;
}

void Branch::write(std::ostream &out) const
{
  out << me_ << "Problem modifications:" << std::endl;
  for (ModificationConstIterator it = pMods_.begin(); it != pMods_.end();
       ++it) {
    (*it)->write(out);
  }
  out << me_ << "Relaxation modifications:" << std::endl;
  for (ModificationConstIterator it = rMods_.begin(); it != rMods_.end();
       ++it) {
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
