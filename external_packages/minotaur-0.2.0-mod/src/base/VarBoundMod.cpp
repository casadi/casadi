// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file VarBoundMod.cpp
 * \brief Implement the Modification class VarBoundMod, that is used to store
 * modifications to a bound on a variable.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <iostream>

#include "MinotaurConfig.h"
#include "Engine.h"
#include "Problem.h"
#include "Relaxation.h"
#include "VarBoundMod.h"
#include "Variable.h"


using namespace Minotaur;


VarBoundMod::VarBoundMod(VariablePtr var, BoundType lu, double new_val)
  : lu_(lu),
    newVal_(new_val),
    var_(var)
{
  switch (lu) {
   case (Lower):
     oldVal_ = var->getLb();
     break;
   case (Upper):
     oldVal_ = var->getUb();
     break;
   default:
     break;
  }
  //std::cout << "varboundmod: old value = " << oldVal_ << std::endl;
}


VarBoundMod::~VarBoundMod()
{
  var_.reset();
}


ModificationPtr VarBoundMod::fromRel(RelaxationPtr rel, ProblemPtr) const
{
  VarBoundModPtr mod = (VarBoundModPtr) new VarBoundMod(
                                        rel->getOriginalVar(var_),
                                        lu_, newVal_);
  mod->oldVal_ = oldVal_;
  return mod;
}


VariablePtr VarBoundMod::getVar() const
{
  return var_;
}


BoundType VarBoundMod::getLU() const
{
  return lu_;
}


double VarBoundMod::getNewVal() const
{
  return newVal_;
}


void VarBoundMod::applyToProblem(ProblemPtr problem) 
{
  problem->changeBound(var_, lu_, newVal_);
}


ModificationPtr VarBoundMod::toRel(ProblemPtr, RelaxationPtr rel) const
{
  VarBoundModPtr mod = (VarBoundModPtr) new VarBoundMod(
                                        rel->getRelaxationVar(var_),
                                        lu_, newVal_);
  mod->oldVal_ = oldVal_;
  return mod;
}


void VarBoundMod::undoToProblem(ProblemPtr problem) 
{
  problem->changeBound(var_, lu_, oldVal_);
}


void VarBoundMod::write(std::ostream &out) const
{
  out << "var bound mod: "
      << "var name = " << var_->getName()
      << " bound type = ";
  if (lu_==Lower) {
    out << "lb"; 
  } else {
    out << "ub";
  }
  out << " old value = " << oldVal_
      << " new value = " << newVal_
      << std::endl;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

VarBoundMod2::VarBoundMod2(VariablePtr var, double new_lb, double new_ub) 
  : newLb_(new_lb),
    newUb_(new_ub),
    var_(var)
{
  oldLb_ = var->getLb();
  oldUb_ = var->getUb();
}


VarBoundMod2::~VarBoundMod2()
{
  var_.reset();
}


ModificationPtr VarBoundMod2::fromRel(RelaxationPtr rel, ProblemPtr) const
{
  VarBoundMod2Ptr mod = (VarBoundMod2Ptr) new VarBoundMod2(
                                          rel->getOriginalVar(var_),
                                          newLb_, newUb_);
  mod->oldLb_ = oldLb_;
  mod->oldUb_ = oldUb_;

  return mod;
}


VariablePtr VarBoundMod2::getVar() const
{
  return var_;
}


double VarBoundMod2::getNewLb() const
{
  return newLb_;
}


double VarBoundMod2::getNewUb() const
{
  return newUb_;
}


void VarBoundMod2::applyToProblem(ProblemPtr problem)
{
  problem->changeBound(var_, newLb_, newUb_);
}


ModificationPtr VarBoundMod2::toRel(ProblemPtr, RelaxationPtr rel) const
{
  VarBoundMod2Ptr mod = (VarBoundMod2Ptr) new VarBoundMod2(
                                          rel->getRelaxationVar(var_),
                                          newLb_, newUb_);
  mod->oldLb_ = oldLb_;
  mod->oldUb_ = oldUb_;

  return mod;
}


void VarBoundMod2::undoToProblem(ProblemPtr problem)
{
  problem->changeBound(var_, oldLb_, oldUb_);
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
