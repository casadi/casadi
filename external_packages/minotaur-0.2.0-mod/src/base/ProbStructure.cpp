//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file ProbStructure.cpp
 * \brief Define base class for Problem structure.
 * \author Serdar Yildiz, Argonne National Laboratory
 */
#include "ProbStructure.h"
#include "LinearFunction.h"

using namespace Minotaur;

ProbStructure::ProbStructure()
{
  // To be filled.
}

ProbStructure::ProbStructure(ProblemPtr p, EnvPtr env)
  : env_(env), p_(p)
{
  // Initialize statistics.
  stats_ = new ProbStructStats();
  stats_->totalcons = 0;
  stats_->totalGUBs = 0; 
  // Create an empty list for GUBs.
  list_    = (ConstConstraintVectorPtr) new ConstConstraintVector();
  // Create an empty list for GUBs corresponding to variables.
  varlist_ = (VarConsPtr) new VarCons();
  // Generate the lists.
  generateLists();
}

ProbStructure::~ProbStructure()
{
  // Deallocate memory.
  // If statistics is created, deallocate the memory.
  if (stats_) {
    delete stats_;
  }
  
}

void ProbStructure::generateLists()
{
  // Iterators for the first and last constraint.
  ConstraintConstIterator it;
  ConstraintConstIterator begin = p_->consBegin();
  ConstraintConstIterator end   = p_->consEnd();
  
  // Current constraint being checked.
  ConstConstraintPtr cons;
  // Shows if the constraint is GUB.
  bool isGub = false;
  // Iterate through each constraint.
  for (it=begin; it!=end; ++it) {
    cons = *it;
    // Check if constraint is GUB.
    isGub = evalConstraint(cons);
    if (isGub) {
      // If it is a GUB, then add it to the GUB lists of corresponding variables.
      addConstraint(cons);
      stats_->totalGUBs += 1;
    }
    stats_->totalcons += 1;
  }
}

/** Assumption is that:
 * All the GUB functions in the form of ax <= b where
 * b != 0, b > 0 and a/b = 1;
 * We do not consider a constraint if b < 0.
 * May be we should do the same thing we do for b>0 case here. 
 * Furthermore, variables can be binary or integer with bounds [0,1].
 */
bool ProbStructure::evalConstraint(ConstConstraintPtr cons)
{
  // The type of the function considered.
  FunctionType type;
  // Function type of constraint.
  type = cons->getFunctionType();
  // If function is linear then evaluate.
  if (type == Linear) {
    // Get upper bound of constraint.
    double b = cons->getUb();
    // If right hand side is negative, then it is not GUB.
    if (b < 0) { // Serdar add tolerance.
      // Do something for that
      return false;
    } else  if (b == 0) { // Serdar add a tolerance here.
      // If rhs is 0, then it is not a GUB.
      return false;
    } else {
      // Get linear function.
      LinearFunctionPtr lf = cons->getLinearFunction();
      // Iterators for variables.
      VariableGroupConstIterator it;
      VariableGroupConstIterator begin = lf->termsBegin();
      VariableGroupConstIterator end   = lf->termsEnd();
      // Current variable.
      ConstVariablePtr var;
      // Coefficient of current variable.
      double coeff = 0.0;
      // Ratio of a_i/b.
      double ratio = 0.0; 
      for (it=begin; it!=end; ++it) {
        var = it->first;
        // If coefficient is not binary, then it is not GUB.
        if (var->getType() != Binary) {
          // Check if it is integer.
          if (var->getType() == Integer) {
            // Check the bounds
            if ((var->getLb() != 0) || (var->getUb() != 1)) {
              return false;
            } // Variable is binary
          } else {
            // Variable is not integer.
            return false;
          } // Variable is not binary
        } // Variable is binary.
        coeff = it->second;
        // If the ratio is a/b!=1, then it is not GUB.
        ratio = coeff / b;
        if (ratio != 1) {
          return false;
        }
      } // for loop for variables ends
    } // else for b > 0  ends.
  } else {
    // Function is not linear so not GUB.
    return false;
  }

  // Constraint is GUB.
  return true;
}

void ProbStructure::addConstraint(ConstConstraintPtr cons)
{
  // Add the constraint to general list.
  list_->push_back(cons);
  // Linear function of the constraint.
  LinearFunctionPtr lf = cons->getLinearFunction();
  // Iterators for the variables in the constraint.
  VariableGroupConstIterator it;
  VariableGroupConstIterator begin = lf->termsBegin();
  VariableGroupConstIterator end   = lf->termsEnd();
  // Current variable.
  ConstVariablePtr var;
  // Map pair for the variable.
  VarConsIterator varcons;
  // GUB list of the variable.
  ConstConstraintVectorPtr constraints;
  // Iterate through all the variables.
  for (it=begin; it!=end; ++it) {
    // Add the constraint for the corresponding GUB lists of variables.
    var = it->first;
    // If does not exist in the map already, insert it.
    varcons = varlist_->find(var);
    if (varcons == varlist_->end()) {
      // Construct constraint list for new variable.
      constraints = (ConstConstraintVectorPtr) new ConstConstraintVector();
      // Add constraint to constraint list.
      constraints->push_back(cons);
      // Create a pair<var,constraints> to add to map.
      VarConsPair newvar(var,constraints);
      varlist_->insert(newvar);
    } else {
      // Variable is already in the map. Just add the constraint to its list.
      constraints = varcons->second;
      constraints->push_back(cons);
    }
  }
}


UInt ProbStructure::getNumVarGUB(ConstVariablePtr var) const
{
  // If variable is not in the map, then there is no GUB that includes the variable.
  VarConsIterator varcons = varlist_->find(var);
  ConstConstraintVectorPtr constraints;
  UInt numcons;
  if (varcons == varlist_->end()) {
    return 0;
  } else {
    constraints = varcons->second;
    numcons = constraints->size();
    return numcons;
  }
}

ConstConstraintVectorPtr ProbStructure::getVarGUBs(ConstVariablePtr var) const
{
  // If a variable is not in the map, then return an empty list.
  VarConsIterator varcons = varlist_->find(var);
  ConstConstraintVectorPtr constraints = 
    (ConstConstraintVectorPtr) new ConstConstraintVector();
  if (varcons == varlist_->end()) {
    return constraints;
  } else {
    constraints = varcons->second;
    return constraints;
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
