//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file KnapsackList.cpp
 * \brief Define base class KnapsackList.
 * \author Serdar Yildiz, Argonne National Laboratory
 */

#include "KnapsackList.h"
#include "Constraint.h"
#include "LinearFunction.h"

using namespace Minotaur;

KnapsackList::KnapsackList()
{
  // To be filled.
}

KnapsackList::KnapsackList(ProblemPtr p)
{
  // Pointer to the problem that Knapsack constraints are being searched.
  p_ = p;
  list_ = (ConstraintVectorPtr) new ConstraintVector(); 
  // This creates the list of Knapsack constraints.
  generateList();
}

KnapsackList::~KnapsackList()
{
  // To be filled.
}

// For now, we add the constraints where a_i > 0 forall i, a_i <= b forall i, 
// and b >= 0.
// Here, we do not check if there is any cover cut in the constraint.
// We identify if the constraint is a knapsack constraint with specified properties.
void KnapsackList::evalConstraint(ConstraintConstIterator itCons)
{
  bool add = true;
  // Check the values coefficients to be suitable for Gu, Nemhauser,
  // and Savelsbergh. 
  double ub = (*itCons)->getUb();
  // Check if the right hand side is positive.
  if (ub < 0) {
    add = false;
  }
  LinearFunctionPtr lf = (*itCons)->getLinearFunction();
  VariableGroupConstIterator it;
  VariableGroupConstIterator begin = lf->termsBegin();
  VariableGroupConstIterator end   = lf->termsEnd();
  for (it=begin; it!=end; ++it) {
    // Check if coefficients are positive.
    // Check if variables are less than b.
    if ((it->second <= 0) || (it->second > ub)) {      
      add = false;
      break;
    }
    // Check if variable is binary.
    if (it->first->getType() != Binary) {
    // If it is integer then check bounds. 
      if ((it->first->getType() == Integer) &&
          (it->first->getLb() == 0) &&
          (it->first->getUb() == 1)) {
        continue;
      }
      add=false;
      break;
    } 
  }
  
  // Decide if the constraint should be added
  if (add == true) {
    addConstraint(itCons);
  }
}

void KnapsackList::addConstraint(ConstraintConstIterator it)
{
  ConstraintPtr cons = (*it);
  list_->push_back(cons);
}

void KnapsackList::generateList()
{  
  // Iterator for first and last constraint is obtained.
  ConstraintConstIterator begin = p_->consBegin();
  ConstraintConstIterator end   = p_->consEnd();
 
  /** For each constraint, the type of function is checked if it is "Linear".
   * If it is, ConstConstraintPtr is added to the knapsack constraint list. 
   */
  FunctionType funType;
  numConsChecked_ = 0;
  for (ConstraintConstIterator it=begin; it!=end; ++it) {
    numConsChecked_ += 1;
    funType = (*it)->getFunctionType();
    if (funType == Linear) {
      evalConstraint(it);
    }    
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
