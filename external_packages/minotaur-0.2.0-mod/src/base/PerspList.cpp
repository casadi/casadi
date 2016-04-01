//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file PerspList.cpp
 * \brief Define base class that determines perspective constraints..
 * \author Serdar Yildiz, Argonne National Laboratory
 */

#include <iostream>
using std::endl;
using std::flush;

#include "PerspList.h"
#include "Function.h"
#include "LinearFunction.h"
#include "QuadraticFunction.h"
#include "NonlinearFunction.h"

# define DEBUG_LEVEL -1

using namespace Minotaur;

/*********************************************************************************/
// TODO Serdar, we have to add the functionality so that we should consider
// the case for complimentary binary, 4x4-5*(1-u)<=0.
// I guess that we will change the partial derivative of constraint function
// by negative of itself. May be other fundamental changes are needed as well.
/*********************************************************************************/

PerspList::PerspList()
{
  // To be filled.
}

PerspList::PerspList(RelaxationPtr rel, EnvPtr env/*, MultilinearTermsHandlerPtr mlh*/)
  : env_(env), rel_(rel)/*, mlh_(mlh)*/
{
  // Initialize perspective constraint list.
  list_ = (PerspConsPtr) new PerspCons();
  // Initialize statistics.
  stats_ = new PerspListStats();
  stats_->totalcons  = 0;
  stats_->totalpersp = 0;
  // Debug preparation.
  if (DEBUG_LEVEL >= 0) {
    outfile_ = "PerspListDebug.txt";
    // This is done just to clean the output debug file.
    output_.open(outfile_.c_str());
    output_.close();
  }
  // Generate the lists.
  generateList();
}

void PerspList::generateList()
{
  // Iterators for the first and last constraint.
  ConstraintConstIterator it;
  ConstraintConstIterator begin = rel_->consBegin();
  ConstraintConstIterator end   = rel_->consEnd();

  // Current constraint being checked.
  ConstConstraintPtr cons;
  // Shows if it is a perspective constraint.
  bool ispersp =  false;
  // Iterate through each constraint.
  for (it=begin; it!=end; ++it) {
    cons = *it;
    // Map for variable bound constraints.
    VarUbLbPtr boundcons = (VarUbLbPtr) new VarUbLb();
    // Binary variable 
    VariablePtr binvar; 
    // If binary variable is included 
    // Check if constraint perspective.
    ispersp = evalConstraint(cons, boundcons, binvar);
    if (ispersp) {
      // If it is a perpsective constraint, then add it to list.
      // Add constraint will be different
      addConstraint(cons, boundcons, binvar);
      if (DEBUG_LEVEL >= 9) {
        printPersp(cons, boundcons, binvar);
      }
      stats_->totalpersp += 1;
    }
    stats_->totalcons += 1;
  }
}

void PerspList::addConstraint(ConstConstraintPtr cons, VarUbLbPtr boundcons, VariablePtr binvar) 
{
  list_->insert(std::pair<ConsVar, VarUbLbPtr> (ConsVar(cons,binvar), boundcons) );
}


bool PerspList::evalConstraint(ConstConstraintPtr cons, VarUbLbPtr boundcons,
                               VariablePtr& binvar)
{
  // Type of function considered.
  FunctionType type;
  // Function type of constraint.
  type = cons->getFunctionType();
  // We do not consider linear constraints.
  if (type == Linear){
    return false;
  }

  // Get function of constraint.
  const FunctionPtr f = cons->getFunction(); 
  // Binary variable in constraint.
  //ConstVariablePtr binvar;
  // add one more parameter that stores the binary variable.
  bool vartypeok = checkVarTypes(f, binvar);
  // If all the variables are not continuos or at most one of them is binary
  // do not consider constraint for perspective cut generation.
  if (vartypeok ==  false) {
    return false;
  }
 
  // Check if we can consider the constraint further.
  // We have to check if the constraint is separable.
  bool isseparable = false;
  if (binvar == NULL) {
    isseparable = true;
  } else {
    isseparable = separable(cons, binvar);
  }
  if (isseparable == false) {
    return false;
  }
  
  // Shows if all variables are bounded by binary.
  bool boundsok = false;
  if (binvar == NULL){
    VarSetPtr binaries = (VarSetPtr) new VarSet();
    // Take the first element of constraint for initial binary search.
    ConstVariablePtr initvar = *(f->varsBegin());
    initialBinary(initvar, binaries);
    // If there is no binary for that then terminate cut generation.
    if (binaries->size() == 0) {
      return false;
    }
    // Iterate through all binaries.
    VarSetConstIterator it;
    VarSetConstIterator begin = binaries->begin();
    VarSetConstIterator end   = binaries->end();
    for (it=begin; it!=end; ++it) {
      binvar = *it;
      boundsok = checkVarsBounds(f, binvar, boundcons);
      if (boundsok == true) {
        return true;
      }
    }
  } else {
    // For each variable check if it is bounded by a binary variable.
    boundsok = checkVarsBounds(f, binvar, boundcons);
  }
  // If any variable is not bounded by binary variable, 
  // constraint is not a perspective constraint.
  if (boundsok ==  false) {
    return false;
  }
  
  // to be continued.  
  return true;
}

bool PerspList::separable(ConstConstraintPtr cons, ConstVariablePtr binvar)
{
  // Quadratic part should not include the binary variable u.
  QuadraticFunctionPtr qf = cons->getQuadraticFunction();
  if ((qf != NULL) && (qf->getNumVars() >= 1)) {
    if (qf->hasVar(binvar) == true) {
      return false;
    }
  }
  // Nonlinear part should not include the binary variable u.
  NonlinearFunctionPtr nlf = cons->getNonlinearFunction();
  if ((nlf != NULL) && (nlf->numVars() >= 1)) {
    if (nlf->hasVar(binvar) == true) {
      return false;
    }
  }    
  // If it comes to here,  
  return true;
}

bool PerspList::checkVarsBounds(const FunctionPtr f, ConstVariablePtr binvar, 
                                VarUbLbPtr boundcons)
{
  // Check if all variables are bounded by binary variable.
  // Iterators for variables.
  VarSetConstIterator it;
  VarSetConstIterator begin = f->varsBegin();
  VarSetConstIterator end   = f->varsEnd();
  // Current variable considered.
  ConstVariablePtr var;
  // For each variable check the bounds from these constraints.
  // We must not consider binvar in this loop, it should not be bounded.
  for (it=begin; (it!=end) && (*it!=binvar); ++it) {
    var = *it;
    // Construct map for variable bound constraints.
    bool varbounded = checkVarBounds(var, binvar, boundcons);
    // If variable is not bounded, then constraint is not a perspective constraint.
    if (varbounded == false) {
      return false;
    }
  }
  return true;
}

bool PerspList::checkVarBounds(ConstVariablePtr var, ConstVariablePtr binvar,
			       VarUbLbPtr boundcons)
{
  // Shows if variable is upper bounded.
  bool ubbounded = false;
  // Shows if variable is lower bounded.
  bool lbbounded = false;

  // First, check if variable is lower bounded or upper bounded by 0.
  double varlb = var->getLb();
  if (varlb == 0) {
    lbbounded = true;
  }
  double varub = var->getUb();
  if (varub == 0) {
    ubbounded = true;
  }

  // If any bound is assigned from variable bounds, 
  // corresponding bound constraint will be uninitialized as below.
  // Lb bounding constraint.
  ConstConstraintPtr lbcons = (ConstConstraintPtr) new Constraint();
  // Ub bounding constraint.
  ConstConstraintPtr ubcons = (ConstConstraintPtr) new Constraint();

  // Get list of constraints that includes current variable.
  ConstrSet::iterator it;
  ConstrSet::iterator begin = var->consBegin();;
  ConstrSet::iterator end   = var->consEnd();
  // Current constraint considered.
  ConstConstraintPtr cons;
  // Iterate through each constraint
  for (it=begin; it!=end; ++it) {
    cons = *it;
    // Function of constraint.
    const FunctionPtr f = cons->getFunction();
    // Type of constraint.
    FunctionType type = cons->getFunctionType();
    // Only consider linear constraints.
    if (type != Linear) {
      continue;
    } 
    // Get linear function.
    const LinearFunctionPtr lf = cons->getLinearFunction();
    
    // Number of variables in the constraint should be two
    // and one of them is current variable and the other one is binvar.
    UInt numvars = f->getNumVars();
    if (numvars != 2) {
      continue;
    }

    // Coefficient of variable.
    double coeffvar = lf->getWeight(var);
    // Coefficient of binary variable.
    double coeffbin = lf->getWeight(binvar);
    // Bounds of constraint.
    double lb = cons->getLb();
    double ub = cons->getUb();
    // Check upper and lower bounds.
    if (ub == 0) {
      if (lbbounded == false) {
        // If the following is correct, then it bounds from lb.
        if ( (coeffvar < 0) && (coeffbin > 0) ) {
          lbcons = cons;
          lbbounded = true;
        }
      }
      // Upper bound of constraint should be zero in order to be considered 
      // to bound variable from up.
      // Check if it bounds variable from up.      
      // If the following is correct, then it bounds from up.
      if (ubbounded == false) {
        if ( (coeffvar > 0) && (coeffbin < 0) ) {
          ubcons = cons;
          ubbounded = true;
        } 
      }
    } // end of if ub is zero.
    
    // Check upper and lower bounds.
    if (lb == 0) {
      // If the following is correct, then it bounds from lb.
      // Lower bound of the constraint should be zero.
      if (lbbounded == false) {
        if ( (coeffvar > 0) && (coeffbin < 0) ) {
          lbcons = cons;
          lbbounded = true;
        }
      }
      // If the following is correct, then it bounds from up.
      if (ubbounded == false) {
        if ( (coeffvar < 0) && (coeffbin > 0)) {
          ubcons = cons;
          ubbounded = true;
        }
      }
    } // end if lb is zero.   

    // If variable is both bounded from up and down from binary variable, then
    // we say it is bounded by binary variable.
    if ( (lbbounded == true) && (ubbounded == true) ) {
      std::pair<ConstConstraintPtr, ConstConstraintPtr> lbub(lbcons, ubcons);
      std::pair<ConstVariablePtr, std::pair<ConstConstraintPtr, ConstConstraintPtr> > varlbub (var, lbub);
      boundcons->insert(varlbub);
      return true;
    }

  } // end of for loop to consider all constraints that has the variable.

  // If it comes to here, variable is not bounded by binary variable.
  return false;
}


bool PerspList::checkVarTypes(const FunctionPtr f, ConstVariablePtr& binvar)
{
  // Check if all variables are binary or only one of them is binary.
  // Iterator for variables.
  VarSetConstIterator it;
  VarSetConstIterator begin = f->varsBegin();
  VarSetConstIterator end   = f->varsEnd();
  // Current variable considered.
  ConstVariablePtr var;
  // Type of variable considered.
  VariableType type;
  // Number of binary variables.
  UInt numbins = 0;
  // Number of integer binary variables.
  UInt numintbins = 0;
  // Iterate through all variables.
  for (it=begin; it!=end; ++it) {
    var = *it;
    type = var->getType();
    // Check the type of variable.
    switch (type) {
    case Binary:
      binvar = var;
      numbins += 1; 
      // If number of binary variables is more than one
      // we do not consider constraint for perspective cuts generation.
      if (numbins + numintbins >= 2) {
	return false;
      }
      break;
    case Integer:
      if ( (var->getLb() == 0) && (var->getUb() == 1) ) {
	binvar = var;
	numintbins += 1; 
      } else {
	// It is a general variable
	// we do not consider constraint.
	return false;
      }
      
      // If number of binary variables is more than one
      // we do not consider constraint for perspective cuts generation.
      if (numintbins + numbins >= 2) {
	return false;
      }
      break;
    case Continuous: /* Do nothing.*/ break;
    default:
      // If it comes to here, we have a variable which was not expected 
      // when algorithm was designed.
      return false;
    }
  }

  // If it comes here, it means all variables are continuous 
  // or there exists only one binary variable. 
  return true;
}

bool PerspList::initialBinary(ConstVariablePtr var, VarSetPtr binaries)
{
  // List of constraints that has variable considered.
  ConstrSet::const_iterator it;
  ConstrSet::const_iterator begin = var->consBegin();
  ConstrSet::const_iterator end   = var->consEnd();
  // Current constraint considered.
  ConstConstraintPtr cons;
  // Iterate through each constraint.
  for (it=begin; it!=end; ++it) {
    cons = *it;
    // Function of constraint.
    const FunctionPtr f = cons->getFunction();
    // Type of function.
    FunctionType type = cons->getFunctionType();
    // Only consider linear constraints.
    if (type != Linear) {
      continue;
    }
    // Get linear function.
    const LinearFunctionPtr lf = cons->getLinearFunction();
    // Number of variables should be two.
    UInt numvars = f->getNumVars();
    if (numvars != 2) {
      continue;
    }
    // Check if the other variable is binary.
    // Iterators for variables in constraint.
    VariableGroupConstIterator itvar;
    VariableGroupConstIterator beginvar = lf->termsBegin();
    VariableGroupConstIterator endvar   = lf->termsEnd();
    ConstVariablePtr curvar;
    VariableType vartype;
    double varlb;
    double varub;
    for (itvar=beginvar; itvar!=endvar; ++itvar) {
      curvar = itvar->first;
      vartype = curvar->getType();
      if (vartype == Binary) {
        binaries->insert(curvar);
      }
      if (vartype == Integer) {
        varlb = curvar->getLb();
        varub = curvar->getUb();
        if ( (varlb == 0) && (varub==1) ) {
          binaries->insert(curvar);
        }
      }
    } // end of for for variables.

  } // end of for loop.

  if (binaries->size() >= 1) {
    return true;
  } else {
    return false;
  }
}


PerspList::~PerspList()
{
  // Deallocate memory.
  // If statistics is created, deallocated the memory.
  if (stats_) {
    delete stats_;
  }
}

void PerspList::printPersp(ConstConstraintPtr cons, VarUbLbPtr boundcons, 
                           ConstVariablePtr binvar)
{
  output_.open(outfile_.c_str());
  output_ << endl;
  output_ << "Perspective constraint is: " << endl;
  cons->write(output_);
  output_ << "Binary variable for this constraint is: " << endl;
  binvar->write(output_);
  // Iterators for the constraints for each variable and its bounding
  // constraints.
  VarUbLb::const_iterator it;
  VarUbLb::const_iterator begin = boundcons->begin();
  VarUbLb::const_iterator end   = boundcons->end();
  // Print out for each variable and bounds
  ConstVariablePtr curvar;
  ConstConstraintPtr lbcons;
  ConstConstraintPtr ubcons;
  for (it=begin; it!=end; ++it) {
    curvar = it->first;
    output_ << "Bound constraints for variable " << curvar->getName()
            << " : " << endl;  
    lbcons = it->second.first;
    output_ << "Lower bounding constraint is: " << endl;
    lbcons->write(output_);
    ubcons = it->second.second;
    output_ << "Upper bounding constraint is: " << endl;
    ubcons->write(output_);
  }
  output_.close();
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
