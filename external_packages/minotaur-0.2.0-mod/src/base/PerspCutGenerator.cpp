//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//


/**
 * \file PerspCutGenerator.cpp 
 * \brief Declare PerspCutGenerator class. 
 * \author Serdar Yildiz, Argonne National Laboratory 
*/

#include <limits>
#include <algorithm>
#include <cmath>
#include <iostream>
using std::endl;
using std::flush;
#include <sstream>
using std::ostringstream;


#include "PerspList.h"
#include "PerspCutGenerator.h"
#include "Function.h"
#include "LinearFunction.h"
#include "Operations.h"
#include "Option.h"

# define DEBUG_LEVEL -1

using namespace Minotaur;

PerspCutGenerator::PerspCutGenerator() {}

PerspCutGenerator::PerspCutGenerator(RelaxationPtr rel, ConstSolutionPtr sol,
				     EnvPtr env)
  : env_(env), rel_(rel), s_(sol)
{
  initialize();
  generateAllCuts();
}

bool PerspCutGenerator::generateAllCuts()
{
  if (persplist_->getNumPersp() >= 1) {
    const double * x = s_->getPrimal();
    ConstPerspConsPtr perspcons = persplist_->getPerspCons();
    // Iterators for perpsepective constraints.
    PerspCons::const_iterator it;
    PerspCons::const_iterator begin = perspcons->begin();
    PerspCons::const_iterator end   = perspcons->end();
    for (it=begin; it!=end; ++it) {
      stats_->perspcons += 1;
      ConstConstraintPtr cons = it->first.first;
      ConstVariablePtr binvar = it->first.second;
      UInt binindex = binvar->getIndex();
      double binsol = x[binindex];
      // Check if binary variable has an integral value.
      bool binint = IsInt(binsol, intTol_);
      if (binint) {
        stats_->perspnotcons += 1;
        continue;
      }
      // Print out constraint for debugging.
      if (DEBUG_LEVEL >= 9) {
        output_.open(outfile_.c_str(), std::ios_base::app);
        output_ << "\n\n Constraint considered is: " << endl;
        cons->write(output_);
        output_.close();
      }
      stats_->perspconsidered += 1;
      generateCut(cons, binvar);
    } // end of for loop.
  } else {
    if (DEBUG_LEVEL >= 9) {
      output_ << "PerspCutGenerator::generateAllCuts. No constraint used for cut generation." << endl; 
    }
  } // end of else.
  
  return true;
}

bool PerspCutGenerator::generateCut(ConstConstraintPtr cons,ConstVariablePtr binvar)
{
  FunctionPtr f = cons->getFunction();
  LinearFunctionPtr lf = cons->getLinearFunction();
  double coeffbin = lf->getWeight(binvar);
  // Number of variables in problem.
  UInt numvars = rel_->getNumVars();
  // Get primal solution.
  const double * x = s_->getPrimal();
  // Get index of binary variable.
  UInt binindex = binvar->getIndex();
  // Generate the solution (p/u,u)
  double * y = new double[numvars];
  std::fill(y, y + numvars, 0);
  UInt indexvar = 0;
  ConstVariablePtr curvar;
  double solbin = x[binindex];
  // If binary variable slution is too small, or if the value is very close to
  // one, then do not proceed.
  double tolerance = 1e-3;
  if (solbin <= tolerance || solbin >=(1-tolerance)) {
    // clean up.
    delete [] y;
    return false;
  }
  double varcursol = 0.0;
  for (VarSetConstIterator it=f->varsBegin(); it!=f->varsEnd(); ++it) {
    curvar = *it;
    indexvar = curvar->getIndex();
    if (indexvar != binindex) {
      varcursol = x[indexvar];
      y[indexvar] = varcursol/solbin; 
    } else {
      y[binindex] = solbin;
    }
  } // end of for loop.
  // In case binary variable is not included in the constraint,
  // we save its solution value here, 
  y[binindex] = solbin;
  // Evaluate function at given soution (p*/u*, u*).
  int error = 0;
  double conseval = f->eval(y, &error);
  // Evaluate gradient at given solution (p*/u*, u*).
  int errorgr = 0;
  double * consgradient = new double[numvars];
  std::fill(consgradient, consgradient + numvars, 0);
  f->evalGradient(y, consgradient, &errorgr);

  // Iterators for variables.
  VarSetConstIterator it;
  VarSetConstIterator begin = f->varsBegin();
  VarSetConstIterator end   = f->varsEnd();
  // Linear function of perspective cut.
  LinearFunctionPtr lfpersp = (LinearFunctionPtr) new LinearFunction();
  // Constant term of perspective cut.
  double feval = conseval - coeffbin*solbin;
  double constpersp = solbin*feval +coeffbin*solbin;
  // Gradient of binary variable.
  double gradu = feval + coeffbin;
  // Gradient of continuous variables in the constraint.
  double gradfi = 0.0;
  // Index of current variable considered.
  UInt varindex = 0;
  // Soluion of variable at vector y.
  double ysolvar = 0.0;
  for (it=begin; it!=end; ++it) {
    ConstVariablePtr var = *it; 
    varindex = var->getIndex();
    gradfi =consgradient[varindex]; 
    // Add terms for each variable.
    lfpersp->addTerm(var, gradfi);
    // Solution value of variable.
    ysolvar = y[varindex];
    // Increment constant of constraint.      
    constpersp -= gradfi * ysolvar;
    // gradient for variable u is calculated.
    gradu -= gradfi * ysolvar; 
  }
  constpersp -= gradu*solbin; 
  // We can add binary variable to linearization, when we considered all
  // variables except the binary variable so that the linearization
  // coefficient of binary variable is calculated.
  lfpersp->addTerm(binvar, gradu);
  
  FunctionPtr fpersp = (FunctionPtr) new Function(lfpersp);
  double infin = std::numeric_limits<double>::infinity();
  CutPtr cut = (CutPtr) new Cut(rel_, fpersp, -infin, -constpersp, false, false);
  
  // add cut to the lists.
  addCut(cut);
  
  // clean up
  delete [] y;
  delete [] consgradient;

  // Return cut from here.
  return true;
}

bool PerspCutGenerator::addCut(CutPtr cut)
{
  // Total number of cuts increased by one.
  stats_->totalcuts += 1;

  // In debug mode, we write all cuts generated.
  if (DEBUG_LEVEL >= 9) {
    // These lines converts number to string.
    int numcuts =  stats_->totalcuts;
    string numc;
    ostringstream convert;
    convert << numcuts;
    numc = convert.str();
    // until here.
    output_.open(outfile_.c_str());
    output_ << "Cut " << numc << ":" << endl;
    cut->write(output_);
    output_.close();
  }

  bool cutexists = checkExists(cut);
  if(cutexists == false) {
    stats_->cuts += 1;
    // Add to the cut list.
    cutList_.push_back(cut);
    double viol = violation(cut);
    if (viol > objtol_) {
      viollist_.push_back(cut);
      viols_.push_back(viol);
      stats_->violated += 1;
    }
  }
  return false;
}

double PerspCutGenerator::violation(CutPtr cut)
{
  FunctionPtr f= cut->getFunction();
  const double * x = s_->getPrimal();
  int error = 0;
  double evaluation = f->eval(x, &error);
  double violub = std::max(0.0, evaluation - cut->getUb());
  
  return violub;
}

bool PerspCutGenerator::checkExists(CutPtr cut)
{
  double rhs = cut->getUb();
  UInt numvars = rel_->getNumVars();
  std::vector<double> coeffs(numvars, 0);
  UInt varindex = 0;
  ConstVariablePtr curvar;
  double coeff = 0.0;
  double dividedcoeff;
  FunctionPtr f = cut->getFunction();
  LinearFunctionPtr lf = f->getLinearFunction();
  // Iterators for variables
  VariableGroupConstIterator it;
  VariableGroupConstIterator begin = lf->termsBegin();
  VariableGroupConstIterator end   = lf->termsEnd();
  for (it=begin; it!=end; ++it) {
    curvar = it->first;
    varindex = curvar->getIndex();
    coeff = it->second;
    if (rhs >= eTol_) {
      dividedcoeff = coeff / double(rhs);
    } else {
      dividedcoeff = coeff;
    }
    coeffs[varindex] = dividedcoeff;
  }
    
  // Check if the cut already exists.
  std::map< std::vector<double>, UInt >::const_iterator found = cutmap_.find(coeffs);
  std::map< std::vector<double>, UInt >::const_iterator endmap = cutmap_.end();
  if (found == endmap) {
    cutmap_.insert(std::pair< std::vector<double>, UInt >(coeffs,stats_->totalcuts));
  return false; 
} else
    return true;
}

PerspCutGenerator::~PerspCutGenerator()
{
  // deallocate heap.
  if (stats_) {
    delete stats_;
  }
}

void PerspCutGenerator::initialize()
{
  stats_ = new PerspGenStats();
  stats_->totalcuts = 0;
  stats_->cuts = 0;
  stats_->violated = 0;
  stats_->noviol = 0;
  stats_->perspcons = 0;
  stats_->perspnotcons = 0;
  stats_->perspconsidered = 0;
  stats_->time = 0.0;
  // Perspective constraints are identified.
  persplist_ = (PerspListPtr) new PerspList(rel_, env_); 
  objtol_ = 1e-6;
  intTol_ = env_->getOptions()->findDouble("int_tol")->getValue();
  eTol_ = 1e-6; 
  if (DEBUG_LEVEL >= 0) {
    outfile_ = "PerspCutGeneratorDebug.txt";
    output_.open(outfile_.c_str());
    // Print out current solution.
    VariableConstIterator it;
    VariableConstIterator begin = rel_->varsBegin();
    VariableConstIterator end   = rel_->varsEnd();
    VarVector solution;
    ConstVariablePtr curvar;
    for (it=begin; it!=end; ++it) {
      curvar = *it;
      solution.push_back(curvar);
    }// end of for loop.
    output_.precision(2);
    output_ << "Given primal solution " << endl;
    s_->writePrimal(output_, &solution);
    output_ << "\nDual solution is written " << endl;
    s_->writeDual(output_);
    output_.close();
  }// end of debug.
}

bool PerspCutGenerator::checkIntegral(RelaxationPtr p, ConstSolutionPtr s)
{
  // Iterators for variables in problem.
  VariableConstIterator it;
  VariableConstIterator begin = p->varsBegin();
  VariableConstIterator end   = p->varsEnd();
  
  // Primal solution.
  const double * x = s->getPrimal();
  // Current variable.
  ConstVariablePtr var;
  // Index of variable.
  UInt index = 0;
  // Value of variable.
  double value = 0.0;
  // Absolute of fractional part of variable value.
  double fraction = 0.0;
  // Type of variable.
  VariableType type;
  // Iterate through all variables in the problem.
  for (it=begin; it!=end; ++it) {
    var = (*it);
    type = var->getType();
    // Check if variable is type integer.
    if (type == Binary || type == Integer) {
      index = var->getIndex();
      value = x[index];
      fraction = fabs(value - floor(value+0.5));
      if (fraction > intTol_) {
        // Check if the variable has fractional value.
        return false;
      }
    }
  } // end of for loop.
  
  // Solution is integral.
  return true;
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
