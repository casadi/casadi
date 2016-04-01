// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

/**
 * \file CxUnivarHandler.cpp
 * \brief Implement the handler for functions of the general univariate convex
 * form y=g(x).
 * \author Jim Luedtke, Jeff Linderoth UW-Madison 
 */


#include <cmath>
#include <iostream>
#include <iomanip>

#include "MinotaurConfig.h"
#include "Branch.h"
#include "BrVarCand.h"
#include "Constraint.h"
#include "CxUnivarHandler.h"
#include "Environment.h"
#include "Function.h"
#include "LinMods.h"
#include "Logger.h"
#include "Node.h"
#include "Objective.h"
#include "Operations.h"
#include "Option.h"
#include "QuadraticFunction.h"
#include "ProblemSize.h"
#include "Relaxation.h"
#include "SolutionPool.h"
#include "Variable.h"

//#define SPEW 1

#undef DEBUG_CXUNIVARHANDLER

using namespace Minotaur;

const std::string CxUnivarHandler::me_ = "CxUnivarHandler: ";

CxUnivarConstraintData::CxUnivarConstraintData(double eTol,
					       double vTol,
					       ConstraintPtr newcon,
                                               ConstVariablePtr ivar,
                                               ConstVariablePtr ovar,
                                               char sense)
: eTol_(eTol), vTol_(vTol), con_(newcon),
  iv_(ivar),
  riv_(VariablePtr()),
  ov_(ovar),
  rov_(VariablePtr()),
  sense_(sense),
  secCon_(ConstraintPtr()),
  linCons_(ConstraintVector())
{

}


/// Creates initial relaxations   
void CxUnivarConstraintData::initRelax(RelaxationPtr rel, DoubleVector& tmpX,
                                       DoubleVector& grad)
{
  // get and set the relaxation relaxation variables that match the original
  // problem variables
  riv_ = rel->getVariable(iv_->getIndex());
  rov_ = rel->getVariable(ov_->getIndex());

  // This isn't used
  ModVector mods;

  // Add secant
  if (sense_ == 'E' || sense_ == 'L') {
    addSecant(rel, riv_, rov_, con_->getFunction(), tmpX, true, mods);
  }

  // Add linearizations 
  // TODO: Make strategy a parameter
  if (sense_ == 'E' || sense_ == 'G') {
    addLin(rel, iv_, ov_, con_->getFunction(), tmpX, grad, true, mods); 
  }

}

/// Update the current relaxation based on current variable bounds
void CxUnivarConstraintData::updateRelax(RelaxationPtr  rel , DoubleVector&
                                         tmpX , DoubleVector& grad,
					 ModVector &mods) 
{
  // Add secant
  if (sense_ == 'E' || sense_ == 'L') {
    addSecant(rel, riv_, rov_, con_->getFunction(), tmpX, false, mods);
  }

  // Add linearizations 
  // TODO: Make strategy a parameter
  if (sense_ == 'E' || sense_ == 'G') {
    addLin(rel, riv_, rov_, con_->getFunction(), tmpX, grad, false, mods); 
  }
}

double CxUnivarConstraintData::getViol(const DoubleVector &x)
{
  int error;
  double fval = con_->getFunction()->eval(x,&error); 

  // TODO: Put in a better (scaled) feasibility check here
  double absViol = 0.0;
  double relViol = 0.0;
  if (fval < con_->getLb() - eTol_) {
    absViol = con_->getLb() - fval;
  }
  if (fval > con_->getUb() + eTol_ && fval - con_->getUb() > absViol) {
    absViol = fval - con_->getUb(); 
  }
  relViol = absViol;
  if (fabs(fval) + absViol > 1.0) {
    relViol = absViol/(fabs(fval) + absViol);
  }

  return relViol;
}

bool CxUnivarConstraintData::isFeasible(const double* x)
{
  bool isfeas = true;
  int error;
  double fval = con_->getFunction()->eval(x,&error); 
  // TODO: Put in a better (scaled) feasibility check here
  if ((fval < con_->getLb() - eTol_ || fval > con_->getUb() + eTol_) && 
      riv_->getUb() - riv_->getLb() > vTol_) {
    isfeas = false;
  }
  return isfeas;
}


CxUnivarHandler::CxUnivarHandler(EnvPtr , ProblemPtr problem)
: eTol_(1e-6),
  vTol_(1e-5)
{
  problem_ = problem; 
  //logger_  = (LoggerPtr) new Logger((LogLevel) 
  //    env->getOptions()->findInt("handler_log_level")->getValue());
  logger_  = (LoggerPtr) new Logger((LogLevel) 0 );

  tmpX_.resize(problem->getNumVars(), 0.0);
  grad_.resize(problem->getNumVars(), 0.0);

}


CxUnivarHandler::~CxUnivarHandler()
{

}

void CxUnivarHandler::relaxInitInc(RelaxationPtr rel, bool *) 
{

  if (tmpX_.size() != problem_->getNumVars()) {
    tmpX_.resize(problem_->getNumVars(), 0.0);
  }
  if (grad_.size() != problem_->getNumVars()) {
    grad_.resize(problem_->getNumVars(), 0.0);
  }
    
  CxUnivarConstraintIterator dit; 
  
  for (dit = cons_data_.begin(); dit != cons_data_.end(); ++dit) {
    (*dit)->initRelax(rel, tmpX_, grad_);
  }
}

void CxUnivarHandler::relaxNodeInc(NodePtr , RelaxationPtr rel, bool *is_inf)
{

#if defined(DEBUG_CXUNIVARHANDLER)
  std::cout << "CxUnivarHandler::relaxNodeInc.  Current relaxation: "
            << std::endl;
  rel->write(std::cout);
  std::cout << "Current node: " << std::endl;
  node->write(std::cout);
#endif

  // You must apply all modifications to the nodes

  // I think this is all done by branches now.
  //  FALSE -- If you don't branch on a variable, you still need to update the relaxation from
  //  this node

  ModVector mods;
  CxUnivarConstraintIterator dit;   
  for (dit = cons_data_.begin(); dit != cons_data_.end(); ++dit) {
    (*dit)->updateRelax(rel, tmpX_, grad_, mods);
  }

  ModificationConstIterator it;
  for (it = mods.begin(); it != mods.end(); ++it) {
    assert(!"add Mod correctly here.");
    // node->addPMod(*it);
  }

 *is_inf = false;
}

void CxUnivarConstraintData::addLin(RelaxationPtr rel, ConstVariablePtr riv,
				    ConstVariablePtr rov, FunctionPtr fn,
                                    DoubleVector& tmpX, DoubleVector& grad,
                                    bool init, ModVector &mods)
{
 
  int error;
  ConstraintPtr cons; 
  double xlb = riv->getLb();
  double xub = riv->getUb();
  double fxlbval=0, fxubval=0, dfxlbval=0, dfxubval=0;
  double tmpxval, fxval, dfxval; 
  LinearFunctionPtr lf; 
  FunctionPtr f;

  // More sophisticated strategies hopefully could be obtained by simply
  // changing this array 
  int npts = 3;
  double xvals[] = {xlb, xub, (xub-xlb)/2.0};

#if defined(DEBUG_CXUNIVARHANDLER)
  std::cout << "Adding linearizations.  rix id: " << riv->getId() 
	    << " rix index: " << riv->getIndex() << " rov id: " << rov->getId() 
	    << " rov index: " << rov->getIndex()
	    << " xlb: " << xlb << " xub: " << xub << std::endl;
#endif
  
  for (int i = 0; i < npts; i++) {

    // Zero out tmpX and grad each time, or else bad things happen
    for (UInt j = 0; j < tmpX.size(); ++j) {
      tmpX[j] = 0.0;
      grad[j] = 0.0;
    }
    
    if (i == 2) {
      // Third linearization point taken to be where first two intersect:
      // x3 = (f'(xub)*xub - f'(xlb)*xlb + f(xlb) - f(xub))/(f'(xub) - f'(xlb))
      // Unless this would put it too close to one of the end points
      if (dfxubval - dfxlbval > 0.0001 || dfxubval - dfxlbval < -0.0001) {
        tmpxval = (dfxubval*xub - dfxlbval*xlb + fxlbval - fxubval)/
                  (dfxubval - dfxlbval);
        if (tmpxval < xlb + (xub-xlb)*0.05) {
          xvals[2] = xlb + (xub-xlb)*0.05;
        }
        else if (tmpxval > xub - (xub-xlb)*0.05) {
          xvals[2] = xub - (xub-xlb)*0.05;
        }
        else {
          xvals[2] = tmpxval;
        }
      }
    }
    tmpX[riv->getIndex()] = xvals[i];
    error = 0;
    fxval =  fn->eval(tmpX, &error);
    fn->evalGradient(&tmpX[0], &grad[0], &error);
#if defined(DEBUG_CXUNIVARHANDLER2)
    for (UInt j = 0; j < tmpX.size(); ++j) {
      std::cout << "x[" << j << "] = " << tmpX[j] << " dfdx[" << j << "] = "
                << grad[j] << std::endl;
    }
#endif
    dfxval = grad[riv->getIndex()];
    if (i == 0) {
       fxlbval = fxval;
       dfxlbval = dfxval; 
    }
    else if (i == 1) {
       fxubval = fxval;
       dfxubval = dfxval; 
    }
    // linearization:  rov >= f(xval) + f'(xval)(riv - xval) 
    //                 rov - f'(xval)*riv >= f(xval) - f'(xval)*xval
    lf = (LinearFunctionPtr) new LinearFunction();
    lf->addTerm(rov, 1.0);
    lf->addTerm(riv, -dfxval);
    if (init) {
        f = (FunctionPtr) new Function(lf);
    	cons = rel->newConstraint(f, fxval - dfxval*xvals[i], INFINITY);
	linCons_.push_back(cons);
    }
    else {
#if defined(DEBUG_CXUNIVARHANDLER)
       std::cout << "Will change 'linearization  ' constraint to have "
                 << "linear function: ";
       lf->write(std::cout);  
       std::cout << std::endl;
#endif

       rel->changeConstraint(linCons_[i], lf, fxval - dfxval*xvals[i], INFINITY); 
       LinConModPtr lcmod = (LinConModPtr) new LinConMod(linCons_[i], lf, 
                                                         fxval -
                                                         dfxval*xvals[i],
                                                         INFINITY); 
       mods.push_back(lcmod);
    }
  }
  tmpX[riv->getIndex()] = 0.0;
  grad[riv->getIndex()] = 0.0;

}


void CxUnivarConstraintData::addSecant(RelaxationPtr rel,
                                       ConstVariablePtr riv,
                                       ConstVariablePtr rov,
                                       FunctionPtr fn, DoubleVector& tmpX,
                                       bool init, ModVector &mods) 
{

  int error;
  double xlb, xub, fxlb, fxub, m, intercept;
  LinearFunctionPtr lf; 
  FunctionPtr f;

  // First add the secant inequalities based on variable bounds
  xlb = riv->getLb();
  xub = riv->getUb();
	
#if defined(DEBUG_CXUNIVARHANDLER)
  std::cout << "Adding secant on variable rix index: " << riv->getIndex() 
	    << " rov index: " << rov->getIndex()
	    << " xlb: " << xlb << " xub: " << xub << std::endl;
#endif 
  // no secant if unbounded either way
  if (xlb <= -0.9*INFINITY || xub >= 0.9*INFINITY) {
    std::cout << "Cannot add secant -- bound is infinite" << std::endl;
    return;
  }

  // TODO: Check the error value!
  tmpX[riv->getIndex()] = xlb;
  fxlb =  fn->eval(tmpX, &error);
  tmpX[riv->getIndex()] = xub;
  fxub =  fn->eval(tmpX, &error);
  tmpX[riv->getIndex()] = 0.0;

  // TODO: check/remedy numerical issues in this division
  if (xub - xlb > 10e-7) {
    m = (fxub - fxlb)/(xub - xlb);
  }
  else {
    m = 0.0;
  }

  intercept = fxlb - m*xlb;
  lf = (LinearFunctionPtr) new LinearFunction();
  lf->addTerm(rov, 1.0);
  lf->addTerm(riv, -m);

  // rovar <= m*rivar + intercept 
  if (init) {
     f = (FunctionPtr) new Function(lf);
     secCon_ = rel->newConstraint(f, -INFINITY, intercept);
  }
  else {
    rel->changeConstraint(secCon_, lf, -INFINITY, intercept);
    LinConModPtr lcmod = (LinConModPtr) new LinConMod(secCon_, lf, -INFINITY,
                                                      intercept);
    mods.push_back(lcmod);
  }

  
}

 
  
void CxUnivarHandler::addConstraint(ConstraintPtr newcon, ConstVariablePtr ivar,
                                    ConstVariablePtr ovar, char sense)
{

  Handler::addConstraint(newcon);
  cons_data_.push_back(CxUnivarConstraintDataPtr(new CxUnivarConstraintData
                                                 (eTol_, vTol_, newcon, ivar,
                                                  ovar, sense)));
}


bool CxUnivarHandler::isFeasible(ConstSolutionPtr sol, RelaxationPtr , bool &,
                                 double &)
{
  bool isfeas = true;
  CxUnivarConstraintIterator dit; 
  
  for (dit = cons_data_.begin(); dit != cons_data_.end(); ++dit) {
      if (!(*dit)->isFeasible(sol->getPrimal())) {
         isfeas = false;
         break;
      }	
  }

  return isfeas;

}


// Eventually, could add additional linearization inequalities for the convex
// side here... but not absolutely necessary
void CxUnivarHandler::separate(ConstSolutionPtr , NodePtr , RelaxationPtr ,
                               CutManager *, SolutionPoolPtr , bool *,
                               SeparationStatus *)
{
}


void CxUnivarHandler::getBranchingCandidates(RelaxationPtr,
                                             const DoubleVector &x,
                                             ModVector &, BrVarCandSet &cands,
                                             BrCandVector &, bool &is_inf)
{
  
  is_inf = true;
  CxUnivarConstraintIterator dit; 
  // Create a map of variables to their weights
  // Weights will be sum of scaled violation of constraints they are argument for
  std::map<ConstVariablePtr, double> allCands;
  std::map<ConstVariablePtr, double>::iterator curc_it;
  double curviol = 0.0;
  
  for (dit = cons_data_.begin(); dit != cons_data_.end(); ++dit) {
    curviol =(*dit)->getViol(x);   
    if (curviol > eTol_) {
      is_inf = false;
      curc_it = allCands.find((*dit)->getRInputVar());
      if (curc_it == allCands.end()) {
        allCands[(*dit)->getRInputVar()] = curviol;
      }
      else {
        (*curc_it).second = (*curc_it).second + curviol;
      }
    }	
  }
 
  // TODO: For now putting in all candidates, eventually probably want to choose a reasonable subset 
  for (curc_it = allCands.begin(); curc_it != allCands.end(); ++curc_it) {
     BrVarCandPtr br_can = (BrVarCandPtr) new BrVarCand((*curc_it).first,
                                                        (*curc_it).first->getIndex(),
                                                        (*curc_it).second,
                                                        (*curc_it).second);
     cands.insert(br_can);
  }


}


// Implement Handler::getBrMod().
ModificationPtr CxUnivarHandler::getBrMod(BrCandPtr cand, DoubleVector &x, 
                                          RelaxationPtr ,
                                          BranchDirection brdir) 
{

  // This method is used in Reliability branching.
  //XXX If you want it to be more accurate, you should add the new secant and linearizations
  //   into lmods.

  LinModsPtr lmods;
  lmods = (LinModsPtr) new LinMods();

  double minFromBds = 0.1;
  BrVarCandPtr vcand = boost::dynamic_pointer_cast <BrVarCand> (cand);
  VariablePtr v = vcand->getVar();

#if defined(DEBUG_CXUNIVARHANDLER)  
  std::cout << "Branching mod candidate (working problem) ID: " << v->getId() 
            << " address: " << (  v.get() ) << std::endl;
#endif  

  // x is a *relaxation* solution, while we have put the *original* (or
  // working) problem variables into the BrCandPtr, so we need to
  // update our value appropriately...
  
  double xval = x[v->getIndex()];
  double value = xval;  // Make sure branch value is not too close to an end point
  double len = v->getUb() - v->getLb();
  if (value < v->getLb() + minFromBds*len) {
    value = v->getLb() + minFromBds*len;
  } else if (value > v->getUb() - minFromBds*len) {
    value = v->getUb() - minFromBds*len; 
  }

  // can't branch on something that is at its bounds.
  if (!(value > v->getLb()+1e-8 && value < v->getUb()-1e-8)) {
    std::cerr << "Warning!  Branching on variable with bounds/value: [" << 
      v->getLb() << " , " << value << "  " << v->getUb() << " ]" << std::endl;
    //assert(value > v->getLb()+1e-8 && value < v->getUb()-1e-8);
  }

  if (brdir == DownBranch) {
    // down branch
    VarBoundModPtr mod = (VarBoundModPtr) new VarBoundMod(v, Upper, value);
    lmods->insert(mod);
  } else if (brdir ==  UpBranch) {
  // up branch
    VarBoundModPtr mod    = (VarBoundModPtr) new VarBoundMod(v, Lower, value);
    lmods->insert(mod);
  }

  return lmods;

}


// Implement Handler::getBranches().
Branches CxUnivarHandler::getBranches(BrCandPtr cand, DoubleVector &x,
                                      RelaxationPtr , SolutionPoolPtr )
{
  double minFromBds = 0.1;
  BrVarCandPtr vcand = boost::dynamic_pointer_cast <BrVarCand> (cand);
  VariablePtr v = vcand->getVar();

  double xval = x[v->getIndex()];
  double value = xval;  // Make sure branch value is not too close to an end point
  double len = v->getUb() - v->getLb();
  if (value < v->getLb() + minFromBds*len) {
    value = v->getLb() + minFromBds*len;
  } else if (value > v->getUb() - minFromBds*len) {
    value = v->getUb() - minFromBds*len; 
  }

  // can't branch on something that is at its bounds.
  if (!(value > v->getLb()+1e-8 && value < v->getUb()-1e-8)) {
    std::cerr << "Warning!  Branching on variable with bounds/value: [" << 
      v->getLb() << " , " << value << "  " << v->getUb() << " ]" << std::endl;
    //assert(value > v->getLb()+1e-8 && value < v->getUb()-1e-8);
  }

  Branches branches = (Branches) new BranchPtrVector();

  BranchPtr branch = (BranchPtr) new Branch();
  VarBoundModPtr mod = (VarBoundModPtr) new VarBoundMod(v, Upper, value);
  assert(!"add Mod correctly here.");
  branch->addPMod(mod);
  branch->setActivity((v->getUb()-value)/len);
  branches->push_back(branch);

  branch = (BranchPtr) new Branch();
  mod = (VarBoundModPtr) new VarBoundMod(v, Lower, value);
  assert(!"add Mod correctly here.");
  branch->addPMod(mod);
  branch->setActivity((value - v->getLb())/len);
  branches->push_back(branch);


  logger_->msgStream(LogDebug2) << "branching on " << v->getName();
  logger_->msgStream(LogDebug2) << " <= " << value << " or " 
    << " >= " << value << std::endl;

#if defined(DEBUG_CXUNIVARHANDLER)  
  std::cout << "branching on " << v->getName();
  std::cout << " <= " << value << " or " << " >= " << value << std::endl;
#endif
  
  return branches;

}


BranchPtr CxUnivarHandler::doBranch_(BranchDirection UpOrDown,
                                     ConstVariablePtr v, double bvalue)
{
  BranchPtr branch;
  LinModsPtr linmods;

#if defined(DEBUG_CXUNIVARHANDLER)
  std::cout << "CxUnivarHandler, Branching: " << (UpOrDown == DownBranch ? "Down" : "Up")
            << " at value: " << bvalue << " on: " << std::endl;
  v->write(std::cout);
#endif

  // Zero out tmpX and grad each time, or else bad things happen
  for (UInt j = 0; j < tmpX_.size(); ++j) {
    tmpX_[j] = 0.0;
    grad_[j] = 0.0;
  }
    
  branch = (BranchPtr) new Branch();
  linmods = (LinModsPtr) new LinMods();

  // Change bounds on the x var (called v here)
  if (UpOrDown == DownBranch) {

    VarBoundModPtr mod = (VarBoundModPtr) new VarBoundMod(v, Upper, bvalue);
    linmods->insert(mod);

    // Find *all* cons_data that has v as an input variable.
    CxUnivarConstraintIterator dit; 
  
    for (dit = cons_data_.begin(); dit != cons_data_.end(); ++dit) {
      if ((*dit)->getRInputVar() == v) {

	ConstVariablePtr rov = (*dit)->getROutVar();
	FunctionPtr fn = (*dit)->getOriginalCon()->getFunction();
	int error;

	// Change the secant constraint
	ConstraintPtr secCon = (*dit)->getSecantCon();

	LinearFunctionPtr lf; 
	FunctionPtr f;

	double xlb = v->getLb();
	double xub = bvalue;

	// TODO: Check the error value!
	tmpX_[v->getIndex()] = xlb;
	double fxlb =  fn->eval(tmpX_, &error);
	tmpX_[v->getIndex()] = xub;
	double fxub =  fn->eval(tmpX_, &error);
	tmpX_[v->getIndex()] = 0.0;

	// TODO: check/remedy numerical issues in this division
	double m = 0.0;
	if (xub - xlb > 10e-7) {
	  m = (fxub - fxlb)/(xub - xlb);
	}
	double intercept = fxlb - m*xlb;
	lf = (LinearFunctionPtr) new LinearFunction();
	lf->addTerm(rov, 1.0);
	lf->addTerm(v, -m);

        LinConModPtr lcmod = (LinConModPtr) new LinConMod(secCon, lf,
                                                          -INFINITY,
                                                          intercept);
	linmods->insert(lcmod);

	// Change all linearization constraints
	ConstraintVector::iterator lin_it;
        for(lin_it = (*dit)->linConsBegin(); lin_it != (*dit)->linConsEnd();
            ++lin_it) {
	  ConstraintPtr c = *lin_it;
	}
      }
    }
  } else {
    VarBoundModPtr mod = (VarBoundModPtr) new VarBoundMod(v, Lower, bvalue);
    linmods->insert(mod);
  }

  assert(!"add Mod correctly here.");
  branch->addPMod(linmods);
  return branch;

}


// presolve.
SolveStatus CxUnivarHandler::presolve(PreModQ *, bool *)
{
  return Finished;
}


// Implement Handler::presolveNode().
bool CxUnivarHandler::presolveNode(RelaxationPtr, NodePtr, SolutionPoolPtr,
                                   ModVector &, ModVector &)
{
  return false;
}


// Write name
std::string CxUnivarHandler::getName() const
{
  return "CxUnivarHandler (Handling univariate convex/concave terms).";
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
