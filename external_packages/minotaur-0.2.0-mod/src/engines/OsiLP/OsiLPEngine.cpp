// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file OsiLPEngine.cpp
 * \brief Implement an interface to the OSI-LP solver.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iostream>
#include <iomanip>

#if MNTROSICLP
#include "coin/OsiClpSolverInterface.hpp"
#endif

#if MNTROSICPX
#include "coin/OsiCpxSolverInterface.hpp"
#endif

#if MNTROSIGRB
#include "coin/OsiGrbSolverInterface.hpp"
#endif
#include "coin/CoinPackedMatrix.hpp"
#include "coin/CoinWarmStart.hpp"

#undef F77_FUNC_
#undef F77_FUNC

#include "MinotaurConfig.h"
#include "Constraint.h"
#include "Environment.h"
#include "Function.h"
#include "LinearFunction.h"
#include "HessianOfLag.h"
#include "Logger.h"
#include "Objective.h"
#include "Option.h"
#include "OsiLPEngine.h"
#include "Problem.h"
#include "Solution.h"
#include "Timer.h"
#include "Variable.h"

using namespace Minotaur;

//#define SPEW 1

const std::string OsiLPEngine::me_ = "OsiLPEngine: ";

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //

OsiLPWarmStart::OsiLPWarmStart()
  : coinWs_(0),
    mustDelete_(true)
{
}


OsiLPWarmStart::~OsiLPWarmStart()
{
  if (coinWs_ && mustDelete_) {
    delete coinWs_;
    coinWs_ = 0;
  }
}


bool OsiLPWarmStart::hasInfo()
{
  if (coinWs_) {
    return true;
  } else {
    return false;
  }
}


CoinWarmStart * OsiLPWarmStart::getCoinWarmStart() const
{
  return coinWs_;
}


void OsiLPWarmStart::setCoinWarmStart(CoinWarmStart *coin_ws, bool must_delete)
{
  if (coinWs_ && mustDelete_) {
    delete coinWs_;
    coinWs_ = 0;
  }

  coinWs_ = coin_ws;
  mustDelete_ = must_delete;
}


void OsiLPWarmStart::write(std::ostream &) const
{
  assert(!"implement me!");
}


// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //

OsiLPEngine::OsiLPEngine()
  : bndChanged_(true),
    consChanged_(true),
    env_(EnvPtr()),
    eName_(OsiUndefEngine),
    maxIterLimit_(10000),
    objChanged_(true),
    stats_(0),
    strBr_(false),
    timer_(0)
{
  logger_ = (LoggerPtr) new Logger(LogInfo);
#if USE_OSILP 
  osilp_ = newSolver_(OsiClpEngine);
#else 
#error Need to set USE_OSILP
#endif
  osilp_->setHintParam(OsiDoReducePrint);
  osilp_->messageHandler()->setLogLevel(0); 
#if MNTROSICLP
  OsiClpSolverInterface *osiclp = (OsiClpSolverInterface *)
    (dynamic_cast<OsiClpSolverInterface*>(osilp_));
  osiclp->setupForRepeatedUse();
#endif
}

  
OsiLPEngine::OsiLPEngine(EnvPtr env)
  : bndChanged_(true),
    consChanged_(true),
    env_(env),
    maxIterLimit_(10000),
    objChanged_(true),
    strBr_(false)
{
#if USE_OSILP
#else 
#error Need to set USE_OSILP
#endif
  std::string etype = env_->getOptions()->findString("lp_engine")->getValue();

  logger_ = (LoggerPtr) new Logger((LogLevel) env->getOptions()->
                                   findInt("engine_log_level")->getValue());
  eName_ = OsiUndefEngine;
  if (etype == "OsiClp") {
    eName_ = OsiClpEngine;
  } else if (etype == "OsiCpx") {
    eName_ = OsiCpxEngine;
  } else if (etype == "OsiGrb") {
    eName_ = OsiGrbEngine;
  }
  osilp_ = newSolver_(eName_);
  assert(osilp_);
  stats_ = new OsiLPStats();
  stats_->calls    = 0;
  stats_->strCalls = 0;
  stats_->time     = 0;
  stats_->strTime  = 0;
  stats_->iters    = 0;
  stats_->strIters = 0;

  timer_ = env->getNewTimer();

#if MNTROSICLP
  OsiClpSolverInterface *osiclp = (OsiClpSolverInterface *)
    (dynamic_cast<OsiClpSolverInterface*>(osilp_));
  osiclp->setupForRepeatedUse();
#endif
}


OsiLPEngine::~OsiLPEngine()
{
  delete osilp_;
  delete stats_;
  delete timer_;
  if (problem_) {
    problem_->unsetEngine();
    problem_.reset();
  }
}


void OsiLPEngine::addConstraint(ConstraintPtr con)
{
  LinearFunctionPtr lf = con->getLinearFunction();
  int nz = lf->getNumTerms();
  int *cols = new int[nz];
  double *elems = new double[nz];
  VariableGroupConstIterator it;
  int i=0;

  for (it = lf->termsBegin(); it != lf->termsEnd(); ++it, ++i){
    cols[i] = it->first->getIndex();
    elems[i] = it->second;
  }

  osilp_->addRow(nz, cols, elems, con->getLb(), con->getUb());
  //assert(!"implement me!");
  delete [] cols;
  delete [] elems;
  consChanged_ = true;
}


void OsiLPEngine::changeBound(ConstraintPtr cons, BoundType lu, double new_val)
{
  if (Upper==lu) {
    osilp_->setRowUpper(cons->getIndex(), new_val);
  } else {
    osilp_->setRowLower(cons->getIndex(), new_val);
  }
  bndChanged_ = true;
}


void OsiLPEngine::changeBound(VariablePtr var, BoundType lu, double new_val)
{
  //XXX: need a better map than the following for mapping variables to indices
  //and vice versa
  int col = var->getIndex();
  switch (lu) {
   case Lower:
     osilp_->setColLower(col, new_val);
     break;
   case Upper:
     osilp_->setColUpper(col, new_val);
     break;
   default:
     break;
  }
  bndChanged_ = true;
}


void OsiLPEngine::changeBound(VariablePtr var, double new_lb, double new_ub)
{
  int col = var->getIndex();
  osilp_->setColBounds(col, new_lb, new_ub);
  bndChanged_ = true;
}


#if MNTROSICLP
void OsiLPEngine::changeConstraint(ConstraintPtr c, LinearFunctionPtr lf, 
                                   double lb, double ub)
#else
void OsiLPEngine::changeConstraint(ConstraintPtr , LinearFunctionPtr , 
                                   double , double )
#endif
{
  if (eName_==OsiClpEngine) {
#if MNTROSICLP
    // OsiLPInterface does not have a modifyCoefficient function. So we have to
    // downcast it to OsiClpSolverInterface. XXX: Clean this code by creating a
    // map from constraints in problem, to those in engine. Then use deleteRow
    // and addRow, instead of modifyCoefficient().
    OsiClpSolverInterface *osiclp = (OsiClpSolverInterface *)
      (dynamic_cast<OsiClpSolverInterface*>(osilp_));
    int row = c->getIndex();
    ConstLinearFunctionPtr clf = c->getFunction()->getLinearFunction();

    // first zero out all the existing coefficients in the row.
    for (VariableGroupConstIterator it = clf->termsBegin(); it !=
         clf->termsEnd(); 
        ++it) {
      osiclp->modifyCoefficient(row, it->first->getIndex(), 0.0);
    }
    
    // assign new coefficients in the row.
    for (VariableGroupConstIterator it = lf->termsBegin(); it != lf->termsEnd(); 
        ++it) {
      osiclp->modifyCoefficient(row, it->first->getIndex(), it->second);
    }
    osiclp->setRowUpper(row, ub);
    osiclp->setRowLower(row, lb);
    consChanged_ = true;
#endif
  } else {
    assert(!"implement me!");
  }
}


void OsiLPEngine::changeConstraint(ConstraintPtr, NonlinearFunctionPtr)
{
    assert(!"Cannot change a nonlinear function in OsiLPEngine");
}


void OsiLPEngine::changeObj(FunctionPtr f, double)
{
  LinearFunctionPtr lf = f->getLinearFunction();
  double *obj = new double[problem_->getNumVars()];
  std::fill(obj, obj+problem_->getNumVars(),0.0);
  if (lf) {
    for (VariableGroupConstIterator it=lf->termsBegin(); it!=lf->termsEnd();
        ++it) {
      obj[it->first->getIndex()]  = it->second;
    }
  } 
  osilp_->setObjective(obj);
  objChanged_ = true;
  delete [] obj;
}


void OsiLPEngine::clear() {

  if (osilp_) {
    osilp_->reset();
    osilp_->setHintParam(OsiDoReducePrint);
    osilp_->messageHandler()->setLogLevel(0); 
  }
  if (problem_) {
    problem_->unsetEngine();
    problem_.reset();
  }
}


void OsiLPEngine::disableStrBrSetup() 
{
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "disabling strong branching." 
                               << std::endl;
#endif
  strBr_ = false;
}


EnginePtr OsiLPEngine::emptyCopy()
{
  if (env_) {
    return (OsiLPEnginePtr) new OsiLPEngine(env_);
  }
  return (OsiLPEnginePtr) new OsiLPEngine();
}


void OsiLPEngine::enableStrBrSetup() 
{
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "enabling strong branching."
                               << std::endl;
#endif
  strBr_ = true;
}


int OsiLPEngine::getIterationCount()
{
  return osilp_->getIterationCount();
}


std::string OsiLPEngine::getName() const
{
  return "OsiLP";
}


double OsiLPEngine::getSolutionValue() 
{
  return sol_->getObjValue();
}


ConstSolutionPtr OsiLPEngine::getSolution() 
{
  return sol_;
}


OsiSolverInterface * OsiLPEngine::getSolver() 
{
  return osilp_;
}


EngineStatus OsiLPEngine::getStatus() 
{
  return status_;
}


WarmStartPtr OsiLPEngine::getWarmStartCopy()
{
  // create a new copy of warm-start information from osilp_
  CoinWarmStart *coin_copy = osilp_->getWarmStart();

  OsiLPWarmStartPtr ws = (OsiLPWarmStartPtr) new OsiLPWarmStart();
  // save it. It is our responsibility to free it.
  ws->setCoinWarmStart(coin_copy, true);

  return ws;
}


void OsiLPEngine::load(ProblemPtr problem)
{
  problem_ = problem;
  int numvars = problem->getNumVars();
  int numcons = problem->getNumCons();
  int i,j;
  double obj_sense = 1.;
  CoinPackedMatrix *r_mat;
  double *conlb, *conub, *varlb, *varub, *obj;
  double *value;
  int *index;
  CoinBigIndex *start;

  ConstraintConstIterator c_iter;
  VariableConstIterator v_iter;
  
  conlb = new double[numcons];
  conub = new double[numcons];
  varlb = new double[numvars];
  varub = new double[numvars];

  VariableGroupConstIterator it;
  /* map the variables in this constraint to the function type (linear here) */

  //XXX Need to count the number of nnz in the problem 
  //     -- maybe add it to class later  
  LinearFunctionPtr lin;
  int nnz = 0;
  for (c_iter = problem->consBegin(); c_iter != problem->consEnd(); ++c_iter) {
    //XXX Don't want assert here, but not sure of eventually calling sequence
    //     and assumptions
    assert((*c_iter)->getFunctionType() == Linear);
    lin = (*c_iter)->getLinearFunction();
    nnz += lin->getNumTerms();
  }

  index = new int[nnz];
  value = new double[nnz];
  start = new CoinBigIndex[numcons+1];
    
  i = 0;
  j=0;
  start[0] = 0;
  for (c_iter = problem->consBegin(); c_iter != problem->consEnd(); ++c_iter) {
    conlb[i] = (*c_iter)->getLb();
    conub[i] = (*c_iter)->getUb();
    lin = (*c_iter)->getLinearFunction();
    for (it = lin->termsBegin(); it != lin->termsEnd(); ++it){
      ConstVariablePtr vPtr = it->first;
      index[j] = vPtr->getIndex();
      value[j] = it->second;
      ++j;
    }
    ++i;
    start[i]=j;
  }
  
  i = 0;
  for (v_iter=problem->varsBegin(); v_iter!=problem->varsEnd(); ++v_iter, 
       ++i) {
    varlb[i] = (*v_iter)->getLb();
    varub[i] = (*v_iter)->getUb();
  }

  // XXX: check if linear function is NULL
  lin = problem->getObjective()->getLinearFunction();
  if (problem->getObjective()->getObjectiveType() == Minotaur::Maximize) {
    obj_sense = -1.;
  }
  obj = new double[numvars];
  i = 0;
  if (lin) {
    for (v_iter=problem->varsBegin(); v_iter!=problem->varsEnd(); ++v_iter, 
         ++i) {
      obj[i]   = obj_sense*lin->getWeight(*v_iter);
    }
  } else {
    memset(obj, 0, numvars * sizeof(double));
  }

  r_mat = new CoinPackedMatrix(false, numvars, numcons, nnz, value, index, 
                               start, NULL);
  osilp_->loadProblem(*r_mat, varlb, varub, obj, conlb, conub);

  sol_ = (SolutionPtr) new Solution(1E20, 0, problem_);

  objChanged_ = true;
  bndChanged_ = true;
  consChanged_ = true;
  delete r_mat;
  delete [] index;
  delete [] value;
  delete [] start;
  delete [] conlb;
  delete [] conub;
  delete [] varlb;
  delete [] varub;
  delete [] obj;

  // osilp_->writeLp("stub");
  // exit(0);
  problem->setEngine(this);

}



void OsiLPEngine::loadFromWarmStart(const WarmStartPtr ws)
{
  ConstOsiLPWarmStartPtr ws2 = 
    boost::dynamic_pointer_cast <const OsiLPWarmStart> (ws);
  assert (ws2);
  CoinWarmStart *coin_ws = ws2->getCoinWarmStart();
  osilp_->setWarmStart(coin_ws);
}


void OsiLPEngine::negateObj()
{
  UInt n = problem_->getNumVars();
  double *obj = new double[n];
  const double* old_obj = osilp_->getObjCoefficients();
  std::copy(old_obj, old_obj+n, obj);
  osilp_->setObjective(obj);
  objChanged_ = true;
  delete [] obj;
}


OsiSolverInterface* OsiLPEngine::newSolver_(OsiLPEngineName ename)
{
  OsiSolverInterface *si = 0;
  switch (ename) {
  case (OsiClpEngine):
#if MNTROSICLP
    si = new OsiClpSolverInterface();
    si->setHintParam(OsiDoReducePrint);
    si->messageHandler()->setLogLevel(0); 
#else
    logger_->errStream()
      << me_ << "Minotaur is not compiled with OsiClp!" << std::endl;
#endif
    break;
  case (OsiCpxEngine):
#if MNTROSICPX
    si = new OsiCpxSolverInterface();
    si->setHintParam(OsiDoReducePrint);
    si->messageHandler()->setLogLevel(0); 
#else
    logger_->errStream()
      << me_ << "Minotaur is not compiled with OsiCpx!" << std::endl;
#endif
    break;
  case (OsiGrbEngine):
#if MNTROSIGRB
    si = new OsiGrbSolverInterface();
#else
    logger_->errStream()
      << me_ << "Minotaur is not compiled with OsiGrb!" << std::endl;
#endif
    break;
  default:
    break;
  }
  return si;
}


void OsiLPEngine::removeCons(std::vector<ConstraintPtr> &delcons)
{
  int num = delcons.size();
  int *inds = new int[num];
  for (int i=0; i<num; ++i) {
    inds[i] = delcons[i]->getIndex();
  }
  osilp_->deleteRows(num, inds);
  consChanged_ = true;
}


void OsiLPEngine::resetIterationLimit()
{
  OsiIntParam key = OsiMaxNumIteration;
  osilp_->setIntParam(key, maxIterLimit_);
}


void OsiLPEngine::setIterationLimit(int limit)
{
  OsiIntParam key = OsiMaxNumIteration;
  osilp_->setIntParam(key, limit);
}
  

EngineStatus OsiLPEngine::solve()
{
  timer_->start();
  if (true==objChanged_ && false==bndChanged_ && false==consChanged_) {
    osilp_->setHintParam(OsiDoDualInResolve, false);
  } else {
    osilp_->setHintParam(OsiDoDualInResolve, true);
  }

  stats_->calls += 1;
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "in call number " << stats_->calls
                               << std::endl;
#endif

  osilp_->resolve();

  if (osilp_->isProvenOptimal()) {
    status_ = ProvenOptimal;  
    sol_->setPrimal(osilp_->getStrictColSolution());
    sol_->setObjValue(osilp_->getObjValue()
        +problem_->getObjective()->getConstant());
    sol_->setDualOfCons(osilp_->getRowPrice());
    sol_->setDualOfVars(osilp_->getReducedCost());
  } else if (osilp_->isProvenPrimalInfeasible()) {
    status_ = ProvenInfeasible;
    sol_->setObjValue(INFINITY);
    sol_->setDualOfCons(osilp_->getRowPrice());
    sol_->setDualOfVars(osilp_->getReducedCost());
  } else if(osilp_->isProvenDualInfeasible()) {
    status_ = ProvenUnbounded;    // primal is not infeasible but dual is.
    sol_->setObjValue(-INFINITY);
    sol_->setDualOfCons(osilp_->getRowPrice());
    sol_->setDualOfVars(osilp_->getReducedCost());
  } else if (osilp_->isIterationLimitReached()) {
    status_ = EngineIterationLimit;
    sol_->setPrimal(osilp_->getStrictColSolution());
    sol_->setObjValue(osilp_->getObjValue()
        +problem_->getObjective()->getConstant());
    sol_->setDualOfCons(osilp_->getRowPrice());
    sol_->setDualOfVars(osilp_->getReducedCost());
  } else if(osilp_->isAbandoned()) {
    status_ = EngineError;
    sol_->setObjValue(INFINITY);
  } else if (osilp_->isPrimalObjectiveLimitReached()||
	   osilp_->isDualObjectiveLimitReached()) {
    status_ = ProvenObjectiveCutOff;
    sol_->setPrimal(osilp_->getStrictColSolution());
    sol_->setObjValue(osilp_->getObjValue()
        +problem_->getObjective()->getConstant());
    sol_->setDualOfCons(osilp_->getRowPrice());
    sol_->setDualOfVars(osilp_->getReducedCost());
  } else {
    status_ = EngineUnknownStatus;
    sol_->setObjValue(INFINITY);
  }

  stats_->iters += osilp_->getIterationCount();
  stats_->time  += timer_->query();
  if (strBr_) {
    ++(stats_->strCalls);
    stats_->strIters += osilp_->getIterationCount();
    stats_->strTime  += timer_->query();
  } 

#if SPEW
  logger_->msgStream(LogDebug) << me_ << "status = " << status_ << std::endl
                               << me_ << "solution value = " 
                               << sol_->getObjValue() << std::endl
                               << me_ << "iterations = " 
                               << osilp_->getIterationCount() << std::endl;
#endif
  timer_->stop();
  if (true==objChanged_ && false==bndChanged_ && false==consChanged_) {
    osilp_->setHintParam(OsiDoDualInResolve, true);
  }
  bndChanged_ = false;
  consChanged_ = false;
  objChanged_ = false;

  return status_;
}


void OsiLPEngine::writeLP(const char *filename) const 
{ 
  osilp_->writeLp(filename);
}


void OsiLPEngine::writeStats(std::ostream &out) const
{
  if (stats_) {
    std::string me = "OsiLP: ";
    out << me << "total calls            = " << stats_->calls << std::endl
      << me << "strong branching calls = " << stats_->strCalls << std::endl
      << me << "total time in solving  = " << stats_->time  << std::endl
      << me << "time in str branching  = " << stats_->strTime << std::endl
      << me << "total iterations       = " << stats_->iters << std::endl
      << me << "strong br iterations   = " << stats_->strIters << std::endl;
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
