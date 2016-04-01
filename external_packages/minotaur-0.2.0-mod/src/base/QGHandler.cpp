// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/** 
 * \file QGHandler.cpp
 * \Briefly define a handler for the textbook type Quesada-Grossmann
 * Algorithm.
 * \Authors Ashutosh Mahajan and Meenarli Sharma,Indian Institute of
 * Technology Bombay 
 */


#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "MinotaurConfig.h"

#include "CNode.h"
#include "Constraint.h"
#include "Engine.h"
#include "Environment.h"
#include "Function.h"
#include "Logger.h"
#include "Node.h"
#include "NonlinearFunction.h"
#include "Objective.h"
#include "Operations.h"
#include "Option.h"
#include "ProblemSize.h"
#include "QGHandler.h"
#include "Relaxation.h"
#include "Solution.h"
#include "SolutionPool.h"
#include "VarBoundMod.h"
#include "Variable.h"
#include "QuadraticFunction.h"

// #define SPEW 1

using namespace Minotaur;


typedef std::vector<ConstraintPtr>::const_iterator CCIter;
const std::string QGHandler::me_ = "QGHandler: ";

QGHandler::QGHandler()
: env_(EnvPtr()),      
  intTol_(1e-6),
  linCoeffTol_(1e-6),
  minlp_(ProblemPtr()),
  nlCons_(0),
  nlpe_(EnginePtr()),
  nlpStatus_(EngineUnknownStatus),
  nlpWs_(WarmStartPtr()),
  numCuts_(0),
  objVar_(VariablePtr()),
  oNl_(false),
  rel_(RelaxationPtr()),
  solAbsTol_(1e-5),
  solRelTol_(1e-5),
  stats_(0)
{
  logger_ = (LoggerPtr) new Logger(LogDebug2);
}

QGHandler::QGHandler(EnvPtr env, ProblemPtr minlp, EnginePtr nlpe) 
: env_(env),
  intTol_(1e-6),
  linCoeffTol_(1e-6),
  minlp_(minlp),
  nlCons_(0),
  nlpe_(nlpe),
  nlpStatus_(EngineUnknownStatus),
  nlpWs_(WarmStartPtr()),
  numCuts_(0),
  objVar_(VariablePtr()),
  oNl_(false),
  rel_(RelaxationPtr()),
  solAbsTol_(1e-5),
  solRelTol_(1e-5)
{
  logger_ = (LoggerPtr) new Logger((LogLevel)env->getOptions()->
                                   findInt("handler_log_level")->getValue());

  stats_   = new QGStats();
  stats_->nlpS = 0;
  stats_->nlpF = 0;
  stats_->nlpI = 0;
  stats_->cuts = 0;
}

QGHandler::~QGHandler()
{ 
  if (stats_) {
    delete stats_;
  }

  env_.reset();
  minlp_.reset();
  nlpe_.reset();
}

void QGHandler::addInitLinearX_(const double *x)
{ 
  ConstraintPtr con, newcon;
  double act = 0;
  double c;
  FunctionPtr f, f2;
  LinearFunctionPtr lf = LinearFunctionPtr();
  ObjectivePtr o;
  int error=0;

  for (CCIter it=nlCons_.begin(); it!=nlCons_.end(); ++it) { 
    con = *it;
    act = con->getActivity(x, &error);
    if(error==0) {

      f = con->getFunction();
      linearAt_(f, act, x, &c, &lf);
      f2 = (FunctionPtr) new Function(lf);
      if (con->getUb() < INFINITY) {
        newcon = rel_->newConstraint(f2, -INFINITY, con->getUb()-c,
                                     "lnrztn_cut");
        ++(stats_->cuts);
#if SPEW
        logger_->msgStream(LogDebug) << me_ << "initial constr. cut: ";
        newcon->write(logger_->msgStream(LogDebug));
#endif
      }

      if (con->getLb() > -INFINITY) {
        newcon = rel_->newConstraint(f2, con->getLb()-c, INFINITY, "lnrztn_cut");
        ++(stats_->cuts);  

#if SPEW
        logger_->msgStream(LogDebug) << me_ << "initial constr. cut: ";
        newcon->write(logger_->msgStream(LogDebug));
#endif
      }
    }	else {
      logger_->msgStream(LogError) << me_ 
        << "Constraint is not defined at this point" << std::endl;
    } 
  }
  error=0;  
  o = minlp_->getObjective();
  if (oNl_ && o) {
    f = o->getFunction();
    act = o->eval(x, &error);
    if(error==0) {
      linearAt_(f, act, x, &c, &lf);
      lf->addTerm(objVar_, -1.0);
      f2 = (FunctionPtr) new Function(lf);
      newcon = rel_->newConstraint(f2, -INFINITY, -1.0*c, "objlnrztn_cut");
      ++(stats_->cuts);
#if SPEW
      logger_->msgStream(LogDebug) << me_ << "initial obj cut: " << std::endl
        << std::setprecision(9);
      newcon->write(logger_->msgStream(LogDebug));
#endif
    }	else {
      logger_->msgStream(LogError) << me_ <<"Objective not defined at this point"
        << std::endl;
    }
  }	else {
    logger_->msgStream(LogDebug) << "QG Handler: Problem does not have a " 
      << "nonlinear objective." << std::endl;
  }
}

void QGHandler::cutIntSol_(ConstSolutionPtr sol, SolutionPoolPtr s_pool, 
                           bool *sol_found, SeparationStatus *status)
{
  const double *nlp_x;
  double nlpval = INFINITY;
  const double *x = sol->getPrimal();
  double lp_obj = (sol) ? sol->getObjValue() : -INFINITY;
  numCuts_=0;
  fixInts_(x);
  solveNLP_();
  unfixInts_();

  switch(nlpStatus_) {
  case (ProvenOptimal):
  case (ProvenLocalOptimal):

    ++(stats_->nlpF);
    updateUb_(s_pool, &nlpval, sol_found); 
#if SPEW
    logger_->msgStream(LogDebug) 
      << me_ << "solved fixed NLP to optimality, "
      << "lp_obj = " << lp_obj 
      << ", nlpval = " << nlpval << std::endl;
#endif
    if (lp_obj > nlpval-solAbsTol_ && 
        lp_obj > nlpval-(fabs(nlpval))*solRelTol_){

#if SPEW
      logger_->msgStream(LogDebug) 
        << me_ << "Pruned" << std::endl;
#endif         
      *status = SepaPrune;
      break;
    } else {
      relobj_=sol->getObjValue();
      nlp_x = nlpe_->getSolution()->getPrimal();
#if SPEW
      logger_->msgStream(LogDebug) 
        << me_ << "Not Pruned in OA" << std::endl;
#endif
      numCuts_=OAFromPoint_(nlp_x, sol->getPrimal(), status);
#if SPEW
      logger_->msgStream(LogDebug) 
        << me_ << "number of cuts from cutInt" << std::endl;
#endif 
      break;
    }
  case (ProvenInfeasible):
  case (ProvenLocalInfeasible): 
  case (ProvenObjectiveCutOff):
    ++(stats_->nlpI);
    nlp_x = nlpe_->getSolution()->getPrimal();
    relobj_ = sol->getObjValue();
    numCuts_=OAFromPointInf_(nlp_x, sol->getPrimal(), status);         
    break;
  case (ProvenUnbounded):
  case (EngineIterationLimit):
  case (ProvenFailedCQFeas):
  case (ProvenFailedCQInfeas):
  case (FailedFeas):
  case (FailedInfeas):
  case (EngineError):
  case (EngineUnknownStatus):
  default:
    logger_->msgStream(LogError) << me_ << "NLP engine status = " 
  << nlpe_->getStatusString() << std::endl
  << me_ << "No cut generated, may cycle!"
  << std::endl;
    *status = SepaError;

  }
  if (numCuts_ != 0) {
    *status = SepaResolve;
  }else{
#if SPEW
    logger_->msgStream(LogDebug) 
      << me_ << "number of cuts from OAFromPoint_ is 0" << std::endl;
#endif 
    *status = SepaPrune;
  }
}

void QGHandler::fixInts_(const double *x)
{
  VariablePtr v;
  double xval;
  VarBoundMod2 *m = 0;
  for (VariableConstIterator vit=minlp_->varsBegin(); vit!=minlp_->varsEnd(); 
       ++vit) {
    v = *vit;
    if (v->getType()==Binary || v->getType()==Integer) {
      xval = x[v->getIndex()];
      xval = floor(xval + 0.5); 
      m = new VarBoundMod2(v, xval, xval);
      m->applyToProblem(minlp_);
      nlpMods_.push(m);
    }
  }
}

void QGHandler::initLinear_(bool *isInf)
{

  const double *x;

  *isInf = false;

  nlpe_->load(minlp_);
  solveNLP_();

  switch (nlpStatus_) {
  case (ProvenOptimal):
  case (ProvenLocalOptimal):
    ++(stats_->nlpF);
    x = nlpe_->getSolution()->getPrimal();
    addInitLinearX_(x); 
    break;
  case (ProvenInfeasible):
  case (ProvenLocalInfeasible): 
  case (ProvenObjectiveCutOff):
    ++(stats_->nlpI);
    *isInf = true;
    break;
  case (ProvenUnbounded):
  case (EngineIterationLimit):
  case (ProvenFailedCQFeas):
  case (ProvenFailedCQInfeas):
  case (FailedFeas):
  case (FailedInfeas):
  case (EngineError):
  case (EngineUnknownStatus):
  default:
    logger_->msgStream(LogError) << me_ << "NLP engine status = " 
  << nlpStatus_ << std::endl;
    assert(!"QGHandler: Cannot proceed further");
  }
}

bool QGHandler::isFeasible(ConstSolutionPtr sol, RelaxationPtr rel, 
                           bool &, double & )
{

  FunctionPtr f;
  ConstraintPtr c;
  double act;
  const double *x = sol->getPrimal();
  int error=0;

  rel_ = rel;

  for (CCIter it=nlCons_.begin(); it!=nlCons_.end(); ++it) { 
    c = *it;
    act = c->getActivity(x, &error);
    if(error==0) {
      if( (act > c->getUb() + solAbsTol_ && 
           act > c->getUb() + (fabs(c->getUb())*solRelTol_))  ||
          ( act < c->getLb() - solAbsTol_ && 
            act < c->getLb() - (fabs(c->getLb())*solRelTol_)) ) {
#if SPEW
        logger_->msgStream(LogDebug) 
          << me_ << "constraint not feasible" << std::endl
          << me_;
        c->write(logger_->msgStream(LogDebug2));
        logger_->msgStream(LogDebug) 
          << me_ << "activity = " << act << std::endl;  
#endif
        return false;
      }
    }	else {
      logger_->msgStream(LogError) << me_ 
        << "Constraint not defined at this point"<< std::endl;
      return false;
    }
  }
#if SPEW
  logger_->msgStream(LogDebug) 
    << me_ << "all nonlinear constraints feasible." << std::endl;
#endif
  error=0;
  if (true == oNl_ ) {
    double alpha = x[objVar_->getIndex()]; 
    relobj_ = alpha;
    act = minlp_->getObjValue(x, &error);
    if(error==0){

      if (alpha < act - solAbsTol_ && alpha < act -(fabs(act)*solRelTol_) ) {
#if SPEW
        logger_->msgStream(LogDebug) << me_ << "objective not feasible" 
          << std::endl;
#endif
        return false;
      }
    }	else {
      logger_->msgStream(LogError) << me_ 
        <<" Objective not defined at this point"<< std::endl;
      return false;
    }
  }
  logger_->msgStream(LogDebug) << me_ << "Looks feasible" 
    << std::endl;
  return true;
}

void QGHandler::linearAt_(FunctionPtr f, double fval, const double *x, 
                          double *c, LinearFunctionPtr *lf)
{

  int n = rel_->getNumVars();
  double *a = new double[n];
  VariableConstIterator vbeg = rel_->varsBegin();
  VariableConstIterator vend = rel_->varsEnd();
  int error=0;

  std::fill(a, a+n, 0.);
  f->evalGradient(x, a, &error);
  if(error==0){
    *lf = (LinearFunctionPtr) new LinearFunction(a, vbeg, vend, linCoeffTol_); 
    *c  = fval - InnerProduct(x, a, numvars_);
    delete [] a;
  }	else {
    logger_->msgStream(LogError) << me_ <<" Gradient not defined at this point "
      <<  std::endl;
  }
}


void QGHandler::linearizeObj_(RelaxationPtr rel)
{
  ObjectivePtr o;
  o = minlp_->getObjective();
  std::string name     = "obj_dummy_var";
  if (!o) {
    assert(!"need objective in QG!");
  } else if (o->getFunctionType() != Linear) {
    VariablePtr vPtr     = rel->newVariable(-INFINITY, INFINITY, 
                                            Continuous, name); 
    LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();

    FunctionPtr f;

    assert(o->getObjectiveType()==Minimize);

    rel->removeObjective();
    lf->addTerm(vPtr, 1.0);
    f = (FunctionPtr) new Function(lf);
    rel->newObjective(f, 0.0, o->getObjectiveType());

    oNl_    = true;
    objVar_ = vPtr;
  }
}

int QGHandler::OAFromPoint_(const double *x, const double *inf_x,
                            SeparationStatus *status)
{
  double act=-INFINITY, nlpact = -INFINITY;
  ConstraintPtr con, newcon;
  double c;
  LinearFunctionPtr lf = LinearFunctionPtr(); 
  FunctionPtr f, f2;
  ObjectivePtr o;
  UInt num_cuts = 0;
  int error=0;
  double vio, lpvio;

  *status=SepaContinue;


  for (CCIter it=nlCons_.begin(); it!=nlCons_.end(); ++it) {
    con = *it; 
    f = con->getFunction();
    act = con->getActivity(x, &error); 
    nlpact =  con->getActivity(inf_x, &error);

    if(error==0){
      if(con->getUb() < INFINITY) {
        vio = std::max(nlpact-con->getUb(), 0.0);

        if (vio>solAbsTol_ && vio > (fabs(con->getUb())*solRelTol_) ) {
          linearAt_(f, act, x, &c, &lf);

          lpvio = std::max(lf->eval(inf_x)-con->getUb()+c, 0.0);

          if (lpvio>1e-4 && lpvio > (fabs(con->getUb()-c)*solRelTol_) ) {
            f2 = (FunctionPtr) new Function(lf);
            newcon = rel_->newConstraint(f2, -INFINITY, con->getUb()-c,
                                         "lnrztn_cut");
            ++num_cuts;
            *status = SepaResolve;
#if SPEW
            logger_->msgStream(LogDebug) << me_ <<" OA cut: " << std::endl
              << std::setprecision(9);
            newcon->write(logger_->msgStream(LogDebug));
#endif
          } else{
#if SPEW
            logger_->msgStream(LogDebug) << me_ <<" no OA cut added " << std::endl;
#endif
          }
        }
      }

      if(con->getLb() > -INFINITY) {
        vio = std::max(con->getLb()-nlpact, 0.0);

        if (vio>solAbsTol_ && vio > (fabs(con->getLb())*solRelTol_) ) {
          linearAt_(f, act, x, &c, &lf);
          lpvio = std::max(con->getLb()-c-lf->eval(inf_x), 0.0);

          if (lpvio>1e-4 && lpvio >(fabs(con->getLb()-c)*solRelTol_)) {
            f2 = (FunctionPtr) new Function(lf);

            newcon = rel_->newConstraint(f2, con->getLb()-c, INFINITY,
                                         "lnrztn_cut");
            ++num_cuts; 
            *status = SepaResolve;
#if SPEW
            logger_->msgStream(LogDebug) << me_ << "OA cut: " << std::endl
              << std::setprecision(9);
            newcon->write(logger_->msgStream(LogDebug));
#endif
          }
        }
      }
    }	else {
      logger_->msgStream(LogError) << me_ << "Constraint not defined at" <<
        " at least one of the two points: "<<  std::endl;
    }
  }
  o = minlp_->getObjective();
  if (oNl_ && o) {
    f = o->getFunction();
    act=o->eval(x,&error);
    nlpact = o->eval(inf_x, &error);
    if (error==0){

      vio = std::max(nlpact-relobj_, 0.0);

      if (vio > solAbsTol_ && vio > (fabs(relobj_)*solRelTol_) ) {

        linearAt_(f, act, x, &c, &lf);
        lpvio = std::max(c+lf->eval(inf_x)-relobj_, 0.0);

        if (lpvio>1e-4 && lpvio >(fabs(relobj_+c)*solRelTol_)) {
          lf->addTerm(objVar_, -1.0);
          f2 = (FunctionPtr) new Function(lf);
          newcon = rel_->newConstraint(f2, -INFINITY, -1.0*c, "objlnrztn_cut"); 
          ++num_cuts;
          *status = SepaResolve;
#if SPEW
          logger_->msgStream(LogDebug) << me_ << "OA cut: " << std::endl
            << std::setprecision(9);
          newcon->write(logger_->msgStream(LogDebug));
#endif
        } else{
#if SPEW
          logger_->msgStream(LogDebug) << me_ << "No objective OA cut: " << std::endl;
#endif

        }
      }
    }	else {
      logger_->msgStream(LogError) << me_ << "Objective not defined at this point" 
        <<  std::endl;
    }
  }

  return num_cuts;

}

int QGHandler::OAFromPointInf_(const double *x, const double *inf_x, 
                               SeparationStatus *status)
{
  int ncuts = 0;
  double act=-INFINITY;
  double lpact;
  ConstraintPtr con, newcon;
  double c;
  LinearFunctionPtr lf = LinearFunctionPtr(); 
  FunctionPtr f, f2;
  ObjectivePtr o;
  int error=0;

  *status=SepaContinue;

  for (CCIter it=nlCons_.begin(); it!=nlCons_.end(); ++it) {
    con = *it; 
    f = con->getFunction();
    act = con->getActivity(x, &error);
    if(error==0){
      if(con->getUb() < INFINITY) {
        linearAt_(f, act, x, &c, &lf);
        f2 = (FunctionPtr) new Function(lf);
        lpact = f2->eval(inf_x, &error);
        if (lpact - con->getUb() + c > solAbsTol_ && 
            lpact - con->getUb() + c >(fabs(con->getUb()-c)*solRelTol_)) {
          newcon = rel_->newConstraint(f2, -INFINITY, con->getUb()-c,
                                       "lnrztn_cut");
          ++ncuts;
          *status = SepaResolve;
#if SPEW
          logger_->msgStream(LogDebug) << me_ << "OA cut: " << std::endl
            << std::setprecision(9);
          newcon->write(logger_->msgStream(LogDebug));
#endif
        }
      }

      if(con->getLb() > -INFINITY) {
        linearAt_(f, act, x, &c, &lf);
        f2 = (FunctionPtr) new Function(lf);
        lpact = f2->eval(inf_x, &error);
        if (lpact - con->getLb() + c < -solAbsTol_  || 
            lpact - con->getLb() + c <-(fabs(con->getLb()-c)*solRelTol_)) {
          newcon = rel_->newConstraint(f2, con->getLb()-c, INFINITY,
                                       "lnrztn_cut");
          ++ncuts; 
          *status = SepaResolve;
#if SPEW
          logger_->msgStream(LogDebug) << me_ << "OA cut: " << std::endl
            << std::setprecision(9);
          newcon->write(logger_->msgStream(LogDebug));
#endif
        }
      }
    }
    else {
      logger_->msgStream(LogError) << me_ 
        << "Objective not defined at this point"<<  std::endl;
    }
  }

  return ncuts;
}

void QGHandler::relaxInitFull(RelaxationPtr , bool *)
{
  //Does nothing
}

void QGHandler::relaxInitInc(RelaxationPtr rel, bool *is_inf)
{
  relax_(rel, is_inf);
}

void QGHandler::relaxNodeFull(NodePtr , RelaxationPtr , bool *)
{

  //Does nothing
}

void QGHandler::relaxNodeInc(NodePtr , RelaxationPtr , bool *)
{
  //Does nothing
}

void QGHandler::relax_(RelaxationPtr rel, bool *is_inf)
{
  ConstraintPtr c;
  rel_ = rel;
  linearizeObj_(rel);
  numvars_ = minlp_->getNumVars();
  for (ConstraintConstIterator it=minlp_->consBegin(); it!=minlp_->consEnd(); 
       ++it) {

    c = *it;
    if (c->getFunctionType()!=Constant && c->getFunctionType() != Linear) {

      nlCons_.push_back(c);
    }
  }

#if SPEW
  logger_->msgStream(LogDebug) << me_ << "Number of nonlinear constraints "
    " = " << nlCons_.size() << std::endl;
  logger_->msgStream(LogDebug) << me_ << "Nonlinear solver used = "
    " = " << nlpe_->getName() << std::endl;
#endif
  initLinear_(is_inf);

#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "Initial relaxation:" 
    << std::endl;
  rel_->write(logger_->msgStream(LogDebug2));
#endif
}

void QGHandler::separate(ConstSolutionPtr sol, NodePtr , RelaxationPtr rel, 
                         CutManager *, SolutionPoolPtr s_pool,
                         bool *sol_found, SeparationStatus *status)
{      
  numvars_ = minlp_->getNumVars();
  VariableConstIterator v_iter;
  VariableType v_type;
  double value;
  const double *x = sol->getPrimal();
  bool is_int_feas = true; 
  rel_=rel;

  for (v_iter=rel_->varsBegin(); v_iter!=rel_->varsEnd(); 
       ++v_iter) {
    v_type = (*v_iter)->getType();
    if (v_type==Binary || v_type==Integer) {
      value = x[(*v_iter)->getIndex()];
      if (fabs(value - floor(value+0.5)) > intTol_) {
        is_int_feas = false;
        break;
      }
    }
  }


  if (is_int_feas) {
#if SPEW
    logger_->msgStream(LogDebug)<< me_ 
      << "solution is integer feasible, may need to add cuts" << std::endl;	  
#endif
    cutIntSol_(sol, s_pool, sol_found, status);
  } else{
#if SPEW
    logger_->msgStream(LogDebug) 
      << me_ << "solution is not integer feasible" << std::endl;	  
#endif
  }
}

void QGHandler::solveNLP_()
{
  nlpStatus_ = nlpe_->solve();
  ++(stats_->nlpS);
}

void QGHandler::updateUb_(SolutionPoolPtr s_pool, double *nlpval, 
                          bool *sol_found)
{
  double     val = nlpe_->getSolutionValue();
  double bestval = s_pool->getBestSolutionValue();


  if (val <= bestval) {
    const double *x = nlpe_->getSolution()->getPrimal();
#if SPEW
    logger_->msgStream(LogDebug) 
      << me_ << "new solution found, value = " << val << std::endl;
#endif
    s_pool->addSolution(x, val);
    *sol_found = true;
  }
  *nlpval = val;
}

void QGHandler::unfixInts_()
{
  Modification *m = 0;
  while(nlpMods_.empty() == false) {
    m = nlpMods_.top();
    m->undoToProblem(minlp_);
    nlpMods_.pop();
    delete m;
  }
}

void QGHandler::writeStats(std::ostream &out) const
{
  out
    << me_ << "number of nlps solved       = " << stats_->nlpS << std::endl
    << me_ << "number of infeasible nlps   = " << stats_->nlpI << std::endl
    << me_ << "number of feasible nlps     = " << stats_->nlpF << std::endl
    << me_ << "number of cuts added        = " << stats_->cuts << std::endl;
}

std::string QGHandler::getName() const
{
  return "QG Handler (Quesada-Grossmann)";
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
