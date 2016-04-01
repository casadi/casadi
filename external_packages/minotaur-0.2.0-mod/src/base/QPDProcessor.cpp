// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file QPDProcessor.cpp
 * \brief Define base class QPDProcessor.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath> // for INFINITY
#include <string.h> // for memset
#include <iomanip> 

#include "MinotaurConfig.h"
#include "Brancher.h"
#include "Constraint.h"
#include "Engine.h"
#include "Environment.h"
#include "Function.h"
#include "QPDProcessor.h"
#include "Handler.h"
#include "HessianOfLag.h"
#include "LinearFunction.h"
#include "Logger.h"
#include "Modification.h"
#include "Node.h"
#include "NonlinearFunction.h"
#include "Objective.h"
#include "Operations.h"
#include "Option.h"
#include "QuadraticFunction.h"
#include "Relaxation.h"
#include "SolutionPool.h"
#include "VarBoundMod.h"

using namespace Minotaur;

const std::string QPDProcessor::me_ = "QPDProcessor: ";

//#define SPEW 0

QPDProcessor::QPDProcessor()
  : branches_(Branches()),
    e_(EnginePtr()),
    env_(EnvPtr()), 
    eta_(VariablePtr()),
    etaL_(VariablePtr()),
    logger_(0),
    negDuals_(false),
    nlCons_(0),
    numSolutions_(0),
    p_(ProblemPtr()), 
    qp_(RelaxationPtr()),
    qpe_(EnginePtr()),
    solAbsTol_(1e-5),
    solRelTol_(1e-5),
    solveQPafNLP_(false),
    ubCon_(ConstraintPtr()),
    ws_(WarmStartPtr())
{
  stats_.bra = 0;
  stats_.inf = 0;
  stats_.nlp = 0;
  stats_.nlpI = 0;
  stats_.nlpU = 0;
  stats_.opt = 0;
  stats_.prob = 0;
  stats_.proc = 0;
  stats_.ub = 0;
}


QPDProcessor::QPDProcessor(EnvPtr env, ProblemPtr p, EnginePtr e, EnginePtr qe,
                           HandlerVector &handlers)
  : branches_(Branches()),
    e_(e),
    env_(env),
    eta_(VariablePtr()),
    etaL_(VariablePtr()),
    handlers_(handlers),
    nlCons_(0),
    numSolutions_(0),
    p_(p),
    qp_(RelaxationPtr()),
    qpe_(qe),
    solAbsTol_(1e-5),
    solRelTol_(1e-5),
    ubCon_(ConstraintPtr()),
    ws_(WarmStartPtr())
{
  stats_.bra = 0;
  stats_.cuts = 0;
  stats_.inf = 0;
  stats_.inf = 0;
  stats_.nlp = 0;
  stats_.nlpI = 0;
  stats_.nlpU = 0;
  stats_.opt = 0;
  stats_.prob = 0;
  stats_.proc = 0;
  stats_.sep = 0;
  stats_.ub = 0;
  logger_ = new Logger((LogLevel)
                       env->getOptions()->findInt("node_processor_log_level")
                       ->getValue());
  solveQPafNLP_ = (env->getOptions()->findString("nlp_engine")->getValue()
                   == "Filter-SQP");
  negDuals_ = solveQPafNLP_;
#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "initialized." << std::endl;
#endif 
}


QPDProcessor::~QPDProcessor()
{
  env_.reset();
  qp_.reset();
  p_.reset();
  nlCons_.clear();
  qpe_.reset();
  e_.reset();
  ws_.reset();
  if (logger_) {
    delete logger_;
    logger_ = 0;
  }
}


bool QPDProcessor::boundTooFar_(ConstSolutionPtr sol, NodePtr node,
                                double best) 
{
  int err = 0;
  double d;
  bool large = false;

  d = p_->getObjective()->eval(sol->getPrimal(), &err);
  assert(0==err);
  if (best < INFINITY) {
    if (d > best-0.001*fabs(best)) {
      large = true;
    } else if (d > 0.6*best + 0.4*node->getLb()) {
      large = true;
    } else if (fabs(node->getLb())>1e-3 && 
               d>node->getLb()+0.2*fabs(node->getLb())) {
      large = true;
    }
  } else if (fabs(node->getLb())<1e-3) {
    large = (d > 0.05); 
  } else {
    large = (d > node->getLb() + 0.2*fabs(node->getLb())); 
  }
#if SPEW
  logger_->msgStream(LogDebug2) << std::setprecision(6) << me_ 
                                << " obj at QPSOL = " << d
                                << " best known = " << best
                                << " node lb = " << node->getLb()
                                << " is large = " << large
                                << std::endl;
#endif
  return large;
}


void QPDProcessor::getLin_(FunctionPtr f, const double *x,
                           UInt n, VariableConstIterator vbeg,
                           VariableConstIterator vend, 
                           LinearFunctionPtr &lf, double &val)
{ 
  double *grad = new double[n+1]; 
  int err=0;

  val = 0.0;
  memset(grad, 0, (n+1) * sizeof(double));
  val = f->eval(x, &err); assert(0==err);
  f->evalGradient(x, grad, &err); assert(0==err);
  val -= InnerProduct(grad, x, n);
  lf = (LinearFunctionPtr) new LinearFunction(grad, vbeg, vend, 1e-7);

  delete [] grad;
  return;
}


void QPDProcessor::getObjLin_(NonlinearFunctionPtr nlf, const double *x,
                              UInt n, VariableConstIterator vbeg,
                              VariableConstIterator vend, 
                              LinearFunctionPtr &lf, double &val)
{ 
  double *grad = new double[n]; // n = no. of vars in qp.
  int err=0;

  // not a linearization at (x,eta)!. it is a linearization at (x, f(x)).

  val = 0.0;
  memset(grad, 0, n * sizeof(double));
  val = nlf->eval(x, &err); assert(0==err);
  nlf->evalGradient(x, grad, &err); assert(0==err);
  val -= InnerProduct(grad, x, n-1);
  lf = (LinearFunctionPtr) new LinearFunction(grad, vbeg, vend, 1e-7);
  lf->addTerm(eta_, -1.0);
  delete [] grad;
  return;
}


#if 0
  // check if violated.
  newact = lf->eval(x);
  if (newact < c->getUb()-fval+minvio && newact > c->getLb()-fval-minvio) {
#if SPEW
    logger_->msgStream(LogDebug2) << me_ << "linearization inequality is "
      << "not sufficiently violated " << newact << " bounds "
      << c->getLb()-fval << " " << c->getUb()-fval << std::endl;
#endif 
    return newc;
  } else {
#if SPEW
    logger_->msgStream(LogDebug2) << me_ << "inequality activity = "
      << newact << " bounds = " << c->getLb()-fval << " " << c->getUb()-fval
      << std::endl;
#endif 
  }


  f2 = (FunctionPtr) new Function(lf);
  if (c->getUb()>=INFINITY) {
    // f >= lb type of constraint.
    newc = rel->newConstraint(f2, c->getLb()-fval, INFINITY, c->getName());
  } else {
    newc = rel->newConstraint(f2, -INFINITY, c->getUb()-fval, c->getName());
  }
#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "added inequality ";
  newc->write(logger_->msgStream(LogDebug2));
#endif 

  assert(0==err);
  return newc;
#endif


void QPDProcessor::fixInts_(const double *x,
                            std::stack<Modification *> *nlp_mods)
{
  VariablePtr v;
  double xval;
  VarBoundMod2 *m = 0;
  for (VariableConstIterator vit=p_->varsBegin(); vit!=p_->varsEnd(); ++vit) {
    v = *vit;
    if (v->getType()==Binary || v->getType()==Integer) {
      xval = x[v->getIndex()];
      xval = floor(xval + 0.5); // round it to integer.
      m = new VarBoundMod2(v, xval, xval);
      m->applyToProblem(p_);
      nlp_mods->push(m);
    }
  }
}


bool QPDProcessor::foundNewSolution()
{
  return numSolutions_>0;
} 


Branches QPDProcessor::getBranches()
{
  ++stats_.bra;
  return branches_;
}


WarmStartPtr QPDProcessor::getWarmStart()
{
  return ws_;
}


bool QPDProcessor::isHFeasible_(ConstSolutionPtr sol, bool &should_prune)
{
  bool is_feas = true;
  HandlerIterator h;
  double inf_meas = 0.0;

  should_prune = false;
  for (h = handlers_.begin(); h != handlers_.end(); ++h) {
    is_feas = (*h)->isFeasible(sol, qp_, should_prune, inf_meas);
    if (is_feas == false || should_prune == true) {
      break;
    }
  }
#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "is integer = " << is_feas
    << std::endl;
#endif 


  return is_feas;
}


bool QPDProcessor::isNLPFeasible_(ConstSolutionPtr sol, double *vio)
{
  const double* x = sol->getPrimal();
  double act, vl, vu;
  int err = 0;
  ConstConstraintPtr c;
  bool is_feas = true;
  UInt i = 0;

  for (std::vector<ConstConstraintPtr>::const_iterator it=nlCons_.begin();
       it!=nlCons_.end(); ++it, ++i) {
    c = *it;
    act = c->getActivity(x, &err);
    assert(0==err);
    vl = c->getLb()-act;
    vu = act-c->getUb();
    vio[i] = (vl>vu) ? (vl>0?vl:0) : (vu>0?vu:0);
    if (vio[i]>1e-5) {
      is_feas = false;
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << std::setprecision(6) 
        << "constraint " << c->getName()
        << " is violated, activity = " << act << " bounds = " << c->getLb() << " " << 
       c->getUb() << std::endl;
#endif 
    }
#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "vio[" << i << "] = " << vio[i]
    << std::endl;
#endif 
  }

  // objective.
  if (eta_) {
    act = p_->getObjective()->getNonlinearFunction()->eval(x, &err);
    assert(0==err);
    vu = act - x[eta_->getIndex()];
    vio[i] = (vu>0.0) ? vu : 0.0;
  } else {
    vio[i] = 0.0;
  }
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "nonlinear constraints satisfied? = "
                               << is_feas << std::endl
                               << me_ << "obj vio = " << vio[i] 
                               << " " << qp_->getConstraint(qp_->getNumCons()-1)->getFunction()->eval(x, &err)
                               << " " << qp_->getConstraint(qp_->getNumCons()-1)->getUb() << " " 
//                               << " eta = " << x[eta_->getIndex()]
                               << " act = " << act << std::endl;
  //if (vio[i]>1.0) {
    //qp_->write(std::cout);
    //p_->getObjective()->getNonlinearFunction()->write(std::cout);
    //sol->write(std::cout);
    //exit(0);
  //}
#endif 
  return is_feas;
}


SolutionPtr QPDProcessor::nlpSol2Qp_(ConstSolutionPtr sol)
{
  double *x = new double[qp_->getNumVars()];
  SolutionPtr sol2;
  double objval = 0.0;
  int err = 0;

  memcpy(x, sol->getPrimal(), p_->getNumVars()*sizeof(double));
  if (eta_) {
    objval =  p_->getObjective()->getFunction()->getNonlinearFunction()
                ->eval(x, &err);
    assert(0==err);
    x[eta_->getIndex()] = objval;
  }
  objval =  qp_->getObjective()->eval(x, &err);
  assert(0==err);
  sol2 = (SolutionPtr) new Solution(objval, x, qp_);
  delete [] x;
  return sol2;
}


bool QPDProcessor::presolveNode_(NodePtr node, SolutionPoolPtr s_pool) 
{
  ModVector n_mods;      // Mods that are applied in this node.
  ModVector t_mods;      // Mods that need to be saved for subsequent nodes 
                         // as well. It is a subset of n_mods;
  ModificationPtr mod2;
  bool is_inf = false;

  // TODO: make this more sophisticated: loop several times until no more
  // changes are possible.
  for (HandlerIterator h = handlers_.begin(); h != handlers_.end() 
      && false==is_inf; ++h) {
    is_inf = (*h)->presolveNode(qp_, node, s_pool, n_mods, t_mods);
    for (ModificationConstIterator m_iter=t_mods.begin(); m_iter!=t_mods.end(); 
        ++m_iter) {
      node->addRMod(*m_iter);
      (*m_iter)->applyToProblem(qp_);
    }
    n_mods.clear();
    t_mods.clear();
  }

  if (is_inf) {
    node->setStatus(NodeInfeasible);
    ++stats_.inf;
  }

  //relaxation_->write(std::cout);
  //std::cout << "*** *** ***\n";
  return is_inf;
}


void QPDProcessor::processNLP_(NodePtr node, ConstSolutionPtr &sol,
                               ConstSolutionPtr qpsol,
                               SolutionPoolPtr s_pool, bool &should_prune)
{
  EngineStatus nlp_status;
  bool int_feas = false;

  should_prune = false;
  if (node->getDepth()>0 && qpsol) {
    p_->setInitialPoint(qpsol->getPrimal());
  }
  solveNLP_(sol, nlp_status);

  if (solveQPafNLP_ && qp_ && stats_.nlp>1) {
    qpe_->clear();
    qp_->setInitialPoint(qpsol->getPrimal());
    qpe_->load(qp_);
    qpe_->solve();
  }

  should_prune = shouldPrune_(node, nlp_status, sol, s_pool);
  if (should_prune) {
    return;
  }

  int_feas = isHFeasible_(sol, should_prune);
  if (should_prune) {
    return;
  }

  if (int_feas) {
    saveSol_(sol, s_pool, node);
    updateObjCons_(sol);
    node->setStatus(NodeHitUb);
    ++stats_.opt;
    ++stats_.ub;
    should_prune = true;
  }
}


void QPDProcessor::processQP_(UInt iter, NodePtr node, ConstSolutionPtr &sol,
                              SolutionPoolPtr s_pool, bool &should_prune,
                              bool &should_resolve)
{
  EngineStatus qp_status;
  SolutionPtr sol2;
  bool nlp_feas = false;
  bool int_feas = false;
  double *vio = 0;
  double tvio, maxvio;
  SeparationStatus status = SepaNone;
  ConstSolutionPtr nlp_sol = SolutionPtr();
  bool large_vio = false;
  bool large_objvio = false;
  ConstSolutionPtr qpsol = SolutionPtr();
  //bool do_obj = false;

  should_resolve = false;

  solveQP_(sol, qp_status);
  sol = qpe_->getSolution();
#if SPEW
  if (0==node->getId()) {
    logger_->msgStream(LogDebug2) << me_ << "actual value of QP = " 
                                << sol->getObjValue() << std::endl;
  }
#endif

  // update pseudo-cost from last branching.
  if (1==iter) {
    brancher_->updateAfterLP(node, sol);
  }


#if SPEW
  int err = 0;
  logger_->msgStream(LogDebug2) << me_ << "obj value of NLP at QP sol = " 
                              << p_->getObjective()->
                                 eval(sol->getPrimal(), &err) 
                              << std::endl;
  assert(0==err);
#endif

  // check if the QP is infeasible. Also set lb of node.
  should_prune = shouldPruneQP_(node, qp_status, sol, s_pool);
  if (should_prune) {
    return;
  }

  
  int_feas = isHFeasible_(sol, should_prune);
  if (should_prune) {
    return;
  }

  nlp_feas = isNLPFeasible_(sol, vio);

  if (nlp_feas && int_feas) {
    // int and nlp feasible.
    sol2 = translateSol_(sol);
    //std::cout << "sol2 value = " << sol2->getObjValue() << std::endl;
    if (sol2->getObjValue()<s_pool->getBestSolutionValue()) {
      saveSol_(sol2, s_pool, node);
      updateObjCons_(sol2);
    }
    if (node->getLb()>s_pool->getBestSolutionValue()-solAbsTol_) {
      node->setStatus(NodeOptimal);
      should_prune = true;
      return;
    } else {
      processNLP_(node, nlp_sol, qpsol, s_pool, should_prune); 
      sol = nlpSol2Qp_(nlp_sol);
      return;
    }
  } 

  if (eta_) {
    chkObjVio_(vio[nlCons_.size()], sol->getPrimal()[eta_->getIndex()],
               large_objvio);
    if (large_objvio) {
      processNLP_(node, nlp_sol, qpsol, s_pool, should_prune); 
      if (should_prune) {
        return;
      }
      separateObj_(sol, nlp_sol, vio[nlCons_.size()], &status);
      if (SepaPrune==status) {
        ++stats_.sep;
        should_prune = true;
        return;
      } else if (SepaResolve == status) {
        ++stats_.sep;
        should_resolve = true;
        return;
      }
    }
  }

  if (boundTooFar_(sol, node, s_pool->getBestSolutionValue())) {
    processNLP_(node, nlp_sol, qpsol, s_pool, should_prune); 
  }
  if (nlp_feas && !int_feas) {
    if (boundTooFar_(sol, node, s_pool->getBestSolutionValue())) {
      processNLP_(node, nlp_sol, qpsol, s_pool, should_prune); 
      sol = nlpSol2Qp_(nlp_sol);
      return;
    } else {
      return;
    }
  }

  chkVio_(node, vio, tvio, maxvio, large_vio);

  if (!large_vio && !int_feas) {
    return;
  }

  if (large_vio && !int_feas) {
    separateECP_(sol, vio, node, s_pool, &status);
    ++stats_.sep;
    if (SepaPrune==status) {
      should_prune = true;
    } else if (SepaResolve == status) {
      should_resolve = true;
    } 
    return;
  }

  assert(int_feas);
  processNLP_(node, nlp_sol, qpsol, s_pool, should_prune); 
  if (should_prune) {
    return;
  } else {
    separate_(true, sol, vio, nlp_sol, node, s_pool, &status);
    if (SepaPrune==status) {
      ++stats_.sep;
      should_prune = true;
    } else if (SepaResolve == status) {
      ++stats_.sep;
      should_resolve = true;
    } else {
      sol = nlpSol2Qp_(nlp_sol);
    }
  }
  return;
}


void QPDProcessor::processQP2_(UInt iter, NodePtr node, ConstSolutionPtr &sol,
                              SolutionPoolPtr s_pool, bool &should_prune,
                              bool &should_resolve)
{
  EngineStatus qp_status = EngineUnknownStatus;
  bool hiobjd = false;
  bool hifrac = false;
  bool hicvio = false;
  bool hietavio = false;
  bool nlp_feas = false;
  bool int_feas = false;
  const double *qpx = 0;
  ConstSolutionPtr nlp_sol = ConstSolutionPtr();
  ConstSolutionPtr sol2 = ConstSolutionPtr();
  SeparationStatus status = SepaNone;
  double *vio = new double[nlCons_.size()+1];

  should_resolve = false;
  solveQP_(sol, qp_status);
  qpx = sol->getPrimal();

#if SPEW
  if (0==node->getId()) {
    logger_->msgStream(LogDebug2) << me_ << "solution value of QP = " 
                                << sol->getObjValue() << std::endl;
  }
#endif

  // check if the QP is infeasible. Also set lb of node.
  should_prune = shouldPruneQP_(node, qp_status, sol, s_pool);
  if (should_prune) {
    goto CLEANUP;
  }

  // update pseudo-cost from last branching.
  if (1==iter) {
    brancher_->updateAfterLP(node, sol);
  }

  int_feas = isHFeasible2_(sol, hifrac, should_prune);
  chkObjVio2_(qpx, node, s_pool->getBestSolutionValue(), hiobjd, hietavio);
  nlp_feas = isNLPFeasible2_(qpx, vio, node, hicvio);
#if SPEW
  logger_->msgStream(LogDebug2) << me_ << std::endl 
                                << " hifrac = " << hifrac << std::endl
                                << " hicvio = " << hicvio << std::endl
                                << " hietavio = " << hietavio << std::endl
                                << " hiobjd = " << hiobjd << std::endl
                                << " int_feas = " << int_feas << std::endl
                                << " nlp_feas = " << nlp_feas << std::endl;
#endif

  if (nlp_feas && int_feas) {
    sol2 = translateSol_(sol);
    //std::cout << "sol2 value = " << sol2->getObjValue() << std::endl;
    if (sol2->getObjValue()<s_pool->getBestSolutionValue()) {
      saveSol_(sol2, s_pool, node);
      updateObjCons_(sol2);
    }
    if (node->getLb()>s_pool->getBestSolutionValue()) {
      node->setStatus(NodeOptimal);
      should_prune = true;
      goto CLEANUP;
    } else {
      //std::cout << "proc1\n";
      processNLP_(node, nlp_sol, sol, s_pool, should_prune); 
      if (!should_prune) {
        sol = nlpSol2Qp_(nlp_sol);
      }
      goto CLEANUP;
    } 
  }

  if (!nlp_feas && int_feas) {
    // TODO: it is a poor strategy, that should be improved.
    if (hiobjd) {
      //std::cout << "proc2\n";
      processNLP_(node, nlp_sol, sol, s_pool, should_prune); 
      if (should_prune) {
        goto CLEANUP;
      }
    }
    separateB_(sol, nlp_sol, vio, node, s_pool, &status);
    if (SepaResolve == status) {
      qp_->setInitialPoint(sol->getPrimal());
      should_resolve = true;
    } else if (!nlp_sol) {
      //std::cout << "proc2.5\n";
      processNLP_(node, nlp_sol, sol, s_pool, should_prune); 
      if (!should_prune) {
        sol = nlpSol2Qp_(nlp_sol);
      }
    } else {
      sol = nlpSol2Qp_(nlp_sol);
    }
    goto CLEANUP;
  }

  if (!hicvio && !hietavio && !hiobjd) {
    // branch
    goto CLEANUP;
  }
  if (!hifrac) {
    // branch
    goto CLEANUP;
  }

  if (hicvio && hietavio) {
    //std::cout << "proc4\n";
    //if (!nlp_sol) processNLP_(node, nlp_sol, sol, s_pool, should_prune);
    //if (should_prune) {
    //  goto CLEANUP;
    //}
    separateB_(sol, nlp_sol, vio, node, s_pool, &status);
    if (SepaResolve == status) {
      qp_->setInitialPoint(sol->getPrimal());
      should_resolve = true;
      goto CLEANUP;
    }
    sol = nlpSol2Qp_(nlp_sol);
    goto CLEANUP;
  }

  if (!hicvio && hietavio) {
    separateO_(sol, nlp_sol, vio, node, s_pool, &status);
    if (SepaResolve == status) {
      qp_->setInitialPoint(sol->getPrimal());
      should_resolve = true;
    }
    goto CLEANUP;
  }

  if (hicvio && !hietavio) {
    separateC_(sol, nlp_sol, vio, node, s_pool, &status);
    if (SepaResolve == status) {
      qp_->setInitialPoint(sol->getPrimal());
      should_resolve = true;
    }
    goto CLEANUP;
  }

  if (hiobjd) {
    goto CLEANUP;
    // do nothing
    //std::cout << "proc3\n";
    //processNLP_(node, nlp_sol, sol, s_pool, should_prune);
    //if (NodeInfeasible!=node->getStatus()) {
    //  separateB_(sol, nlp_sol, vio, node, s_pool, &status);
    //}
    //if (should_prune) {
    //  goto CLEANUP;
    //}
    //if (SepaResolve == status) {
    //  should_resolve = true;
    //  goto CLEANUP;
    //}
    //goto CLEANUP;
  }

  assert(!"should not be here!");
CLEANUP:
  delete [] vio;
  return;
}


void QPDProcessor::chkObjVio_(double vio, double etaval, bool &large_vio)
{
  double obj_large_thresh = 0.05; // important parameter.

  large_vio = false;

  if (eta_) {
    double d = vio;
    if (d<0.0) {
      d = 0.0;
    } else {
      d = d/fabs(d+etaval);
      if (d>obj_large_thresh) {
        large_vio = true;
      }
    }
  }
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "etaval " << etaval
                              << " obj vio = " << vio
                              << " large_vio = " << large_vio << std::endl;
#endif 
}


void QPDProcessor::chkVio_(NodePtr node, double *vio, double &tvio,
                           double &maxvio, bool &large_vio)
{
  double d;
  ConstConstraintPtr c;
  double large_thresh = (node->getDepth()<5)?0.05:0.2;
  UInt i = 0;

  large_vio = false;
  tvio = 0.0; maxvio = 0.0;
  for (i=0; i<nlCons_.size(); ++i) {
    c = nlCons_[i];
    d = vio[i];
    if (d<0.0) {
      d = 0.0;
    }
    tvio += d;
    if (d>maxvio) {
      maxvio = d;
    }
    if (c->getUb()<INFINITY) {
      d = (fabs(c->getUb())>1) ? d/fabs(c->getUb()) : d;
    } else {
      d = (fabs(c->getLb())>1) ? d/fabs(c->getLb()) : d;
    }
    if (d > large_thresh) {
      large_vio = true;
    }
  }

#if SPEW
  logger_->msgStream(LogDebug) << me_ << "id " << node->getId()
                              << " depth " << node->getDepth()
                              << " tvio " << tvio 
                              << " max vio = " << maxvio
                              << " large_vio = " << large_vio << std::endl;
#endif 
}


void QPDProcessor::processRootNode(NodePtr node, RelaxationPtr rel, 
                                   SolutionPoolPtr s_pool)
{
  // solve NLP first.
  ConstSolutionPtr sol;
  ConstSolutionPtr qpsol = SolutionPtr(); // NULL
  bool should_prune = false;
  BrancherStatus br_status;
  ModVector mods;
  UInt iter = 0;
  bool should_resolve = false;

#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "processing root node." << std::endl;
#endif 

  ++stats_.proc;
  numSolutions_ = 0;

  // important to initialize qp_ variables before processing nlp.
  qp_ = rel;
  qp_->newVariables(p_->varsBegin(), p_->varsEnd());

  processNLP_(node, sol, qpsol, s_pool, should_prune); 
  if (should_prune) {
    return;
  }

#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "branching in root node." << std::endl;
#endif 
  // get QP approx. for branching.
  setupQP_(sol);
  qpe_->load(qp_);

  while (true) {
    ++iter;

#if SPEW
  logger_->msgStream(LogDebug) <<  me_ << "iteration " << iter 
                               << std::endl;
#endif

    processQP2_(iter, node, sol, s_pool, should_prune, should_resolve);
    if (should_prune) {
      break;
    }

    if (false == should_resolve) {
      // save warm start information before branching. This step is expensive.
      ws_ = qpe_->getWarmStartCopy();
      branches_ = brancher_->findBranches(qp_, node, sol, s_pool, 
                                          br_status, mods);
      if (br_status==PrunedByBrancher) {
        should_prune = true;
        node->setStatus(NodeInfeasible);
        ++stats_.inf;
        break;
      } else if (br_status==ModifiedByBrancher) {
        for (ModificationConstIterator miter=mods.begin(); miter!=mods.end();
             ++miter) {
          // modify the relaxation only
          (*miter)->applyToProblem(qp_);
          node->addRMod(*miter);
        }
        should_prune = presolveNode_(node, s_pool);
        if (should_prune) {
          break;
        }
        should_resolve = true;
#if SPEW
        logger_->msgStream(LogDebug2) << me_ << "brancher returned a mod" 
                                      << std::endl;
#endif
      }
    } else {
      ++stats_.sep;
    }

    if (false == should_resolve) {
      break;
    }
  }
  return;
}


void QPDProcessor::process(NodePtr node, RelaxationPtr rel, 
                           SolutionPoolPtr s_pool)
{
  bool should_prune = true;
  bool should_resolve;
  BrancherStatus br_status;
  ConstSolutionPtr sol;
  ModVector mods;
  int iter = 0;

#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "processing node " << node->getId()
    << std::endl;
#endif 

  ++stats_.proc;
  qp_ = rel;
  numSolutions_ = 0;

  // TODO: is this check repeated?
  if (node->getParent()) {
    node->setLb(node->getParent()->getLb());
  }

  // main loop
  while (true) {
    ++iter;

#if SPEW
    logger_->msgStream(LogDebug) <<  me_ << "iteration " << iter 
                                 << std::endl;
#endif

    should_prune = presolveNode_(node, s_pool);
    if (should_prune) {
      break;
    }

    processQP2_(iter, node, sol, s_pool, should_prune, should_resolve);
    if (should_prune) {
      break;
    }

    if (false == should_resolve) {
      // save warm start information before branching. This step is expensive.
      ws_ = qpe_->getWarmStartCopy();
      branches_ = brancher_->findBranches(qp_, node, sol, s_pool, 
                                          br_status, mods);
      if (br_status==PrunedByBrancher) {
        should_prune = true;
        node->setStatus(NodeInfeasible);
        ++stats_.inf;
        break;
      } else if (br_status==ModifiedByBrancher) {
        for (ModificationConstIterator miter=mods.begin(); miter!=mods.end();
             ++miter) {
          (*miter)->applyToProblem(qp_);
        }
        should_resolve = true;
#if SPEW
        logger_->msgStream(LogDebug2) << me_ << "brancher returned a mod" 
                                      << std::endl;
#endif
      }
    } else {
      ++stats_.sep;
    }
    if (should_resolve == false) {
      break;
    }
  }
#if 0
  if ((true==should_prune || node->getLb() >-4150) && true==xfeas) {
    std::cout << "problem here!\n";
    std::cout << node->getStatus() << "\n";
    rel->write(std::cout);
    exit(0);
  }
#endif

  return;
}


void QPDProcessor::saveSol_(ConstSolutionPtr sol, SolutionPoolPtr s_pool,
                            NodePtr)
{
  s_pool->addSolution(sol);
  ++numSolutions_;
#if SPEW
  logger_->msgStream(LogDebug2) << me_
    << "found solution with value " << sol->getObjValue()
    << std::endl;
#endif 
}


#if 0
void QPDProcessor::saveQPSol_(ConstSolutionPtr sol, SolutionPoolPtr s_pool,
                              NodePtr node)
{
  SolutionPtr sol2;

  sol2 = translateSol_(sol);
  if (sol2->getObjValue() > s_pool->getBestSolutionValue()) {
    node->setStatus(NodeHitUb);
    ++stats_.ub;
#if SPEW
    logger_->msgStream(LogDebug2) << me_ << "found solution with value " <<
      sol2->getObjValue() << " sub-optimal!" << std::endl;
#endif 
    updateObjCons_(sol2);
  } else {
    s_pool->addSolution(sol2);
    ++numSolutions_;
    node->setStatus(NodeOptimal);
    ++stats_.opt;
#if SPEW
    logger_->msgStream(LogDebug2) << me_ << "found solution with value " <<
      sol2->getObjValue() << std::endl;
#endif 
    updateObjCons_(sol2);
  }
#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "pruning node " << std::endl;
#endif 
}
#endif


void QPDProcessor::separate_(bool is_nec, ConstSolutionPtr sol,
                             const double *vio,
                             ConstSolutionPtr nlp_sol,
                             NodePtr , SolutionPoolPtr ,
                             SeparationStatus *status) 
{
  double val;
  ConstConstraintPtr c;
  FunctionPtr f;
  UInt i = 0;
  UInt n = qp_->getNumVars();
  VariableConstIterator vbeg, vend;
  LinearFunctionPtr lf;
  double lfvio;
  double minlfvio = 1e-4;
  ConstraintPtr cnew;

  *status = SepaContinue;
  vbeg = qp_->varsBegin();
  vend = qp_->varsEnd();

  for (std::vector<ConstConstraintPtr>::const_iterator it=nlCons_.begin();
       it!=nlCons_.end(); ++it, ++i) {
    c = *it; 
    if (!shouldSep_(is_nec, vio[i], c)) {
      continue;
    }

    f = c->getFunction();
    getLin_(c->getFunction(), nlp_sol->getPrimal(), n, vbeg, vend, lf, val);
    lfvio = lf->eval(sol->getPrimal()) + val;
    if (c->getUb()<INFINITY) {
      lfvio = lfvio - c->getUb();
    } else {
      lfvio = c->getLb() - lfvio;
    }
    if (lfvio <= 0.0) {
      lfvio = 0.0;
    }


    if ((is_nec && lfvio > minlfvio) || (!is_nec && lfvio > 1e-4)) {
      f = (FunctionPtr) new Function(lf);
      cnew = qp_->newConstraint(f, c->getLb()-val, c->getUb()-val,
                                c->getName());
      ++stats_.cuts;
      *status = SepaResolve;
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "new inequality: ";
      cnew->write(logger_->msgStream(LogDebug2));
      logger_->msgStream(LogDebug2) << me_ << "lfvio = " << lfvio << std::endl;
#endif 
    } else {
#if SPEW
      //logger_->msgStream(LogDebug2) << me_ << "linear inequality not violated"
      //  << " activity = " << lf->eval(sol->getPrimal()) << " " 
      //  << c->getLb()-val << " " << c->getUb()-val << std::endl;
      //lf->write(std::cout);
      //std::cout << "val = " << val << std::endl;
      //c->write(std::cout);
      //for (VarSetConstIterator it=c->getFunction()->varsBegin(); it!=c->getFunction()->varsEnd(); ++it) {
      //  std::cout << (*it)->getName() << " " << nlp_sol->getPrimal()[(*it)->getIndex()] << " " << sol->getPrimal()[(*it)->getIndex()] << " rhs = " << c->getUb() - val << std::endl;
      //}
#endif 
      //getLin_(c->getFunction(), sol->getPrimal(), n, vbeg, vend, lf, val);
      //lfvio = lf->eval(sol->getPrimal()) + val;
      //std::cout << " new violation = " << lfvio << std::endl;
      //f = (FunctionPtr) new Function(lf);
      //cnew = qp_->newConstraint(f, c->getLb()-val, c->getUb()-val,
      //                          c->getName());
      //++stats_.cuts;
      //*status = SepaResolve;
    }
  }

#if 0
  if (do_obj) {
    getObjLin_(p_->getObjective()->getFunction()->getNonlinearFunction(),
               nlp_sol->getPrimal(), n, vbeg, vend, lf, val);
    lfvio = lf->eval(sol->getPrimal()) + val;
    if (lfvio > 1e-1) {
      f = (FunctionPtr) new Function(lf);
      cnew = qp_->newConstraint(f, -INFINITY, -val, "obj_lnrztn");
      ++stats_.cuts;
      *status = SepaResolve;
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "new obj inequality: ";
      cnew->write(logger_->msgStream(LogDebug2));
      logger_->msgStream(LogDebug2) << me_ << "lfvio = " << lfvio << std::endl;
#endif 
    } else {
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "obj linear inequality not violated"
        << " activity = " << lf->eval(sol->getPrimal())
        << " val = " << val << std::endl;
#endif 
    }
  }
#endif
}


void QPDProcessor::separateECP_(ConstSolutionPtr sol,
                                const double *vio,
                                NodePtr , SolutionPoolPtr ,
                                SeparationStatus *status) 
{
  double val;
  ConstConstraintPtr c;
  FunctionPtr f;
  UInt i = 0;
  UInt n = qp_->getNumVars();
  VariableConstIterator vbeg, vend;
  LinearFunctionPtr lf;
  double lfvio;
  double minlfvio = 1e-1;
  ConstraintPtr cnew;

  *status = SepaContinue;
  vbeg = qp_->varsBegin();
  vend = qp_->varsEnd();

  for (std::vector<ConstConstraintPtr>::const_iterator it=nlCons_.begin();
       it!=nlCons_.end(); ++it, ++i) {
    c = *it; 
    if (!shouldSep_(false, vio[i], c)) {
      continue;
    }

    f = c->getFunction();
    getLin_(c->getFunction(), sol->getPrimal(), n, vbeg, vend, lf, val);
    lfvio = lf->eval(sol->getPrimal()) + val;
    if (c->getUb()<INFINITY) {
      lfvio = lfvio - c->getUb();
    } else {
      lfvio = c->getLb() - lfvio;
    }
    if (lfvio <= 0.0) {
      lfvio = 0.0;
    }
    if (lfvio > minlfvio) {
      f = (FunctionPtr) new Function(lf);
      cnew = qp_->newConstraint(f, c->getLb()-val, c->getUb()-val,
                                c->getName());
      ++stats_.cuts;
      *status = SepaResolve;
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "new ECP inequality: ";
      cnew->write(logger_->msgStream(LogDebug2));
      logger_->msgStream(LogDebug2) << me_ << "lfvio = " << lfvio << std::endl;
#endif 
    } else {
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "ECP linear inequality not "
        << "violated. activity = " << lf->eval(sol->getPrimal()) << " " 
        << c->getLb()-val << " " << c->getUb()-val << std::endl;
      lf->write(logger_->msgStream(LogDebug2));
      logger_->msgStream(LogDebug2) << "val = " << val << std::endl;
      c->write(logger_->msgStream(LogDebug2));
      for (VarSetConstIterator vit=c->getFunction()->varsBegin();
           vit!=c->getFunction()->varsEnd(); ++vit) {
        logger_->msgStream(LogDebug2) << (*vit)->getName() << " " 
          << sol->getPrimal()[(*vit)->getIndex()] << " rhs = " 
          << c->getUb() - val << std::endl;
      }
#endif 
    }
  }
}


void QPDProcessor::separateObj_(ConstSolutionPtr sol, ConstSolutionPtr nlp_sol,
                                double vio, SeparationStatus *status) 
{
  double val;
  FunctionPtr f;
  UInt n = qp_->getNumVars();
  VariableConstIterator vbeg, vend;
  LinearFunctionPtr lf;
  double lfvio;
  ConstraintPtr cnew;

  *status = SepaContinue;
  vbeg = qp_->varsBegin();
  vend = qp_->varsEnd();

  if (eta_ && shouldSepObj_(false, vio, sol->getPrimal()[eta_->getIndex()])) {
    getObjLin_(p_->getObjective()->getFunction()->getNonlinearFunction(),
               nlp_sol->getPrimal(), n, vbeg, vend, lf, val);
    lfvio = lf->eval(sol->getPrimal()) + val;
    if (lfvio > 1e-1) {
      f = (FunctionPtr) new Function(lf);
      cnew = qp_->newConstraint(f, -INFINITY, -val, "obj_lnrztn");
      ++stats_.cuts;
      *status = SepaResolve;
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "new obj inequality: ";
      cnew->write(logger_->msgStream(LogDebug2));
      logger_->msgStream(LogDebug2) << me_ << "lfvio = " << lfvio << std::endl;
#endif 
    } else {
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "obj linear inequality not violated"
        << " activity = " << lf->eval(sol->getPrimal())
        << " val = " << val << std::endl;
#endif 
    }
  }
}


void QPDProcessor::setupQP_(ConstSolutionPtr sol)
{
  ConstConstraintPtr c, cnew;
  FunctionPtr f;
  LinearFunctionPtr lf, lf2;
  VariableConstIterator vbeg, vend;
  ObjectivePtr obj = ObjectivePtr();
  int err = 0;
  UInt n = p_->getNumVars();
  double *grad = new double[n+1];
  double *hess = 0;
  double *hdotx = 0;
  UInt *jcol = 0;
  UInt *irow = 0;
  QuadraticFunctionPtr qf;
  UInt hess_nz = 0;
  double pred = 0; // for debug.
  double obj_const = 0;
  UInt zduals = 0;
  NonlinearFunctionPtr nlf;
  double val;
  SolutionPtr qpsol;
  double *ddd = new double[p_->getNumCons()];

  //sol->write(std::cout);
#if SPEW
  logger_->msgStream(LogInfo) << me_ << "nlp optimal value = " 
                              << sol->getObjValue() << std::endl;
  for (ConstraintConstIterator it=p_->consBegin(); it!=p_->consEnd(); ++it) {
    c = *it;
    if (Linear == c->getFunction()->getType() || 
        Constant == c->getFunction()->getType()) {
    } else {
      logger_->msgStream(LogDebug2) << c->getName() << " " << c->getActivity(sol->getPrimal(), &err) << std::endl
       << c->getName() << " " << c->getFunction()->getNonlinearFunction()->eval(sol->getPrimal(), &err) << std::endl;
    }
  }
#endif

  vbeg = qp_->varsBegin();
  vend = qp_->varsEnd();

  obj = p_->getObjective();
  if (obj) {
    f = obj->getFunction();
    if (f) {
      lf2 = f->getLinearFunction();
      if (lf2) {
        lf = lf2->cloneWithVars(vbeg);
        pred += lf2->eval(sol->getPrimal());
      } else {
        lf = (LinearFunctionPtr) new LinearFunction();
      }

      // qf not allowed yet.
      assert(!f->getQuadraticFunction());

      nlf = f->getNonlinearFunction();
      if (nlf) {
        eta_ = qp_->newVariable(-INFINITY, INFINITY, Continuous, "eta");
        lf->addTerm(eta_, 1.0);
        pred += nlf->eval(sol->getPrimal(), &err); assert(0==err);
        vbeg = qp_->varsBegin();
        vend = qp_->varsEnd();
      } 

      hess_nz = p_->getHessian()->getNumNz();
      hess = new double[hess_nz];
      hdotx = new double[n];
      jcol = new UInt[hess_nz];
      irow = new UInt[hess_nz];
      memset(hess, 0, hess_nz * sizeof(double));
      memset(hdotx, 0, n * sizeof(double));
      p_->getHessian()->fillRowColIndices(irow, jcol) ;

      memcpy(ddd, sol->getDualOfCons(), p_->getNumCons()*sizeof(double));
      for (UInt iii=0; iii<p_->getNumCons(); ++iii) {
        if (negDuals_) {
          ddd[iii] *= -1.0;
        }
        if (fabs(ddd[iii])<1e-4) {
          if (p_->getConstraint(iii)->getUb()<INFINITY) {
            ddd[iii] = 0.5; // c(x) <= ub type constraint
          } else {
            ddd[iii] = -0.5;// c(x) >= lb type constraint
          }
        }
      }
      p_->getHessian()->fillRowColValues(sol->getPrimal(), 1.0,
                                         ddd, hess, &err);
      assert(0==err);
      qf = (QuadraticFunctionPtr) new QuadraticFunction(hess_nz, hess, irow,
                                                        jcol, vbeg);
      if (qf->getNumTerms()==0) {
        qf.reset();
#if SPEW
        logger_->msgStream(LogInfo) << me_ << "qf is empty!" << std::endl; 
#endif
      }
      symMatDotV(hess_nz, hess, irow, jcol, sol->getPrimal(), hdotx);
      obj_const += 0.5*InnerProduct(hdotx, sol->getPrimal(), n);;
      pred -= 0.5*InnerProduct(hdotx, sol->getPrimal(), n);
      for (UInt i=0; i<n; ++i) {
        lf->incTerm(*(vbeg+i),-hdotx[i]);
      }

      f = (FunctionPtr) new Function(lf, qf);
      qp_->newObjective(f, 0, Minimize, "qp_obj");
      delete [] hess;
      delete [] hdotx;
      delete [] irow;
      delete [] jcol;
      hess = hdotx = 0;
      irow = jcol = 0;
    }
  }


  for (ConstraintConstIterator it=p_->consBegin(); it!=p_->consEnd(); ++it) {
    c = *it;
    if (Linear == c->getFunction()->getType() || 
        Constant == c->getFunction()->getType()) {
      f = c->getFunction()->cloneWithVars(vbeg, &err);
      assert(err==0);
      cnew = qp_->newConstraint(f, c->getLb(), c->getUb(), c->getName());
    } else {
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "dual multiplier for constraint "
                                    << c->getName() << " = " 
                                    << sol->getDualOfCons()[c->getIndex()]
                                    << std::endl;
#endif
      if (fabs(sol->getDualOfCons()[c->getIndex()])<1e-4) {
        ++zduals;
      }
      nlCons_.push_back(c);
      getLin_(c->getFunction(), sol->getPrimal(), n, vbeg, vend, lf, val);
      f = (FunctionPtr) new Function(lf);
      cnew = qp_->newConstraint(f, c->getLb()-val, c->getUb()-val,
                                c->getName());
    }
  }

  if (eta_) {
    val = 0;
    memset(grad, 0, (n+1) * sizeof(double));
    p_->getObjective()->getFunction()->getNonlinearFunction()->
      evalGradient(sol->getPrimal(), grad, &err);
    assert(0==err);
    lf = (LinearFunctionPtr) new LinearFunction(grad, vbeg, vend, 1e-7);
    lf->addTerm(eta_, -1.0);
    f = (FunctionPtr) new Function(lf);
    val = nlf->eval(sol->getPrimal(), &err); assert(0==err);
    qp_->newConstraint(f, -INFINITY, InnerProduct(grad,
                                                  sol->getPrimal(),n)-val,
                       "obj_cons");
  }

  // set initial point.
  qpsol = nlpSol2Qp_(sol);
  qp_->setInitialPoint(qpsol->getPrimal());

  qp_->setNativeDer();

#if SPEW
  logger_->msgStream(LogInfo) << me_ << "predicted value of f(qp) = " 
                              << pred << std::endl
                              << "obj_const = " << obj_const << std::endl
                              << "zduals = " << zduals << std::endl;
#endif
  //p_->write(std::cout);
  //sol->write(std::cout);
  //qp_->write(std::cout);
  //exit(0);
  
  delete [] grad;
  delete [] ddd;
}


bool QPDProcessor::shouldPrune_(NodePtr node, EngineStatus status, 
                                ConstSolutionPtr sol, SolutionPoolPtr s_pool)
{
  bool should_prune = false;
  switch (status) {
   case (FailedInfeas):
   case (ProvenFailedCQInfeas):
     ++stats_.prob;
   case (ProvenInfeasible):
   case (ProvenLocalInfeasible):
     node->setStatus(NodeInfeasible);
     ++stats_.inf;
     ++stats_.nlpI;
     should_prune = true;
     break;

   case (ProvenObjectiveCutOff):
     // update stats
     node->setStatus(NodeHitUb);
     should_prune = true;
     ++stats_.ub;
     break;

   case (ProvenUnbounded):
     should_prune = false;
     logger_->msgStream(LogDebug2) << me_ << "relaxation is unbounded!" 
                                   << std::endl;
     break;

   case (FailedFeas):
     node->setLb(node->getParent()->getLb());
     node->setStatus(NodeContinue);
     ++stats_.prob;
     break;

   case (ProvenFailedCQFeas):
     if (node->getParent()) {
       node->setLb(node->getParent()->getLb());
     } else {
       node->setLb(-INFINITY);
     }
     node->setStatus(NodeContinue);
     ++stats_.prob;
     break;
   case (ProvenLocalOptimal):
   case (ProvenOptimal):
     node->setLb(sol->getObjValue());
     if (sol->getObjValue() > s_pool->getBestSolutionValue()-solAbsTol_) {
       node->setStatus(NodeHitUb);
       should_prune = true;
       ++stats_.ub;
     } else {
       should_prune = false;
       node->setStatus(NodeContinue);
     }
     break;
   case (EngineIterationLimit):
     node->setStatus(NodeInfeasible);
     should_prune = true;
     ++stats_.inf;
     ++stats_.nlpI;
     break;
   case (EngineError):
     node->setStatus(NodeContinue);
     if (node->getParent()) {
       node->setLb(node->getParent()->getLb());
     } else {
       node->setLb(-INFINITY);
     }
     ++stats_.prob;
     break;
   default:
     break;
  }

  return should_prune;
}


bool QPDProcessor::shouldPruneQP_(NodePtr node, EngineStatus status, 
                                  ConstSolutionPtr sol, SolutionPoolPtr )
{
  bool should_prune = false;
  int err=0;

  switch (status) {
   case (FailedInfeas):
   case (ProvenFailedCQInfeas):
     ++stats_.prob;
   case (ProvenInfeasible):
   case (ProvenLocalInfeasible):
     node->setStatus(NodeInfeasible);
     ++stats_.inf;
     should_prune = true;
     break;

   case (ProvenObjectiveCutOff):
     // update stats
     node->setStatus(NodeHitUb);
     should_prune = true;
     ++stats_.ub;
     break;

   case (ProvenUnbounded):
     should_prune = false;
     logger_->msgStream(LogDebug2) << me_ << "relaxation is unbounded!" 
                                   << std::endl;
     break;

   case (FailedFeas):
     node->setStatus(NodeContinue);
     ++stats_.prob;
     if (node->getParent()) {
       node->setTbScore(node->getParent()->getTbScore());
     }
     break;

   case (ProvenFailedCQFeas):
     node->setStatus(NodeContinue);
     ++stats_.prob;
     if (node->getParent()) {
       node->setTbScore(node->getParent()->getTbScore());
     }
     break;
   case (ProvenLocalOptimal):
   case (ProvenOptimal):
     should_prune = false;
     node->setStatus(NodeContinue);
     //node->setTbScore(sol->getObjValue());
     node->setTbScore(p_->getObjective()->eval(sol->getPrimal(), &err));
     assert(0==err);
     break;
   case (EngineIterationLimit):
     node->setStatus(NodeInfeasible);
     should_prune = true;
     if (node->getParent()) {
       node->setTbScore(node->getParent()->getTbScore());
     }
     ++stats_.prob;
     break;
   case (EngineError):
     node->setStatus(NodeContinue);
     ++stats_.prob;
     break;
   default:
     break;
  }
  return should_prune;
}


bool QPDProcessor::shouldSep_(bool is_nec, double vio, ConstConstraintPtr c)
{
  if (is_nec) {
    return vio > 1e-5;
  } else {
    double d = (c->getUb() < INFINITY) ? c->getUb() : c->getLb();
    assert(c->getUb() >= INFINITY || c->getLb() <= INFINITY);
    d = (fabs(d)>1.0) ? fabs(d) : 1.0;
    return vio/d > 2e-2;
  }
}


bool QPDProcessor::shouldSepObj_(bool is_nec, double vio, double etaval)
{
  if (is_nec) {
    return vio > 1e-5;
  } else {
    double d = vio/fabs(vio+etaval);
    return vio/d > 0.2;
  }
}


void QPDProcessor::solveNLP_(ConstSolutionPtr &sol, EngineStatus &nlp_status)
{
  // first copy bounds from qp to nlp
  Modification *m = 0;
  std::stack<Modification *> nlpMods;
  VariablePtr var, qpvar;

  for (VariableConstIterator vit=p_->varsBegin(); vit!=p_->varsEnd(); 
       ++vit) {
    var = *vit;
    qpvar = qp_->getVariable(var->getIndex());
    m = new VarBoundMod2(var, qpvar->getLb(), qpvar->getUb());
    m->applyToProblem(p_);
    nlpMods.push(m);
  }

  nlp_status = e_->solve();
  sol = e_->getSolution();
  ++(stats_.nlp);
#if SPEW
  logger_->msgStream(LogDebug) 
    << me_ << "NLP status = " << e_->getStatusString() << std::endl 
    << me_ << "solution value = " << sol->getObjValue() << std::endl; 
#endif 

  // undo all bound changes. // can be made more efficient.
  while(nlpMods.empty() == false) {
    m = nlpMods.top();
    m->undoToProblem(p_);
    nlpMods.pop();
    delete m;
  }
}


void QPDProcessor::solveQP_(ConstSolutionPtr &sol, EngineStatus &qp_status)
{
  qp_status = qpe_->solve();
  sol = qpe_->getSolution();
#if SPEW
  logger_->msgStream(LogDebug) 
    << me_ << "QP status = " << qpe_->getStatusString() << std::endl 
    << me_ << "solution value = " << sol->getObjValue() << std::endl; 
#endif 

}


SolutionPtr QPDProcessor::translateSol_(ConstSolutionPtr sol)
{
  // TODO: replace with a clone function.
  SolutionPtr sol2 = (SolutionPtr) new Solution(sol);
  int err = 0;

  sol2->setObjValue(p_->getObjective()->eval(sol->getPrimal(), &err)); 
  assert(0==err);
  return sol2;
}


void QPDProcessor::unfixInts_(std::stack<Modification *> *nlp_mods)
{
  Modification *m = 0;
  while(nlp_mods->empty() == false) {
    m = nlp_mods->top();
    m->undoToProblem(p_);
    nlp_mods->pop();
    delete m;
  }
}


void QPDProcessor::updateObjCons_(ConstSolutionPtr sol)
{

  if (ubCon_) {
    if (sol->getObjValue() < ubCon_->getUb()-1e-5) {
      qp_->changeBound(ubCon_, Upper, sol->getObjValue()-1e-5);
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "updated ub constraint: ";
      ubCon_->write(logger_->msgStream(LogDebug2));
#endif 
    }
  } else {
    LinearFunctionPtr lf;
    FunctionPtr f = p_->getObjective()->getFunction();
    if (!f) {
      return;
    }

    lf = f->getLinearFunction();
    if (!lf && !eta_) {
      return;
    }

    if (lf) {
      lf = lf->cloneWithVars(qp_->varsBegin());
      if (eta_) {
        lf->addTerm(eta_, 1.0);
      }
      f = (FunctionPtr) new Function(lf);
      ubCon_ = qp_->newConstraint(f, -INFINITY, sol->getObjValue()-1e-5, "ub_cons");
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "new ub constraint: ";
      ubCon_->write(logger_->msgStream(LogDebug2));
#endif 
    } else {
      qp_->changeBound(eta_, Upper, sol->getObjValue()-1e-5);
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "new ub constraint: ";
      eta_->write(logger_->msgStream(LogDebug2));
#endif 
    }
  }
}


bool QPDProcessor::isHFeasible2_(ConstSolutionPtr sol, bool &ishigh,
                                 bool &should_prune)
{
  bool is_feas = true;
  HandlerIterator h;
  double inf_meas = 0.0;

  should_prune = false;
  for (h = handlers_.begin(); h != handlers_.end(); ++h) {
    is_feas = (*h)->isFeasible(sol, qp_, should_prune, inf_meas);
    if (is_feas == false || should_prune == true) {
      break;
    }
  }
#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "is integer = " << is_feas
    << std::endl;
#endif 

  ishigh = true;
  return is_feas;
}


void QPDProcessor::chkObjVio2_(const double *qpx, NodePtr node, double best,
                               bool &hiobjd, bool &hietavio)
{
  int err = 0;
  double d = 0.0;
  double etaval = 0.0;
  double vio = 0.0;
  double large_thresh = (node->getDepth()<5)?0.05:0.1;

  hiobjd = false;
  hietavio = false;

  d = p_->getObjective()->eval(qpx, &err);
  assert(0==err);
  if (best < INFINITY) {
    if (d > best-0.001*fabs(best)) {
      hiobjd = true;
    } else if (d > 0.999*best + 0.001*node->getLb()) {
      hiobjd = true;
    } 
    //} else if (fabs(node->getLb())>1e-3 && 
    //           d>node->getLb()+0.2*fabs(node->getLb())) {
    //  hiobjd = true;
    //}
  } else if (fabs(node->getLb())<1e-3) {
    hiobjd = (d > 0.05); 
  } else {
    hiobjd = (d > node->getLb() + 1.0*fabs(node->getLb())); 
  }

#if SPEW
  logger_->msgStream(LogDebug2) << std::setprecision(6) << me_ 
                                << " NLP obj at QPSOL = " << d
                                << " best known = " << best
                                << " node lb = " << node->getLb()
                                << " hiobjd = " << hiobjd
                                << std::endl;
#endif

  d = 0.0;
  if (eta_) {
    d = p_->getObjective()->getNonlinearFunction()->eval(qpx, &err);
    assert(0==err);
    etaval = qpx[eta_->getIndex()];
    vio = d-etaval;
    if (vio<0.0) {
      vio = 0.0;
    } else {
      vio = (fabs(etaval) > 1) ? vio/fabs(etaval) : vio;
    }
    if (vio > large_thresh) {
      hietavio = true;
    }
  }
#if SPEW
  logger_->msgStream(LogDebug2) << std::setprecision(6) << me_ 
                                << " obj nlf value at QPSOL = " << d
                                << " eta value at QPSOL = " << etaval
                                << " vio = " << vio
                                << " hietavio = " << hietavio
                                << std::endl;
#endif
}


bool QPDProcessor::isNLPFeasible2_(const double *qpx, double *vio,
                                   NodePtr node, bool &hicvio)
{
  double act, vl, vu;
  int err = 0;
  ConstConstraintPtr c;
  bool is_feas = true;
  UInt i = 0;
  UInt depth = node->getDepth();

  hicvio = false;

  for (std::vector<ConstConstraintPtr>::const_iterator it=nlCons_.begin();
       it!=nlCons_.end(); ++it, ++i) {
    c = *it;
    act = c->getActivity(qpx, &err);
    assert(0==err);
    vl = c->getLb()-act;
    vu = act-c->getUb();
    vio[i] = (vl>vu) ? (vl>0?vl:0) : (vu>0?vu:0);
    if (vio[i]>1e-5) {
      is_feas = false;
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << std::setprecision(6) 
        << "constraint " << c->getName()
        << " is violated, activity = " << act << " bounds = " << c->getLb() << " " << 
       c->getUb() << std::endl;
#endif 
    }
#if SPEW
    logger_->msgStream(LogDebug2) << me_ << "vio[" << i << "] = " << vio[i]
      << std::endl;
#endif 
    if (isLargeCVio_(c, vio[i], depth)) {
      hicvio = true;
    }
  }

#if SPEW
  logger_->msgStream(LogDebug) << me_ << "nonlinear constraints satisfied? = "
                               << is_feas << std::endl
                               << me_ << "hicvio = " << hicvio << std::endl;
#endif 
  return is_feas;
}


bool QPDProcessor::isLargeCVio_(ConstConstraintPtr c, double vio, UInt depth)
{
  double large_thresh = (depth<5)?0.05:0.1;
  double d;

  if (c->getUb()<INFINITY) {
    d = (fabs(c->getUb())>1) ? vio/fabs(c->getUb()) : vio;
  } else {
    d = (fabs(c->getLb())>1) ? vio/fabs(c->getLb()) : vio;
  }
  if (d > large_thresh) {
    return true;
  }
  return false;
}


void QPDProcessor::separateB_(ConstSolutionPtr sol, ConstSolutionPtr nlp_sol,
                              double *vio, NodePtr node, SolutionPoolPtr s_pool,
                              SeparationStatus *status)
{
  SeparationStatus s = SepaContinue;

  *status = SepaContinue;

  if (eta_) {
    separateO_(sol, nlp_sol, vio, node, s_pool, &s);
    *status = s;
    if (SepaPrune==s) {
      return;
    } 
  }

  s = SepaContinue;
  separateC_(sol, nlp_sol, vio, node, s_pool, &s);
  if (SepaPrune==s) {
    *status = SepaPrune;
  } else if (SepaResolve==s) { 
    *status = SepaResolve;
  }
}


void QPDProcessor::separateC_(ConstSolutionPtr sol, ConstSolutionPtr nlp_sol,
                              double *vio, NodePtr node, SolutionPoolPtr ,
                              SeparationStatus *status)
{
  UInt depth = node->getDepth();
  ConstConstraintPtr c;
  VariableConstIterator vbeg, vend;
  UInt i = 0;
  UInt n = qp_->getNumVars();
  LinearFunctionPtr lf;
  double lfvio;
  double minlfvio = 0.01;
  FunctionPtr f;
  double val;
  ConstraintPtr cnew;

  *status = SepaContinue;
  vbeg = qp_->varsBegin();
  vend = qp_->varsEnd();

  for (std::vector<ConstConstraintPtr>::const_iterator it=nlCons_.begin();
       it!=nlCons_.end(); ++it, ++i) {
    c = *it; 
    if (!isLargeCVio_(c, vio[i], depth)) {
      continue;
    }
    f = c->getFunction();
    getLin_(c->getFunction(), sol->getPrimal(), n, vbeg, vend, lf, val);
    lfvio = lf->eval(sol->getPrimal()) + val;
    if (c->getUb()<INFINITY) {
      lfvio = lfvio - c->getUb();
    } else {
      lfvio = c->getLb() - lfvio;
    }
    if (lfvio <= 0.0) {
      lfvio = 0.0;
    }
    if (lfvio > minlfvio) {
      f = (FunctionPtr) new Function(lf);
      cnew = qp_->newConstraint(f, c->getLb()-val, c->getUb()-val,
                                c->getName());
      ++stats_.cuts;
      *status = SepaResolve;
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "new ECP inequality: ";
      cnew->write(logger_->msgStream(LogDebug2));
      logger_->msgStream(LogDebug2) << me_ << "lfvio = " << lfvio << std::endl;
#endif 
    } else {
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "ECP linear inequality not "
        << "violated. activity = " << lf->eval(sol->getPrimal()) << " " 
        << c->getLb()-val << " " << c->getUb()-val << std::endl;
#endif 
    }
  }

  if (!nlp_sol) {
    return;
  }

  i = 0;
  for (std::vector<ConstConstraintPtr>::const_iterator it=nlCons_.begin();
       it!=nlCons_.end(); ++it, ++i) {
    c = *it; 
    if (vio[i]<1e-3) {
      continue;
    }
    f = c->getFunction();
    getLin_(c->getFunction(), nlp_sol->getPrimal(), n, vbeg, vend, lf, val);
    lfvio = lf->eval(sol->getPrimal()) + val;
    if (c->getUb()<INFINITY) {
      lfvio = lfvio - c->getUb();
    } else {
      lfvio = c->getLb() - lfvio;
    }
    if (lfvio <= 0.0) {
      lfvio = 0.0;
    }
    if (lfvio > minlfvio) {
      f = (FunctionPtr) new Function(lf);
      cnew = qp_->newConstraint(f, c->getLb()-val, c->getUb()-val,
                                c->getName());
      ++stats_.cuts;
      *status = SepaResolve;
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "new OA inequality: ";
      cnew->write(logger_->msgStream(LogDebug2));
      logger_->msgStream(LogDebug2) << me_ << "lfvio = " << lfvio << std::endl;
#endif 
    } else {
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "OA linear inequality not "
        << "violated. activity = " << lf->eval(sol->getPrimal()) << " " 
        << c->getLb()-val << " " << c->getUb()-val << std::endl;
#endif 
    }
  }
}


void QPDProcessor::separateO_(ConstSolutionPtr sol, ConstSolutionPtr nlp_sol,
                              double *, NodePtr , SolutionPoolPtr ,
                              SeparationStatus *status)
{
  double val;
  FunctionPtr f;
  UInt n = qp_->getNumVars();
  VariableConstIterator vbeg, vend;
  LinearFunctionPtr lf;
  double lfvio;
  ConstraintPtr cnew;
  double minlfvio = 0.1;

  *status = SepaContinue;
  vbeg = qp_->varsBegin();
  vend = qp_->varsEnd();

  assert(eta_);

  getObjLin_(p_->getObjective()->getFunction()->getNonlinearFunction(),
             sol->getPrimal(), n, vbeg, vend, lf, val);
  lfvio = lf->eval(sol->getPrimal()) + val;
  if (lfvio > minlfvio) {
    f = (FunctionPtr) new Function(lf);
    cnew = qp_->newConstraint(f, -INFINITY, -val, "obj_lnrztn");
    ++stats_.cuts;
    *status = SepaResolve;
#if SPEW
    logger_->msgStream(LogDebug2) << me_ << "new ECP obj inequality: ";
    cnew->write(logger_->msgStream(LogDebug2));
    logger_->msgStream(LogDebug2) << me_ << "lfvio = " << lfvio << std::endl;
#endif 
  } else {
#if SPEW
    logger_->msgStream(LogDebug2) << me_ << "obj ECP linear inequality not violated"
      << " activity = " << lf->eval(sol->getPrimal())
      << " val = " << val << std::endl;
#endif 
  }
  if (!nlp_sol) {
    return;
  }


  getObjLin_(p_->getObjective()->getFunction()->getNonlinearFunction(),
             nlp_sol->getPrimal(), n, vbeg, vend, lf, val);
  lfvio = lf->eval(sol->getPrimal()) + val;
  if (lfvio > minlfvio) {
    f = (FunctionPtr) new Function(lf);
    cnew = qp_->newConstraint(f, -INFINITY, -val, "obj_lnrztn");
    ++stats_.cuts;
    *status = SepaResolve;
#if SPEW
    logger_->msgStream(LogDebug2) << me_ << "new OA obj inequality: ";
    cnew->write(logger_->msgStream(LogDebug2));
    logger_->msgStream(LogDebug2) << me_ << "lfvio = " << lfvio << std::endl;
#endif 
  } else {
#if SPEW
    logger_->msgStream(LogDebug2) << me_ << "obj OA linear inequality not violated"
      << " activity = " << lf->eval(sol->getPrimal())
      << " val = " << val << std::endl;
#endif 
  }
}



void QPDProcessor::writeStats(std::ostream &out) const
{
  out << me_ << "nodes processed = " << stats_.proc << std::endl 
      << me_ << "nodes branched = " << stats_.bra << std::endl 
      << me_ << "nodes infeasible = " << stats_.inf << std::endl 
      << me_ << "nodes optimal = " << stats_.opt << std::endl 
      << me_ << "nodes hit ub = " << stats_.ub << std::endl 
      << me_ << "nodes with problems = " << stats_.prob << std::endl 
      << me_ << "nlps solved = " << stats_.nlp << std::endl 
      << me_ << "nlps infeasible = " << stats_.nlpI << std::endl
      << me_ << "cuts added = " << stats_.cuts << std::endl
      << me_ << "times cuts added = " << stats_.sep << std::endl
      ;
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
