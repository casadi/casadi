//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2011 - 2014 The MINOTAUR Team.
//

/**
 * \file qpOASESEngine.cpp
 * \brief Define the class qpOASESEngine.
 * \author Christian Kirches, Interdisciplinary Center for Scientific
 *         Computing (IWR), Heidelberg University, GERMANY
 *
 * Implement the interface to qpOASES engine that provides
 * parametric hot-starting.
 */

#ifndef USE_QPOASES
	#error Need to set USE_QPOASES
#endif

#include <cmath>
#include <iomanip>
#include <iostream>

#include <boost/scoped_array.hpp>

#include "MinotaurConfig.h"
#include "qpOASESEngine.h"
#include "Constraint.h"
#include "Environment.h"
#include "Function.h"
#include "HessianOfLag.h"
#include "LinearFunction.h"
#include "Logger.h"
#include "Objective.h"
#include "Option.h"
#include "QuadraticFunction.h"
#include "Problem.h"
#include "ProblemSize.h"
#include "Solution.h"
#include "Timer.h"
#include "Variable.h"

using namespace Minotaur;

/*
 * qpOASES data structure default constructor
 */
qpOASESData::qpOASESData (int n_, int m_, int nnzh_, int nnza_)
: n (n_),
  m (m_),
  g (new double [n_]),
  lb (new double [n_]),
  ub (new double [n_]),
  lbA (new double [m_]),
  ubA (new double [m_])
{
  qpOASES::Options qp_opt;

  qp   = new qpOASES::SQProblem (n_, m_);

  // set initial options
  qp_opt = qp->getOptions ();
  /*
     qp_opt.enableEqualities  = qpOASES::BT_TRUE;
     qp_opt.enableFullLITests = qpOASES::BT_TRUE;
     qp_opt.enableNZCTests    = qpOASES::BT_TRUE;
     qp_opt.printLevel        = qpOASES::PL_MEDIUM;
     qp_opt.printLevel        = qpOASES::PL_NONE;
     qp_opt.epsFlipping       = 1.0e-12;
     qp_opt.enableCholeskyRefactorisation = 0;
     qp_opt.initialFarBounds  = 1.0e+4;
     */
  qp_opt.setToDefault ();
  qp_opt.enableFullLITests             =  qpOASES::BT_FALSE;
  qp_opt.enableNZCTests                = qpOASES::BT_FALSE;
  qp_opt.enableCholeskyRefactorisation = 1;
  qp_opt.numRefinementSteps            = 0;
  qp_opt.enableEqualities              = qpOASES::BT_FALSE;
  qp_opt.initialFarBounds              = 1.0E+20;
  qp_opt.epsLITests                    = 1.0e8 * 1e-16;
  //	qp_opt.enableFlippingBounds  = qpOASES::BT_FALSE;
  qp_opt.enableFarBounds  = qpOASES::BT_FALSE;

  qp_opt.printLevel        = qpOASES::PL_HIGH;
  qp_opt.printLevel        = qpOASES::PL_TABULAR;
  qp_opt.printLevel        = qpOASES::PL_MEDIUM;
  qp_opt.printLevel        = qpOASES::PL_NONE;
  qp->setOptions (qp_opt);

  Hir  = new qpOASES::sparse_int_t[nnzh_];
  Hjc  = new qpOASES::sparse_int_t[n+1];
  Hval = new qpOASES::real_t[nnzh_];
  H    = new qpOASES::SymSparseMat (n, n, Hir, Hjc, Hval);

  Air  = new qpOASES::sparse_int_t[m+1];
  Ajc  = new qpOASES::sparse_int_t[nnza_];
  Aval = new qpOASES::real_t[nnza_];
  A    = new qpOASES::SparseMatrixRow (m, n, Air, Ajc, Aval);
}

/*
 * qpOASES data structure default destructor
 */
qpOASESData::~qpOASESData ()
{
  if (A != 0) {
    A->free ();
    delete A;
    A = 0;
    Air = 0;
    Ajc = 0;
    Aval = 0;
  }
  if (H != 0) {
    H->free ();
    delete H;
    H = 0;
    Hir = 0;
    Hjc = 0;
    Hval = 0;
  }
  if (ubA != 0) {
    delete[] ubA;
    ubA = 0;
  }
  if (lbA != 0) {
    delete[] lbA;
    lbA = 0;
  }
  if (ub != 0) {
    delete[] ub;
    ub = 0;
  }
  if (lb != 0) {
    delete[] lb;
    lb = 0;
  }
  if (g != 0) {
    delete[] g;
    g = 0;
  }
  if (qp != 0) {
    delete qp;
    qp = 0;
  }
}



qpOASESEngine::qpOASESEngine()
: env_(EnvPtr ()),
  problem_ (ProblemPtr ()),
  sol_ (SolutionPtr ()),
  consModed_ (false),
  hessModed_ (false),
  data_ (0),
  resolveError_ (true),
  iterLimit_ (10000),
  maxIterLimit_ (10000),
  infty_ (1E+20),
  bTol_ (1E-9),
  stats_ (),
  strBr_ (false),
  ldWarm_ (false),
  n_changed_ (0)
{
  status_ = EngineError;
  logger_ = (LoggerPtr) new Logger (LogInfo);

  memset (&stats_, 0, sizeof (stats_));
}


qpOASESEngine::qpOASESEngine(EnvPtr env)
: env_ (env),
  problem_ (ProblemPtr ()),
  sol_ (SolutionPtr ()),
  consModed_ (false),
  hessModed_ (false),
  data_ (0),
  resolveError_ (true),
  iterLimit_ (10000),
  maxIterLimit_ (10000),
  infty_ (1E+20),
  bTol_ (1E-9),
  stats_ (),
  strBr_ (false),
  ldWarm_ (false),
  n_changed_ (0)
{
  status_ = EngineError;
  logger_ = (LoggerPtr) new Logger ((LogLevel) (env->getOptions ()->findInt
                                                ("bqpd_log_level")->getValue
                                                ()));

  memset (&stats_, 0, sizeof (stats_));
}


qpOASESEngine::~qpOASESEngine()
{
  if (sol_)
    sol_.reset ();
  clear ();
}


EnginePtr qpOASESEngine::emptyCopy ()
{
  if (env_)
    return (qpOASESEnginePtr) new qpOASESEngine (env_);
  else
    return (qpOASESEnginePtr) new qpOASESEngine ();
}


void qpOASESEngine::load (ProblemPtr problem)
{
  problem_ = problem;
  consModed_ = true;
  hessModed_ = true;
  problem->setEngine (this);
}


void qpOASESEngine::clear ()
{
  if (problem_) {
    problem_->unsetEngine ();
    problem_.reset ();
  }
  if (data_) {
    delete data_;
    data_ = 0;
  }
}

void qpOASESEngine::load_ ()
{
  problem_->prepareForSolve ();

  int n     = problem_->getNumVars ();
  int m     = problem_->getNumCons ();
  int nnzh  = problem_->getNumHessNnzs () * 2;
  int nnza  = problem_->getNumJacNnzs ();

  data_     = new qpOASESData (n, m, nnzh, nnza);
  sol_      = (SolutionPtr) new Solution (INFINITY, 0, problem_);

  logger_->MsgStream (LogDebug1) << "qpOASES: Loaded problem (n = " << n
    << ", m = " << m << ", nnzh = " << nnzh << ", nnza = " << nnza
    << ")." << std::endl;
}


void qpOASESEngine::setGradient_ ()
{
  qpOASES::sparse_int_t *Air  = data_->Air;
  qpOASES::sparse_int_t *Ajc  = data_->Ajc;
  qpOASES::real_t       *Aval = data_->Aval;
  size_t row, cnt;
  int error = 0;
  double *x;
  LinearFunctionPtr lf;

  // objective gradient (dense)
  ObjectivePtr oPtr = problem_->getObjective ();
  assert(oPtr);

  x = new double[data_->n];
  memset (x, 0, data_->n * sizeof (double));
  objOff_ = oPtr->eval (x, &error);
  oPtr->evalGradient (x, data_->g, &error);
  delete[] x;

  // constraints jacobian (sparse row-compressed)
  cnt = 0;
  row = 0;
  Air[row] = cnt;
  for (ConstraintConstIterator it2 = problem_->consBegin ();
       it2 != problem_->consEnd(); ++it2) {
    lf = (*it2)->getLinearFunction ();
    if (lf) {
      for (VariableGroupConstIterator it = lf->termsBegin (); it != lf->termsEnd(); ++it) {
        Ajc[cnt] = it->first->getIndex ();
        Aval[cnt] = it->second;
        ++cnt;
      }
    }
    ++row;
    Air[row] = cnt;
  }
}


void qpOASESEngine::setHessian_ ()
{
  UInt n, m, hess_nnzs;
  UInt *iRow, *jCol;
  double *values, *con_mult, *x;
  double obj_mult = 1.0;
  int error = 0;

  n         = problem_->getNumVars();
  m         = problem_->getNumCons();
  hess_nnzs = problem_->getHessian ()->getNumNz ();

  iRow      = new UInt[hess_nnzs];
  jCol      = new UInt[hess_nnzs];
  values    = new double[hess_nnzs];
  con_mult  = new double[m];
  x         = new double[n];

  // as for Hessian in a triplets format, the first index vector passed
  // to fillRowColIndices will be sorted in ascending order
  problem_->getHessian ()->fillRowColIndices (jCol, iRow);
  std::fill (x, x+n, 0);
  std::fill (con_mult, con_mult+m, 0);
  problem_->getHessian()->fillRowColValues (x, obj_mult, con_mult, values, &error);
  assert (0 == error);

  // covert it to qpOASES column compressed format. make sure to write _both_ triangles
  qpOASES::sparse_int_t *Hir  = data_->Hir;
  qpOASES::sparse_int_t *Hjc  = data_->Hjc;
  qpOASES::real_t       *Hval = data_->Hval;
  size_t col, cnt;

  cnt = 0;
  col = 0;
  Hjc[col] = cnt;

  for (UInt ii = 0; ii < hess_nnzs; ++ii) {
    while (jCol[ii] > col) {
      // upper tri part of column ended, now complete the part
      // that sits in the lower triangle, qpOASES needs it :(
      for (UInt jj = ii; jj < hess_nnzs; ++jj)
        if (iRow[jj] == col) {
          Hir[cnt] = jCol[jj];
          Hval[cnt] = values[jj];
          ++cnt;
        }
      // move on to next column
      ++col;
      Hjc[col] = cnt;
    }
    // make a column entry that sits in the upper triangle
    Hir[cnt] = iRow[ii];
    Hval[cnt] = values[ii];
    ++cnt;
  }
  // fill remaining empty columns
  while (col < n) {
    ++col;
    Hjc[col] = cnt;
  }

  data_->H->createDiagInfo();

  // aufraeumen
  delete [] x;
  delete [] con_mult;
  delete [] values;
  delete [] jCol;
  delete [] iRow;
}


EngineStatus qpOASESEngine::solve ()
{
  double f;
  int mode;

  if (consModed_) {
    delete data_;
    data_ = 0;
    load_ ();
    setGradient_ ();
    setHessian_ ();
    setVarBounds_();
    setConsBounds_();
    mode = 0;
  } else {
    mode = 2;
  }

  consModed_ = false;
  hessModed_ = false;

  solve_ (mode, &f);

  if (EngineError == status_ && resolveError_) {
    logger_->MsgStream (LogInfo)
      << "qpOASES: failed to solve in mode " << mode << std::endl;
    mode = 0;
    logger_->MsgStream(LogInfo)
      << "qpOASES: now solving in mode " << mode << std::endl;
    solve_ (mode, &f);
    if (EngineError == status_) {
      logger_->MsgStream(LogInfo)
        << "qpOASES: catastrophic failure in mode" << mode << " as well." << std::endl;
    }
  }

  // store the solution
  storeSol_ (f);

  return status_;
}


void qpOASESEngine::solve_ (int mode, double *f)
{
  qpOASES::returnValue rv = qpOASES::SUCCESSFUL_RETURN;

  iterLimit_ = strBr_ ? env_->getOptions ()->findInt ("strbr_pivot_limit")->getValue () : 10000;
  int    nWSR    = iterLimit_;
  double cputime = 1E+20;

  status_ = EngineUnknownStatus;
  *f      = 1E+20;


  switch (mode) {
  case 0:
  default:
    rv = data_->qp->init (data_->H, data_->g, data_->A, data_->lb, data_->ub,
                          data_->lbA, data_->ubA, nWSR, &cputime);
    ++stats_.n_cold;
    stats_.piv_cold += nWSR;
    stats_.time_cold += cputime;
    break;

  case 2:
    bool uws = env_->getOptions ()->findBool ("use_warmstart")->getValue ();
    if (ldWarm_ != 0 && uws)
      rv = data_->qp->init (data_->H, data_->g, data_->A, data_->lb, data_->ub,
                            data_->lbA, data_->ubA, nWSR, &cputime,
                            ws.pt.x, ws.pt.y, &ws.pt.bounds, &ws.pt.constraints);
    else 
      rv = data_->qp->hotstart (data_->g, data_->lb, data_->ub,
                                data_->lbA, data_->ubA, nWSR, &cputime);
    /*		else
                rv = data_->qp->init (data_->H, data_->g, data_->A, data_->lb, data_->ub,
                data_->lbA, data_->ubA, nWSR, &cputime);
                */
    if (ldWarm_) {
      ++stats_.n_warm;
      stats_.piv_warm += nWSR;
      stats_.time_warm += cputime;
    }
    else if (strBr_) {
      ++stats_.n_strbr;
      stats_.piv_strbr += nWSR;
      stats_.time_strbr += cputime;
    }
    else {
      ++stats_.n_dive;
      stats_.piv_dive += nWSR;
      stats_.time_dive += cputime;
    }
    break;
  }

  ++stats_.n_total;
  stats_.piv_total += nWSR;
  stats_.time_total += cputime;


  if (strBr_) {
    stats_.strCalls += 1;
    stats_.strTime  += cputime;
    stats_.strIters += nWSR;
  }
  stats_.calls += 1;
  stats_.time  += cputime;
  stats_.iters += nWSR;

  *f = data_->qp->getObjVal ();

  const char *rvmsg;

  switch (rv) {
  case qpOASES::SUCCESSFUL_RETURN:
    status_ = ProvenLocalOptimal;
    rvmsg = "optimal";
    break;

  case qpOASES::RET_MAX_NWSR_REACHED:
    status_ = EngineIterationLimit;
    rvmsg = "maxiter";
    break;

  case qpOASES::RET_HOTSTART_STOPPED_INFEASIBILITY:
  case qpOASES::RET_INIT_FAILED_INFEASIBILITY:
  case qpOASES::RET_QP_INFEASIBLE:
    status_ = ProvenLocalInfeasible;
    rvmsg = "infsble";
    break;

  case qpOASES::RET_HOTSTART_STOPPED_UNBOUNDEDNESS:
  case qpOASES::RET_INIT_FAILED_UNBOUNDEDNESS:
  case qpOASES::RET_QP_UNBOUNDED:
    status_ = ProvenUnbounded;
    rvmsg = "unbnded";
    break;

  case qpOASES::RET_INIT_FAILED:
  case qpOASES::RET_INIT_FAILED_HOTSTART:
  case qpOASES::RET_HOTSTART_FAILED:
  case qpOASES::RET_INIT_FAILED_TQ:
  case qpOASES::RET_STEPDIRECTION_FAILED_TQ:
  case qpOASES::RET_ENSURELI_FAILED_TQ:
  case qpOASES::RET_INIT_FAILED_CHOLESKY:
  case qpOASES::RET_STEPDIRECTION_FAILED_CHOLESKY:
    status_ = EngineError;
    rvmsg = "LAbrken";
    break;

  default:
    fprintf (stderr, "CAUGHT ERROR CODE %d IN qpOASESEngine.cpp!\n", rv);
    status_ = EngineError;
    rvmsg = "fatal";
    break;
  }

  /*
     if (mode==0)
     action = "cold  ";
     else
     if (ldWarm_)
     action = "warm  ";
     else
     if (strBr_)
     action = "strbr ";
     else
     action = "dive  ";


     fprintf (stderr, "%s  %3d  %s  %3d\n", action, n_changed_, rvmsg, nWSR);
     */
  ldWarm_=false;
  n_changed_=0;

#ifdef DEBUG
  logger_->MsgStream(LogDebug)
    << "qpOASES: returnvalue = " << rv
    << " message = " << qpOASES::MessageHandling::getErrorCodeMessage (rv)
    << "status = " << getStatusString () << std::endl;
  logger_->MsgStream(LogDebug)
    << "qpOASES: objective = " << f << std::endl;
#endif
}


void qpOASESEngine::storeSol_ (double f)
{
  UInt n = data_->n;

  sol_->setObjValue (f + objOff_);
  sol_->setPrimal (data_->qp->x);
  sol_->setDualOfVars (data_->qp->y);
  sol_->setDualOfCons (data_->qp->y + n);
}


void qpOASESEngine::setInitialPoint_ ()
{
  // qpOASES keeps an internal copy of the previous optimal solution.
  // This should suffice.
}


void qpOASESEngine::setVarBounds_ ()
{
  VariablePtr vPtr;
  VariableConstIterator vIter;
  double l, u;
  double *bl = data_->lb;
  double *bu = data_->ub;

  for (vIter = problem_->varsBegin (); vIter != problem_->varsEnd (); ++vIter, ++bl, ++bu) {
    vPtr = *vIter;
    l = std::max (-infty_, vPtr->getLb ());
    u = std::min ( infty_, vPtr->getUb ());
    *bl = l;
    *bu = u;
  }
}


void qpOASESEngine::setConsBounds_ ()
{
  ConstraintPtr cPtr;
  ConstraintConstIterator cIter;
  double *bl = data_->lbA;
  double *bu = data_->ubA;
  double l, u;

  for (cIter = problem_->consBegin (); cIter != problem_->consEnd (); ++cIter, ++bl, ++bu) {
    cPtr = *cIter;
    l = std::max (-infty_, cPtr->getLb ());
    u = std::min ( infty_, cPtr->getUb ());
    *bl = l;
    *bu = u;
  }
}


double qpOASESEngine::getSolutionValue ()
{
  return (sol_) ? sol_->getObjValue() : INFINITY;  // throw exception instead
}


ConstSolutionPtr qpOASESEngine::getSolution ()
{
  return (sol_) ? sol_ : SolutionPtr (); // NULL
}


EngineStatus qpOASESEngine::getStatus ()
{
  return status_;
}


void qpOASESEngine::changeBound (ConstraintPtr, BoundType, double)
{
  // TODO: This will trigger mode 0, is that necessary?
  consModed_ = true;
  fprintf(stderr,"changeBound 1\n");

  assert(!"implement me");
}


void qpOASESEngine::changeBound (VariablePtr vPtr, BoundType lu, double new_val)
{
  if (data_ == 0) return;

  UInt index = vPtr->getIndex ();
  if (lu == Lower) {
    double old = data_->lb[index];
    data_->lb[index] = std::max (-infty_, new_val);
    if (old != data_->lb[index]) ++n_changed_;
  }
  if (lu == Upper) {
    double old = data_->ub[index];
    data_->ub[index] = std::min ( infty_, new_val);
    if (old != data_->ub[index]) ++n_changed_;
  }
}


void qpOASESEngine::changeBound (VariablePtr, double, double)
{
  // empty
  fprintf(stderr,"changeBound 3\n");

  assert(!"implement me");
}


void qpOASESEngine::changeObj (FunctionPtr, double)
{
  // TODO: This will trigger mode 0, is that necessary?
  consModed_ = true;
}


void qpOASESEngine::negateObj ()
{
  // TODO: This will trigger mode 0, is that necessary?
  consModed_ = true;
}


void qpOASESEngine::changeConstraint (ConstraintPtr, LinearFunctionPtr,
                                 const double & )
{
  // TODO: This will trigger mode 0, is that necessary?
  consModed_ = true;
}


void qpOASESEngine::enableStrBrSetup ()
{
  strBr_ = true;
}


void qpOASESEngine::disableStrBrSetup()
{
  strBr_ = false;
}


void qpOASESEngine::setIterationLimit(int)
{
  //	iterLimit_ = (limit >= 1) ? limit : maxIterLimit_;
}


void qpOASESEngine::resetIterationLimit ()
{
  iterLimit_ = maxIterLimit_;
}


void qpOASESEngine::writeStats ()
{
  std::string me = "qpOASES:  ";
  logger_->MsgStream(LogInfo)
    << me << "total calls            = " << stats_.calls << std::endl
    << me << "strong branching calls = " << stats_.strCalls << std::endl
    << me << "total time in solving  = " << stats_.time  << std::endl
    << me << "time in str branching  = " << stats_.strTime << std::endl
    << me << "time in copying data   = " << stats_.cTime << std::endl
    << me << "total iterations       = " << stats_.iters << std::endl
    << me << "strong br iterations   = " << stats_.strIters << std::endl;

  fprintf (stderr, "\nQP mode         count   pivots       time\n");
  fprintf (stderr, "QP cold start  %6d %8d %10.3e\n", stats_.n_cold, stats_.piv_cold, stats_.time_cold);
  fprintf (stderr, "QP str. br.    %6d %8d %10.3e\n", stats_.n_strbr, stats_.piv_strbr, stats_.time_strbr);
  fprintf (stderr, "QP warm start  %6d %8d %10.3e\n", stats_.n_warm, stats_.piv_warm, stats_.time_warm);
  fprintf (stderr, "QP diving      %6d %8d %10.3e\n", stats_.n_dive, stats_.piv_dive, stats_.time_dive);
  fprintf (stderr, "QP total       %6d %8d %10.3e\n\n", stats_.n_total, stats_.piv_total, stats_.time_total);
}


std::string qpOASESEngine::getName() const
{
  return "qpOASES";
}

ConstWarmStartPtr qpOASESEngine::getWarmStart ()
{
  qpOASESWarmStart ws;
  return (WarmStartPtr) &ws;
}

WarmStartPtr qpOASESEngine::getWarmStartCopy ()
{
  qpOASESWarmStart *ws = new qpOASESWarmStart ();
  savePoint (&ws->pt);
  return (WarmStartPtr) ws;
}

void qpOASESEngine::loadFromWarmStart(WarmStartPtr wp)
{
  ws = * static_cast<qpOASESWarmStart*>(wp.get());
  ldWarm_ = true;
}

////////////////////////////////////////////////////////////////////////

void qpOASESEngine::savePoint (qpOASESPoint *pt)
{	
  pt->x   = new double[data_->n];
  pt->y   = new double[data_->n + data_->m];
  pt->lb  = new double[data_->n];
  pt->ub  = new double[data_->n];
  pt->lbA = new double[data_->m];
  pt->ubA = new double[data_->m];

  data_->qp->getBounds (pt->bounds);
  data_->qp->getConstraints (pt->constraints);

  memcpy (pt->x, data_->qp->x, sizeof (double) * data_->n);
  memcpy (pt->y, data_->qp->y, sizeof (double) * (data_->n + data_->m));
  memcpy (pt->lb, data_->qp->lb, sizeof (double) * data_->n);
  memcpy (pt->ub, data_->qp->ub, sizeof (double) * data_->n);
  memcpy (pt->lbA, data_->qp->lbA, sizeof (double) * data_->m);
  memcpy (pt->ubA, data_->qp->ubA, sizeof (double) * data_->m);
}

void qpOASESEngine::restorePoint (qpOASESPoint *pt)
{
  data_->qp->bounds = pt->bounds;
  data_->qp->constraints = pt->constraints;
  memcpy (data_->qp->x, pt->x, sizeof (double) * data_->n);
  memcpy (data_->qp->y, pt->y, sizeof (double) * (data_->n + data_->m));
  memcpy (data_->qp->lb, pt->lb, sizeof (double) * data_->n);
  memcpy (data_->qp->ub, pt->ub, sizeof (double) * data_->n);
  memcpy (data_->qp->lbA, pt->lbA, sizeof (double) * data_->m);
  memcpy (data_->qp->ubA, pt->ubA, sizeof (double) * data_->m);

  data_->qp->A->times(1, 1.0, data_->qp->x, data_->n, 0.0, data_->qp->Ax, data_->m);
  for (UInt ii=0; ii < data_->m; ++ii) {
    data_->qp->Ax_l[ii] = data_->qp->Ax[ii] - pt->lbA[ii];
    data_->qp->Ax_u[ii] = pt->ubA[ii] - data_->qp->Ax[ii];
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
