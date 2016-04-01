// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

/**
 * \file FilterSQPEngine.cpp
 * \brief Define the class FilterSQPEngine.
 * \author Sven Leyffer, Argonne National Laboratory
 * 
 * Implement the interface to filterSQP engine and provide warm-starting
 * capability.
 */

#include <cmath>
#include <iostream>
#include <iomanip>

#include "MinotaurConfig.h"
#include "Constraint.h"
#include "Environment.h"
#include "FilterSQPEngine.h"
#include "FilterSQPEngineTypes.h"
#include "HessianOfLag.h"
#include "Jacobian.h"
#include "Logger.h"
#include "Objective.h"
#include "Option.h"
#include "Problem.h"
#include "Solution.h"
#include "Timer.h"
#include "Variable.h"

//#define SPEW 1

using namespace Minotaur;

const std::string FilterSQPEngine::me_ = "FilterSQPEngine: ";

int * convertPtrToInt(uintptr_t u)
{
  const int MAX_SIZE_PER_INT=4;
  const uintptr_t AND_MASK=0x7FFF;
  const int SHIFT_WIDTH=15;

  // *2.0 because its unsigned. +1 because first element stores the number of
  // integer-blocks.
  int size = (int) ceil(((double)sizeof(uintptr_t)/MAX_SIZE_PER_INT*2.0))+1; 
  int *iarray = new int[size]; 
  int r = 1;
  uintptr_t lsbs; // # SHIFT_WIDTH least significant bits.

  while (u>0) {
    lsbs = (u & AND_MASK);
    u = u>>SHIFT_WIDTH;
    iarray[r] = (int) lsbs;
    ++r;
  }
  iarray[0] = r;
  return iarray;
}


uintptr_t convertIntToPtr(int *iarray)
{
  const int SHIFT_WIDTH=15;
  uintptr_t u = 0;
  int size = iarray[0];
  for (int r=size-1; r>0; --r) {
    u = u<<SHIFT_WIDTH;
    u = (u | (uintptr_t) iarray[r]);
  }
  return u;
}


// objective function evaluation. This routine is required by
// filtersqp-library. We can not put this inside a namespace.
void objfun(real *x, fint *, real *f, real *, fint *iuser, 
    fint *errflag) 
{
  // first get pointer to interface from iuser. 
  assert (iuser);
  uintptr_t u = convertIntToPtr(iuser);
  FilterSQPEngine *f_engine = (FilterSQPEngine *) u;
  assert(f);
  assert (f_engine);

  // call the interface to get objective value
  *f = f_engine->evalObjValue(x, errflag);
}


void objgrad (fint *, fint *, fint *, real *, real *, fint *, fint
    *, real *, fint *, fint *)
{
  assert(!"implement me!");
}


// declaration of constraint functions evaluation
void confun(real *x, fint *, fint *, real *c, real *, fint *,
    real *, fint *iuser, fint *errflag)
{
  // first get pointer to interface from iuser. 
  uintptr_t u = convertIntToPtr(iuser);
  FilterSQPEngine *f_engine = (FilterSQPEngine *) u;

  // call the interface to evaluate constraints
  f_engine->evalCons(x, c, errflag);
}


// evaluation of gradients (objective & constraints)
void gradient(fint *, fint *, fint *, real *x, real *a, fint *,
    fint *, real *, fint *iuser, fint *errflag)
{
  // first get pointer to interface from iuser. 
  uintptr_t u = convertIntToPtr(iuser);
  FilterSQPEngine *f_engine = (FilterSQPEngine *) u;

  // call the interface to evaluate gradients.
  // user, M, mxa, la, maxa are not used?
  f_engine->evalGrad(x, a, errflag);
}


void hessian(real *x, fint *, fint *, fint *phase, real *lam,
    real *ws, fint *lws, real *, fint *iuser,
    fint *l_hess, fint *li_hess, fint *errflag)
{
  // first get pointer to interface from iuser. 
  uintptr_t u = convertIntToPtr(iuser);
  FilterSQPEngine *f_engine = (FilterSQPEngine *) u;

  // call the interface to evaluate hessian
  f_engine->evalHessian(x, lam, *phase, ws, lws, l_hess, li_hess, errflag);
}


FilterSQPWarmStart::FilterSQPWarmStart()
  : sol_(SolutionPtr()) // NULL
{
}


FilterSQPWarmStart::~FilterSQPWarmStart()
{
  sol_.reset();
}


// copy
FilterSQPWarmStart::FilterSQPWarmStart(ConstFilterWSPtr warmSt)
{
  if (warmSt && warmSt->sol_) {
    sol_ = (SolutionPtr) new Solution(warmSt->sol_);
  } else {
    sol_ = SolutionPtr(); // NULL
  }
}
 

SolutionPtr FilterSQPWarmStart::getPoint()
{
  return sol_;
}


bool FilterSQPWarmStart::hasInfo()
{
  if (sol_ && sol_->getPrimal()) {
    return true;
  }
  return false;
}


void FilterSQPWarmStart::setPoint(SolutionPtr sol)
{
  sol_ = sol;
}


void FilterSQPWarmStart::write(std::ostream &out) const
{
  out << "FilterSQP warm start information:"  << std::endl;
  sol_->write(out);
}


FilterSQPEngine::FilterSQPEngine()
: a_(0),
  bl_(0),
  bTol_(1e-9),
  bu_(0),
  c_(0),
  consChanged_(true),
  cstype_(0),
  env_(EnvPtr()),
  feasTol_(1e-6),
  istat_(0),
  la_(0),                 // NULL
  lam_(0),
  lws_(0),                // NULL
  lws2_(0),
  maxIterLimit_(1000),
  mlam_(0),
  prepareWs_(false),
  rstat_(0),
  s_(0),
  saveSol_(true),
  sol_(SolutionPtr()),    // NULL
  stats_(0),
  strBr_(false),
  timer_(0),
  useWs_(false),
  warmSt_(FilterWSPtr()),
  ws_(0),
  x_(0)
{
#ifndef USE_FILTERSQP 
#error Need to set USE_FILTERSQP
#endif
  problem_ = ProblemPtr(); // NULLstatus_ = EngineError;
  logger_ = (LoggerPtr) new Logger(LogInfo);
  iterLimit_ = maxIterLimit_;
}


FilterSQPEngine::FilterSQPEngine(EnvPtr env)
: a_(0),
  bl_(0),
  bTol_(1e-9),
  bu_(0),
  c_(0),
  consChanged_(true),
  cstype_(0),
  env_(env),
  feasTol_(1e-6),
  istat_(0),
  la_(0),
  lam_(0),
  lws_(0),
  lws2_(0),
  maxIterLimit_(1000),
  mlam_(0),
  rstat_(0),
  s_(0),
  saveSol_(true),
  sol_(SolutionPtr()),
  strBr_(false),
  warmSt_(FilterWSPtr()),
  ws_(0),
  x_(0)
{
  problem_ = ProblemPtr(); // NULL
  status_ = EngineUnknownStatus;
  logger_ = (LoggerPtr) new Logger((LogLevel) (env->getOptions()
        ->findInt("engine_log_level")->getValue()));
  iterLimit_ = maxIterLimit_;
  if (env->getOptions()->findBool("use_warmstart")->getValue()==true) {
    prepareWs_ = true;
    useWs_ = true;
  } else {
    prepareWs_ = false;
    useWs_ = false;
  }

  timer_ = env->getNewTimer();

  stats_ = new FilterSQPStats();
  stats_->calls    = 0;
  stats_->strCalls = 0;
  stats_->time     = 0;
  stats_->strTime  = 0;
  stats_->iters    = 0;
  stats_->strIters = 0;
}


FilterSQPEngine::~FilterSQPEngine()
{
  //delete ;
  if (c_) {
    freeStorage_();
    c_ = 0;
  }
  if (sol_) {
    sol_.reset();
  }
  if (timer_) {
    delete timer_;
  }
  if (stats_) {
    delete stats_;
  }
  if (problem_) {
    problem_->unsetEngine();
    problem_.reset();
  }
}


void FilterSQPEngine::addConstraint(ConstraintPtr)
{
  consChanged_ = true;
}


void FilterSQPEngine::changeBound(ConstraintPtr cons, BoundType lu, 
                                  double new_val)
{
  if (bl_) {
    if (Lower == lu) {
      bl_[problem_->getNumVars()+cons->getIndex()] = new_val;
    } else {
      bu_[problem_->getNumVars()+cons->getIndex()] = new_val;
    }
  }
}


void FilterSQPEngine::changeBound(VariablePtr v, BoundType lu, double val)
{
  if (bl_) {
    if (Lower == lu) {
      bl_[v->getIndex()] = val;
    } else {
      bu_[v->getIndex()] = val;
    }
  }
}


void FilterSQPEngine::changeBound(VariablePtr v, double lb, double ub)
{
  if (bl_) {
    bl_[v->getIndex()] = lb;
    bu_[v->getIndex()] = ub;
  }
}


void FilterSQPEngine::changeConstraint(ConstraintPtr, LinearFunctionPtr, 
                                       double , double)
{
  // no need to do anything because the 'solve' function reloads constraints
  // from problem.
  consChanged_ = true;
}


void FilterSQPEngine::changeConstraint(ConstraintPtr, NonlinearFunctionPtr)
{
  consChanged_ = true;
}


void FilterSQPEngine::changeObj(FunctionPtr, double)
{
  consChanged_ = true;
}


void FilterSQPEngine::clear() 
{
  if (c_) {
    freeStorage_();
    c_ = 0;
  }
  if (problem_) {
    problem_->unsetEngine();
    problem_.reset();
  }
  if (sol_) {
    sol_.reset();
  }
  if (warmSt_) {
    warmSt_.reset();
  }
  strBr_       = false;
  consChanged_ = true;
}


void FilterSQPEngine::disableStrBrSetup()
{
  saveSol_ = true;
  strBr_   = false;
}


EnginePtr FilterSQPEngine::emptyCopy()
{
  if (env_) {
    return (FilterSQPEnginePtr) new FilterSQPEngine(env_);
  }
  return (FilterSQPEnginePtr) new FilterSQPEngine();
}


void FilterSQPEngine::enableStrBrSetup()
{
  saveSol_ = false;
  strBr_   = true;
}


void FilterSQPEngine::evalCons(const double *x, double *c, int *error) 
{
  ConstraintConstIterator cIter;
  ConstraintPtr cPtr;
  UInt i=0;
  int e;
  *error = 0;
  //problem_->write(std::cout);
  for (cIter=problem_->consBegin(); cIter!=problem_->consEnd(); ++cIter) {
    e = 0;
    cPtr = *cIter;
    c[i] = cPtr->getActivity(x, &e);
    if (e!=0) {
      *error = 1;
    }
    ++i;
  }
#if SPEW
  if (logger_->getMaxLevel() > LogDebug) {
    VariableConstIterator vIter;
    i=0;
    logger_->msgStream(LogDebug2) << me_ << std::endl;
    for (vIter=problem_->varsBegin(); vIter!=problem_->varsEnd(); ++vIter) {
      logger_->msgStream(LogDebug2) << (*vIter)->getName() 
        << " = " << x[i] << std::endl;
      ++i;
    }
    i=0;
    for (cIter=problem_->consBegin(); cIter!=problem_->consEnd(); ++cIter) {
      logger_->msgStream(LogDebug2) << (*cIter)->getName() 
        << " = " << c[i] << std::endl;
      ++i;
    }
  }
  logger_->msgStream(LogDebug2) << me_ << "error in evalCons = " 
    << *error << std::endl;
#endif
}


void FilterSQPEngine::evalGrad(const double *x, double *a, int *error) 
{
  UInt n         = problem_->getNumVars();
  UInt jac_nnz   = problem_->getJacobian()->getNumNz();
  double *values = NULL;
  ObjectivePtr o = problem_->getObjective();
  int e2         = 0;

  *error = 0;
  // first zero out all values of 'a'
  std::fill(a, a+n+jac_nnz, 0);


  // compute the gradient of the objective function f(x)
  // dense, it fills up a[0,1,2 ..., (n-1)]
  if (o) {
    o->evalGradient(x, a, &e2);
  }
#if SPEW
  if (logger_->getMaxLevel() > LogDebug) {
    logger_->msgStream(LogDebug2) << me_ << "obj gradient" << std::endl;
    for (UInt i=0; i<n; ++i) {
      logger_->msgStream(LogDebug2) << "  f'[" << i << "] = " << a[i] 
        << std::endl;
    }
  }
#endif

  // Gradient of each constraint function sits together in 'a'. First get the
  // starting position where such a gradient should go.
  values = a+n;
  problem_->getJacobian()->fillRowColValues(x, values, error);
  if (e2>0) {
    *error = e2;
  }
#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "error in evalGrad = " 
    << *error << std::endl;
  if (logger_->getMaxLevel() > LogDebug2) {
    logger_->msgStream(LogDebug2) << me_ << "obj gradient" << std::endl;
    for (UInt i=0; i<n+jac_nnz; ++i) {
      logger_->msgStream(LogDebug2) << std::setprecision(15)
        << "  jac[" << i << "] = " << a[i] << std::endl;
    }
  }
#endif
}


void FilterSQPEngine::evalHessian(const double *x, double *lam, 
                                  const int phase, double *ws, int *lws, 
                                  int *l_hess, int *li_hess, int *error)
{
  double obj_mult = 0;
  double *values = ws;
  *error = 0;

  // don't know why im using these values.
  *l_hess  = problem_->getHessian()->getNumNz();
  *li_hess = 3+problem_->getHessian()->getNumNz()+problem_->getNumVars();

  if (phase==2) {
    obj_mult = 1.0;
  }

  // first 'n' multipliers are for variable bounds. skip them.
  lam += problem_->getNumVars();
  // need to negate the rest of the multipliers. I don't know why.
  for (UInt i=0; i<problem_->getNumCons(); ++i) {
    mlam_[i] = -1.0*lam[i];
  }
  std::copy(lws_, lws_+*li_hess, lws);
  problem_->getHessian()->fillRowColValues(x, obj_mult, mlam_, values, error);
#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "l_hess = " << *l_hess 
    << std::endl;
  logger_->msgStream(LogDebug2) << me_ << "obj lam = " << obj_mult
    << std::endl;
  for (UInt i=0; i<problem_->getNumCons(); ++i) {
    logger_->msgStream(LogDebug2) << me_ << "lam[" << i << "] = "
                                  << lam[i] << std::endl;
  }
  for (UInt i=0; i<problem_->getNumVars(); ++i) {
    logger_->msgStream(LogDebug2) << std::setprecision(8)
                                  << me_ << problem_->getVariable(i)->getName()
                                  << " = " << x[i] << std::endl;
  }
  for (int i=0; i<*l_hess; ++i) {
    logger_->msgStream(LogDebug2) << std::setprecision(8) 
                                  << me_ << "hess[" << i << "] = "
                                  << values[i] << std::endl;
  }
  logger_->msgStream(LogDebug2) << me_ << "error in evalHessian = " 
    << *error << std::endl;
#endif
}


double FilterSQPEngine::evalObjValue(const double *x, int *err)
{
  double objval = problem_->getObjValue(x, err);
#if SPEW
  logger_->msgStream(LogDebug2) << me_
    << "  objective value = " << objval << " error = " << *err << std::endl;
#endif
  return objval;
}


void FilterSQPEngine::freeStorage_()
{
  delete [] c_;
  delete [] ws_;
  delete [] lws2_;
  delete [] s_;
  delete [] x_;
  delete [] lam_;
  delete [] mlam_;
  delete [] a_;
  delete [] rstat_;
  delete [] istat_;
  delete [] bl_;
  delete [] bu_;
  delete [] cstype_;
  delete [] lws_;
  delete [] la_;
}


std::string FilterSQPEngine::getName() const
{
  return "Filter-SQP";
}


ConstSolutionPtr FilterSQPEngine::getSolution() 
{
  return sol_;
}


double FilterSQPEngine::getSolutionValue() 
{
  return sol_->getObjValue();
}


EngineStatus FilterSQPEngine::getStatus() 
{
  return status_;
}
  

ConstWarmStartPtr FilterSQPEngine::getWarmStart()
{
  return warmSt_;
}


WarmStartPtr FilterSQPEngine::getWarmStartCopy()
{
  FilterWSPtr warm_st;
  if (warmSt_) {
    warm_st = (FilterWSPtr) new FilterSQPWarmStart(warmSt_); // copy
  } else {
    warm_st = FilterWSPtr(); // NULL
  }
  return warm_st;
}


void FilterSQPEngine::load(ProblemPtr problem)
{
  problem_ = problem;
  problem->setEngine(this);
}


void FilterSQPEngine::loadFromWarmStart(const WarmStartPtr warm_st)
{
  if (warm_st) {
    // Two important points:
    // 1. dynamic cast can't seem to be avoided.
    // 2. we need to use boost::dynamic_pointer_cast instead of dynamic_cast.
    ConstFilterWSPtr warm_st2 = 
      boost::dynamic_pointer_cast <const FilterSQPWarmStart> (warm_st);

    // now create a full copy.
    warmSt_ = (FilterWSPtr) new FilterSQPWarmStart(warm_st2);
    if (!useWs_) {
      logger_->msgStream(LogInfo) << "setWarmStart() method is called but"
        " warm-start is not enabled." << std::endl;
    }
  } else {
    warmSt_ = FilterWSPtr(); //NULL
  }
}


void FilterSQPEngine::negateObj()
{
  consChanged_ = true;
}


void FilterSQPEngine::removeCons(std::vector<ConstraintPtr> &)
{
  consChanged_ = true;
}


void FilterSQPEngine::resetIterationLimit()
{
  iterLimit_ = maxIterLimit_;
}


void FilterSQPEngine::setBounds_()
{
  VariablePtr vPtr;
  VariableConstIterator vIter;
  UInt i=0;
  UInt n=problem_->getNumVars();
  double l,u;
  for (vIter=problem_->varsBegin(); vIter!=problem_->varsEnd(); ++vIter) {
    vPtr = *vIter;
    l = vPtr->getLb();
    u = vPtr->getUb();
    if (l>u && l<u+bTol_) {u=l;}
    bl_[i] = l;
    bu_[i] = u;
    ++i;
  }
  ConstraintPtr cPtr;
  ConstraintConstIterator cIter;
  i=0;
  for (cIter=problem_->consBegin(); cIter!=problem_->consEnd(); ++cIter) {
    cPtr = *cIter;
    bl_[n+i] = cPtr->getLb();
    bu_[n+i] = cPtr->getUb();
    ++i;
  }
}


void FilterSQPEngine::setIterationLimit(int limit)
{
  if (limit<1) {
    limit = maxIterLimit_;
  }
  iterLimit_ = limit;
}


void FilterSQPEngine::setStorage_(int mxwk, int maxa)
{
  UInt n = problem_->getNumVars();
  UInt m = problem_->getNumCons();
  UInt cnt;
  FunctionType ftype;

  if (consChanged_) {
    sol_ = (SolutionPtr) new Solution(1E20, 0, problem_);
    if (prepareWs_) {
      warmSt_ = (FilterWSPtr) new FilterSQPWarmStart();
      warmSt_->setPoint(sol_);
    }
    if (c_) {
      freeStorage_();
      c_ = 0;
    }
  }
  if (!c_) {
    c_     = new real[m];     // values of the constraint functions.
    ws_    = new real[mxwk];  // workspace.
    lws2_  = new fint[mxwk];  // workspace.
    s_     = new real[n+m];   // scale.
    x_     = new real[n];     // primal solution.
    lam_   = new real[n+m];   // lagrange mults.
    mlam_  = new real[m];   // lagrange mults.
    a_     = new real[maxa];
    rstat_ = new real[7];     // statistics.
    istat_ = new fint[14];    // statistics.
    bl_    = new real[n+m]; // lower bounds for variables and constraints.
    bu_    = new real[n+m]; // upper bounds for variables and constraints.
    cstype_= new char[m];   // is constraint ('L')linear or ('N')nonlinear
    lws_   = new fint[n+problem_->getHessian()->getNumNz()+3];
    la_    = new int[n + problem_->getJacobian()->getNumNz() + m + 3];
    std::fill(c_, c_+m, 0.);
    std::fill(ws_, ws_+mxwk, 0.);
    std::fill(lws2_, lws2_+mxwk, 0);
    std::fill(lam_, lam_+n+m, 0.);
    std::fill(a_, a_+maxa, 0.);
  }

  if (consChanged_) {
    for (UInt i=0; i<n+m; ++i) {
      s_[i] = 1.;
    }
    cnt = 0;
    for (ConstraintConstIterator cIter=problem_->consBegin(); 
         cIter!=problem_->consEnd(); ++cIter, ++cnt) {
      ftype = (*cIter)->getFunctionType();
      if (ftype==Linear || ftype==Constant) {
        cstype_[cnt] = 'L';     // mark constraints linear
      } else {
        cstype_[cnt] = 'N';     // mark constraints nonlinear
      }
    }
  }
}


void FilterSQPEngine::setStructure_()
{
  UInt n        = 0,           // number of variables in problem.
       jac_nnzs = 0,           // number of enteries in Jacobian.
       cons_pos = 0;           // The position where non-zeros of current 
                               // constraint start getting filled up.
  UInt *iRow, *jCol;

  n = problem_->getNumVars();
  jac_nnzs = problem_->getJacobian()->getNumNz();

  // la[0] = no. of spaces in 'a' that are reserved for jacobian values.
  //         If not warm-starting, this is n + (nonzeros in Jac) + 1. The
  //         first 'n' is for gradient of objective. Next (nonzeros in Jac)
  //         store the Jacobian of constraints. The last space ...
  //        
  cons_pos = n + jac_nnzs + 1; // position where we start filling pointers.
  la_[0] = cons_pos;

  // values in la_ need to filled only once.
  for(UInt j = 1; j <= n; ++j) {
    la_[j] = j; // objective gradient is assumed dense
  }
  la_[cons_pos] = 1; // where the obj gradient entries begin in la_.
  ++cons_pos;
  // now, la_[0,1, ..., n] are set

  // Column 'i' of 'A' matrix denotes the gradient of i-th constraint. These
  // have to be stored in sparse format: values in 'a' and variable-index in 
  // 'la'.
  iRow = new UInt[jac_nnzs];
  jCol = new UInt[jac_nnzs];
  problem_->getJacobian()->fillRowColIndices(iRow, jCol);
  // for (UInt i=0; i<jac_nnzs; ++i) {
  //   std::cout << "irow[" << i << "] = " << iRow[i] << " jcol[" << i << "] = " 
  //     << jCol[i] << std::endl;
  // }

  // check if iRow is non-decreasing. Complain otherwise and stop.
#if DEBUG
  if (jac_nnzs>0) {
    for (UInt i=0; i<jac_nnzs-1; ++i) {
      assert(iRow[i] <= iRow[i+1]);
    }
  }
#endif

  // fill up la_
  // assumption: all rows are non-empty. If there is a row that does not have
  // any variables, this setup will fail.
  UInt fill_pos = 1+n;
  UInt old_cons = -1;
  for (UInt i=0; i<jac_nnzs; ++i) {
    la_[fill_pos] = jCol[i] + 1; // fortran expects off by 1?
    if (old_cons == iRow[i]) {
    } else {
       old_cons = iRow[i];
       la_[cons_pos] = fill_pos;
       ++cons_pos;
    }
    ++fill_pos;
  }
  la_[cons_pos] = fill_pos;
  assert(fill_pos==1+n+jac_nnzs);
  delete [] iRow;
  delete [] jCol;
  // la_ is all done.
  

  // fill up lws_
  // phl = 0;
  UInt hess_nnzs = 0;
  hess_nnzs = problem_->getHessian()->getNumNz();
  lws_[0] = 1+hess_nnzs;
  //std::cout << "\nhess_nnz = " << hess_nnzs << std::endl;
  iRow = new UInt[hess_nnzs];
  jCol = new UInt[hess_nnzs];

  // we get a lower triangular matrix which is ordered by row. This is same as
  // upper triangular matrix which is ordered by column. This is what filter
  // wants! yay! Hence we switch iRow and jCol, but only in the next statement
  // and nowhere else.
  problem_->getHessian()->fillRowColIndices(jCol, iRow) ;

  for (UInt i=0; i<hess_nnzs; ++i) {
    //std::cout << "irow[" << i << "] = " << iRow[i] << " jcol[" << i << "] = " 
    //  << jCol[i] << " " << problem_->getVariable(iRow[i])->getName() << " " << 
    //  problem_->getVariable(jCol[i])->getName() << std::endl;
    ++iRow[i];
    ++jCol[i];
    //std::cout << iRow[i] << ", " << jCol[i] << std::endl;
  }
#if DEBUG
  if (hess_nnzs>0) {
    for (UInt i=0; i<hess_nnzs-1; ++i) {
      assert(jCol[i] <= jCol[i+1]);
    }
  }
#endif

  //std::copy(iRow, iRow + hess_nnzs, lws_+1);
  //std::copy(jCol, jCol + hess_nnzs, lws_+1+hess_nnzs);


  UInt curr_col = 1;
  fill_pos = 1;
  lws_[hess_nnzs+curr_col] = fill_pos;
  for (UInt i=0; i<hess_nnzs; ++i) {
    while (jCol[i] > curr_col) {
      ++curr_col;
      lws_[hess_nnzs+curr_col] = fill_pos;
    }
    lws_[fill_pos] = iRow[i];
    ++fill_pos;
  }
  while (curr_col < n) {
    ++curr_col;
    lws_[hess_nnzs+curr_col] = fill_pos;
  }
  ++curr_col;
  lws_[hess_nnzs+curr_col] = fill_pos;

  //for (UInt i=0; i<3+hess_nnzs+n; ++i) {
    //std::cout << "lws_[" << i << "] = " << lws_[i] << std::endl;
  //}

  delete [] iRow;
  delete [] jCol;

  //std::cout << "n = " << n << std::endl;
  //std::cout << "m = " << m << std::endl;
  //std::cout << "non-zeros in Jacobian = " << jac_nnzs << std::endl;
}


EngineStatus FilterSQPEngine::solve()
{

  if (consChanged_==true) {
    problem_->prepareForSolve();
  }
 
  // it is safe to use filter-defined types like fint and real.
  // see "FilterSQPEngineTypes.h"
  fint n       = problem_->getNumVars();
  fint m       = problem_->getNumCons();
  fint maxa    = problem_->getNumJacNnzs() + n; // Jacobian also stores 
                                                // objective grad.
  //fint kmax    = (n<500)?n:500;               // Set dim of null-space 
                                                // to min(500,n)
  fint kmax    = 500;                           // size of reduced Hessian <= n.
  fint maxf    = 100;                           // Max size of filter.
  fint mlp     = 1000;                          // Max level of degeneracy.
  fint Lhess   = problem_->getNumHessNnzs();    // nnz(Hessian).
  fint lh1     = Lhess + 8 + 2*n + m;           // Hessian storage space.
  fint mxwk0   = 2000000;                       // Initial workspace.
  fint mxiwk0  = 500000;                        // 
  fint mxwk    = 21*n + 8*m + mlp + 8*maxf + lh1 + kmax*(kmax+9)/2 + mxwk0;
  fint mxiwk   = 13*n + 4*m + mlp + lh1 + kmax + 113 + mxiwk0;
  fint iprint  = 0;                               // level of output
  fint nout    = 6;  // 6 for output to stdout, 7 for no output
  fint ifail   = 0;                               // error message
  real rho     = 10.0;                            // initial trust-region radius
  real f       = 1E20;                            // function value.
  real fmin    = -1.E20;                // lower bnd on objective; unboundedness
  real user;           // not used.
  fint *iuser;         // used to store pointer back to this class.
  long cstype_len = m;
  const double *initial_point = 0;

  setStorage_(mxwk, maxa);
  setStructure_();

  if (useWs_ && warmSt_ && warmSt_->hasInfo()) {
    // load warm start.
    initial_point = warmSt_->getPoint()->getDualOfVars();
    if (initial_point) {
      std::copy(initial_point, initial_point+n, lam_);
    }
    initial_point = warmSt_->getPoint()->getDualOfCons();
    if (initial_point) {
      std::copy(initial_point, initial_point+m, lam_+n);
    }
    initial_point = warmSt_->getPoint()->getPrimal();
    ifail = 0; 
    //if (consChanged_==false) {
    //  ifail = 0; // should be -1 for warm starting?
    //} else {
    //  ifail = 0; 
    //}
  } else {
    initial_point = problem_->getInitialPoint();
    std::fill(lam_, lam_+m+n, 0.);
  }

  if (initial_point) {
    std::copy(initial_point, initial_point + n, x_);
  } else {
    std::fill(x_, x_+n, 0.);
  }

  // reload bounds if necessary
  if (true==consChanged_) {
    setBounds_();
  }

  // set the status of the engine to unknown.
  status_ = EngineUnknownStatus;

  // initialize iuser
  uintptr_t u = (uintptr_t) this;
  iuser = convertPtrToInt(u);

#if SPEW
  logger_->msgStream(LogDebug2) << me_ << std::endl
    << "  call = " << stats_->calls+1 << std::endl
    << "  n = " << n << std::endl
    << "  m = " << m << std::endl
    << "  maxa = " << maxa << std::endl
    << "  kmax = " << kmax << std::endl
    << "  Lhess = " << Lhess << std::endl;
#endif
  
  // solve NLP by calling filter. x contains the final solution. f contains
  // the objective value.
  timer_->start();
  filtersqp_(&n, &m, &kmax, &maxa, &maxf, &mlp, &mxwk, &mxiwk,
	       &iprint, &nout, &ifail, &rho, x_, c_, &f, &fmin, bl_,
	       bu_, s_, a_, la_, ws_, lws2_, lam_, cstype_, &user, iuser,
	       &iterLimit_, istat_, rstat_, cstype_len);
  consChanged_ = false;

#if SPEW
  logger_->msgStream(LogDebug) << me_ << "fail = " << ifail << std::endl;
  logger_->msgStream(LogDebug) << me_ << "rho = " << rho << std::endl;
#endif

  // set return status from filter
  switch (ifail) {
   case(0):
    status_ = ProvenLocalOptimal;
    break;
   case(1):
    status_ = ProvenUnbounded;
    break;
   case(2):
    status_ = ProvenInfeasible;
    break;
   case(3):
    status_ = ProvenLocalInfeasible;
    break;
   case(4):
    if (rstat_[4]<feasTol_) {
      status_ = ProvenFailedCQFeas;
    } else {
      status_ = ProvenFailedCQInfeas;
    }
    break;
   case(5):
    if (rstat_[4]<feasTol_) {
      status_ = FailedFeas;
    } else {
      status_ = FailedInfeas;
    }
    break;
   case(6):
    status_ = EngineIterationLimit;
    break;
   default:
    status_ = EngineError;
  }

#if SPEW
  logger_->msgStream(LogDebug) << me_ << "status = " << getStatusString() 
    << std::endl;
  logger_->msgStream(LogDebug) << me_ << "obj = " << f << std::endl;
  logger_->msgStream(LogDebug) << me_ << "time taken = " << timer_->query() 
    << std::endl;
  logger_->msgStream(LogDebug) << me_ << "iterations = " << istat_[1]
    << std::endl;
  //logger_->msgStream(LogNone) << me_ << "number of iterations = " 
  //<< istat_[0] << " " << istat_[1] << " " << istat_[2] << " " << istat_[3] 
  //<< " " << istat_[4] << " " << istat_[5] << " " << istat_[6] << " " 
  //<< istat_[7] << " " << istat_[8] << " " << istat_[9] << std::endl;
#endif

  if (true == strBr_) {
    stats_->strCalls += 1;
    stats_->strTime  += timer_->query();
    stats_->strIters += istat_[1];
  } 
  stats_->calls += 1;
  stats_->time  += timer_->query();
  stats_->iters += istat_[1];
  timer_->stop();

  // store the solution
  sol_->setObjValue(f);
  if (false == strBr_ && saveSol_) {
    sol_->setPrimal(x_);
    sol_->setDualOfCons(lam_+n);
    sol_->setDualOfVars(lam_);
  }

  delete [] iuser;
  return status_;
}


void FilterSQPEngine::writeStats(std::ostream &out) const
{
  if (stats_) {
    out << me_ << "total calls            = " << stats_->calls << std::endl
      << me_ << "strong branching calls = " << stats_->strCalls << std::endl
      << me_ << "total time in solving  = " << stats_->time  << std::endl
      << me_ << "time in str branching  = " << stats_->strTime << std::endl
      << me_ << "total iterations       = " << stats_->iters << std::endl
      << me_ << "strong br iterations   = " << stats_->strIters << std::endl;
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
