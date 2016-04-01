//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file NlPresHandler.cpp
 * \brief Presolve nonlinear contraints.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string.h> // for memset

#include "MinotaurConfig.h"
#include "Constraint.h"
#include "Environment.h"
#include "Function.h"
#include "HessianOfLag.h"
#include "LinearFunction.h"
#include "Logger.h"
#include "Node.h"
#include "Problem.h"
#include "VarBoundMod.h"
#include "NlPresHandler.h"
#include "NonlinearFunction.h"
#include "Objective.h"
#include "Option.h"
#include "PreAuxVars.h"
#include "ProblemSize.h"
#include "Relaxation.h"
#include "SolutionPool.h"
#include "Timer.h"
#include "Variable.h"


using namespace Minotaur;
const std::string NlPresHandler::me_ = "NlPresHandler: ";

NlPresHandler::NlPresHandler()
  : doPersp_(false),
    doQuadCone_(false),
    env_(EnvPtr()),
    eTol_(1e-6),
    logger_(LoggerPtr()),
    p_(ProblemPtr()),
    zTol_(1e-6)
{
  stats_.cBnd = 0;
  stats_.cImp = 0;
  stats_.conDel = 0;
  stats_.pRefs = 0;
  stats_.iters = 0;
  stats_.nMods = 0;
  stats_.qCone = 0;
  stats_.time = 0;
  stats_.varDel = 0;
  stats_.vBnd = 0;
}


NlPresHandler::NlPresHandler(EnvPtr env, ProblemPtr p)
  : env_(env),
    eTol_(1e-6),
    p_(p),
    zTol_(1e-6)
{
  logger_ = (LoggerPtr) new Logger((LogLevel)(env->getOptions()->
                                              findInt("handler_log_level")->
                                              getValue()));
  doPersp_ = env->getOptions()->findBool("persp_ref")->getValue();
  doQuadCone_ = env->getOptions()->findBool("quad_cone_ref")->getValue();
  stats_.cBnd = 0;
  stats_.cImp = 0;
  stats_.conDel = 0;
  stats_.pRefs = 0;
  stats_.iters = 0;
  stats_.nMods = 0;
  stats_.qCone = 0;
  stats_.time = 0;
  stats_.varDel = 0;
  stats_.vBnd = 0;
}


NlPresHandler::~NlPresHandler()
{
  env_.reset();
}


void  NlPresHandler::chkRed_(bool *changed)
{
  ConstraintPtr c;
  LinearFunctionPtr lf;
  NonlinearFunctionPtr nlf;
  double lfu, lfl;
  double nlfu, nlfl;
  int error = 0;

  for (ConstraintConstIterator cit=p_->consBegin(); cit!=p_->consEnd();
       ++cit) {
    c = *cit;
    if (c->getFunctionType()!=Linear && c->getFunctionType()!=Constant) {
      lf = c->getFunction()->getLinearFunction();
      nlf = c->getFunction()->getNonlinearFunction();
      lfu = lfl = 0.0;
      nlfu = nlfl = 0.0;
      if (lf) {
        lf->computeBounds(&lfl, &lfu);
      }
      error = 0;
      nlf->computeBounds(&nlfl, &nlfu, &error);
      assert(error==0);
      if (c->getLb()>-INFINITY && nlfl+lfl-eTol_ > c->getLb()) {
        p_->changeBound(c, Lower, -INFINITY);
        *changed = true;
      }
      if (c->getUb()<INFINITY && nlfu+lfu+eTol_ < c->getUb()) {
        p_->changeBound(c, Upper, INFINITY);
        *changed = true;
      }
      if (c->getUb()>=nlfu+lfu-eTol_ && c->getLb()<=nlfl+lfl+eTol_) {
        p_->markDelete(c);
        *changed = true;
        ++stats_.conDel;
      }
    }
  }
}


// TODO: make it work even when binaries are absent from a constraint.
void  NlPresHandler::coeffImpr_(bool *changed)
{
  ConstraintPtr c;
  LinearFunctionPtr lf;
  NonlinearFunctionPtr nlf;
  double ll, uu;
  double cu, cl;
  double a0;
  VariablePtr z;

  if (p_->getSize()->bins<1) {
    return;
  }
  for (ConstraintConstIterator cit=p_->consBegin(); cit!=p_->consEnd();
       ++cit) {
    c = *cit;
    if (c->getFunctionType()==Linear || c->getFunctionType()==Constant) {
      continue;
    }
    if (c->getLb()>-INFINITY && c->getUb()<INFINITY) {
      continue;
    }
    lf   = c->getFunction()->getLinearFunction();
    if (!lf) {
      //c->write(std::cout);
      continue;
    }
    nlf   = c->getFunction()->getNonlinearFunction();
    cu = c->getUb();
    cl = c->getLb();
    for (VariableGroupConstIterator it2=lf->termsBegin();
         it2 != lf->termsEnd(); ++it2) {
      z = it2->first;
      a0 = it2->second;
      if ((z->getType()!=Binary && z->getType()!=ImplBin)
          || fabs(z->getUb()-z->getLb()-1) > eTol_
          || p_->isMarkedDel(z)
          || nlf->hasVar(z)) {
        continue;
      }
      computeImpBounds_(c, z, 0.0, &ll, &uu);
      if (uu < cu-eTol_ && uu+a0 >= cu) {
        // constraint is redundant when z=0
        assert(a0 > 0);
        lf->incTerm(z, uu-cu);
        p_->changeBound(c, Upper, uu);
        *changed = true;
        ++stats_.cImp;
        break;
      } else if (ll > cl+eTol_ && ll+a0 < cl) {
        // constraint is redundant when z=0
        assert(a0 < 0);
        lf->incTerm(z, ll-cl);
        p_->changeBound(c, Lower, ll);
        *changed = true;
        ++stats_.cImp;
        break;
      }
      computeImpBounds_(c, z, 1.0, &ll, &uu);
      if (uu < cu-eTol_ && uu-a0 > cu) {
        // constraint is redundant when z=1
        assert(a0 < 0);
        lf->incTerm(z, cu-uu);
        *changed = true;
        ++stats_.cImp;
        break;
      } else if (ll > cl+eTol_ && ll-a0 < cl) {
        // constraint is redundant when z=1
        assert(a0 > 0);
        lf->incTerm(z, cl-ll);
        *changed = true;
        ++stats_.cImp;
        break;
      }
    }
  }
}


void NlPresHandler::bin2LinF_(ProblemPtr p, LinearFunctionPtr lf,
                              UInt nz, const UInt *irow,
                              const UInt *jcol,
                              const double *values, PreAuxVarsPtr mod)
{
  VariablePtr v1, v2, v3;
  LinearFunctionPtr lf3;
  FunctionPtr f;
  double lb, ub;

  for (UInt i=0; i<nz; ++i) {
    if (fabs(values[i])>1e-12) {
      v1 = p->getVariable(irow[i]);
      v2 = p->getVariable(jcol[i]);
      if (v1==v2) {
        assert(Binary == v1->getType() || ImplBin == v1->getType());
        lf->incTerm(v1, 0.5*values[i]);
      } else if ((Binary == v1->getType() || ImplBin == v1->getType()) &&
                 (Binary == v2->getType() || ImplBin == v2->getType())) {
        v3 = p->newVariable(0.0, 1.0, ImplBin);
        lf->incTerm(v3, values[i]);
        lf3 = (LinearFunctionPtr) new LinearFunction();
        lf3->addTerm(v3, 1.0);
        lf3->addTerm(v1, -1.0);
        f = (FunctionPtr) new Function(lf3);
        p->newConstraint(f, -INFINITY, 0.0);

        lf3 = (LinearFunctionPtr) new LinearFunction();
        lf3->addTerm(v3, 1.0);
        lf3->addTerm(v2, -1.0);
        f = (FunctionPtr) new Function(lf3);
        p->newConstraint(f, -INFINITY, 0.0);

        lf3 = (LinearFunctionPtr) new LinearFunction();
        lf3->addTerm(v1, 1.0);
        lf3->addTerm(v2, 1.0);
        lf3->addTerm(v3, -1.0);
        f = (FunctionPtr) new Function(lf3);
        p->newConstraint(f, -INFINITY, 1.0);
        mod->insert(v3);
      } else {
        if (Binary==v1->getType() || ImplBin==v1->getType()) {
          // exchange
          v3 = v1;
          v1 = v2;
          v2 = v3;
        }
        // v1 is continuous, v2 is binary.
        assert (Binary==v2->getType() || ImplBin==v2->getType());
        lb = (v1->getLb()<0.0) ? v1->getLb() : 0.0;
        ub = (v1->getUb()>0.0) ? v1->getUb() : 0.0;
        v3 = p->newVariable(lb, ub, Continuous);
        lf->incTerm(v3, values[i]);

        lf3 = (LinearFunctionPtr) new LinearFunction();
        lf3->addTerm(v1, -1.0);
        lf3->addTerm(v3, 1.0);
        lf3->addTerm(v2, -v1->getLb());
        f = (FunctionPtr) new Function(lf3);
        p->newConstraint(f, -INFINITY, -v1->getLb());

        lf3 = (LinearFunctionPtr) new LinearFunction();
        lf3->addTerm(v1, 1.0);
        lf3->addTerm(v3, -1.0);
        lf3->addTerm(v2, v1->getUb());
        f = (FunctionPtr) new Function(lf3);
        p->newConstraint(f, -INFINITY, v1->getUb());

        lf3 = (LinearFunctionPtr) new LinearFunction();
        lf3->addTerm(v3, 1.0);
        lf3->addTerm(v2, -v1->getUb());
        f = (FunctionPtr) new Function(lf3);
        p->newConstraint(f, -INFINITY, 0);

        lf3 = (LinearFunctionPtr) new LinearFunction();
        lf3->addTerm(v3, -1.0);
        lf3->addTerm(v2, v1->getUb());
        f = (FunctionPtr) new Function(lf3);
        p->newConstraint(f, -INFINITY, 0);

        mod->insert(v3);
      }
    }
  }
}


void  NlPresHandler::bin2Lin_(ProblemPtr p, PreModQ *mods, bool *changed)
{
  ConstraintPtr c;
  LinearFunctionPtr lf;
  FunctionPtr f;
  HessianOfLagPtr hess;
  double *mult = 0;
  int err = 0;
  double *x = 0;
  double *values = 0;
  double *grad = 0;
  double nlconst = 0;
  UInt *irow = 0;
  UInt *jcol = 0;
  UInt nz = 0;
  PreAuxVarsPtr mod = (PreAuxVarsPtr) new PreAuxVars();

  p->calculateSize();
  if (0==p->getSize()->quadCons && Quadratic!=p->getSize()->objType) {
    return;
  }

  p->setNativeDer(); // TODO: avoid setting it up repeatedly.
  hess = p->getHessian();
  nz = hess->getNumNz();
  mult = new double[p->getNumCons()];
  values = new double[nz];
  irow = new UInt[nz];
  jcol = new UInt[nz];
  x = new double[p->getNumVars()];
  grad = new double[p->getNumVars()];
  memset(mult, 0, p->getNumCons()*sizeof(double));
  memset(values, 0, nz*sizeof(UInt));
  memset(irow, 0, nz*sizeof(UInt));
  memset(jcol, 0, nz*sizeof(UInt));
  memset(x, 0, p->getNumVars()*sizeof(double));

  hess->fillRowColIndices(irow, jcol);

  // objective
  f = p->getObjective()->getFunction();
  if (f) {
    if (Quadratic==f->getType()) {
      hess->fillRowColValues(x, 1.0, mult, values, &err); assert(0==err);
      if (canBin2Lin_(p, nz, irow, jcol, values)) {
        memset(grad, 0, p->getNumVars()*sizeof(double));
        f->evalGradient(x, grad, &err); assert(0==err);
        nlconst = f->eval(x, &err); assert(0==err);
        lf = (LinearFunctionPtr) new LinearFunction(grad, p->varsBegin(),
                                                    p->varsEnd(), 1e-12);
        bin2LinF_(p, lf, nz, irow, jcol, values, mod);
        f = (FunctionPtr) new Function(lf);
        p->newObjective(f, p->getObjective()->getConstant()+nlconst, Minimize);
        *changed = true;
      }
    }
  }
  

  // constraints
  for (ConstraintConstIterator cit=p->consBegin(); cit!=p->consEnd();
       ++cit) {
    c = *cit;
    f = c->getFunction(); 
    if (!f) {
      continue;
    }
    if (Quadratic==f->getType()) {
      mult[c->getIndex()] = 1.0;
      memset(values, 0, nz*sizeof(UInt));
      hess->fillRowColValues(x, 0.0, mult, values, &err); assert(0==err);
      mult[c->getIndex()] = 0.0;
      if (canBin2Lin_(p, nz, irow, jcol, values)) {
        lf = f->getLinearFunction();
        bin2LinF_(p, lf, nz, irow, jcol, values, mod);
        f = (FunctionPtr) new Function(lf);
        p->newConstraint(f, c->getLb(), c->getUb());
        p->markDelete(c);
        *changed = true;
        cit = p->consBegin()+c->getIndex();
      }
    }
  }

  if (mod->getSize()>0) {
    mods->push_back(mod);
  }
  p->delMarkedCons();
  delete [] mult;
  delete [] x;
  delete [] values;
  delete [] grad;
  delete [] irow;
  delete [] jcol;

}


bool NlPresHandler::canBin2Lin_(ProblemPtr p, UInt nz, const UInt *irow,
                                const UInt *jcol, const double *values)
{
  VariablePtr v1, v2;

  for (UInt i=0; i<nz; ++i) {
    if (fabs(values[i])>1e-12) {
      v1 = p->getVariable(jcol[i]);
      v2 = p->getVariable(irow[i]);
      // if neither binary, return false
      if ((v1->getType()==Continuous || v1->getType()==Integer) &&
          (v2->getType()==Continuous || v2->getType()==Integer)) {
        return false;
      }
      // if either one is unbounded, return false
      if ((v1->getUb()>=INFINITY || v1->getLb()<=-INFINITY) ||
          (v2->getUb()>=INFINITY || v2->getLb()<=-INFINITY)) {
        return false;
      }
    }
  }
  return true;
}


#if 0
bool  NlPresHandler::canLin_(FunctionPtr f)
{
  CGraphPtr cg;
  NonlinearFunctionPtr nlf;
  LTHessStor stor;
  UInt i;
  UInt nz;
  UInt *cols;
  std::deque<UInt> *indq;
  VariablePtr x1;
  VariablePtr x2;

  if (f->getQuadraticFunction()) {
    assert(!"can not do explicit quadratic function yet!");
  }
  
  nlf = f->getNonlinearFunction();
  if (!nlf) {
    assert(!"Need nonlinear function!");
  }

  assert(Quadratic==nlf->getType());

  stor.nlVars = nlf->numVars();
  i = 0;
  stor.rows   = new VariablePtr[stor.nlVars];
  stor.colQs  = new std::deque<UInt>[stor.nlVars];
  stor.starts = new UInt[stor.nlVars+1];
  stor.nz = 0;
  for(VariableSet::iterator it=nlf->varsBegin(); it!=nlf->varsEnd(); ++it) {
    stor.rows[i] = (*it);
    ++i;
  }

  nlf->fillHessStor(&stor);

  nz = 0;
  for (i=0; i<stor.nlVars; ++i) {
    stor.starts[i] = nz;
    nz += (stor.colQs+i)->size();
  }
  stor.starts[i] = nz;
  stor.nz = nz;
  stor.cols = new UInt[nz];
  cols = stor.cols;
  indq = stor.colQs;
  nz = 0;
  for (i=0; i<stor.nlVars; ++i, ++indq) {
    x1 = stor.rows[i];
    for (std::deque<UInt>::iterator it2=indq->begin(); it2!=indq->end(); 
         ++it2,++cols) {
      *cols = *it2;
      ++nz;
      x2 = p_->getVar(*cols);
      if (Binary!=x1->getType() && ImplBin!=x1->getType() && 
          Binary!=x2->getType() && ImplBin!=x2->getType()) {
        return false;
      }
    }
  }
  assert(nz == stor.nz);
  nlf->finalHessStor(&stor);
}
#endif


void  NlPresHandler::computeImpBounds_(ConstraintPtr c, VariablePtr z,
                                       double zval, double *lb, double *ub)
{
  VariablePtr v;
  double ll = 0.;
  double uu = 0.;
  double l1, u1, a2, b2;
  ConstraintPtr c2;
  LinearFunctionPtr lf, lf2;
  ModStack mods;
  VarBoundModPtr m;
  ModificationPtr m2;
  int error = 0;

  if (zval<0.5) {
    m = (VarBoundModPtr) new VarBoundMod(z, Upper, 0.0);
  } else {
    m = (VarBoundModPtr) new VarBoundMod(z, Lower, 1.0);
  }
  m->applyToProblem(p_);
  mods.push(m);

  for (VarSet::iterator it=c->getFunction()->varsBegin();
       it!=c->getFunction()->varsEnd(); ++it) {
    v = *it;
    if (v==z) {
      continue;
    }

    l1 = v->getLb();
    u1 = v->getUb();
    for (ConstrSet::iterator it2=v->consBegin(); it2!=v->consEnd(); ++it2) {
      c2 = *it2;
      if (c2->getFunctionType()!=Linear) {
        continue;
      }
      lf2 = c2->getLinearFunction();
      if (lf2->getNumTerms()==2 && lf2->hasVar(z)) {
        b2 = lf2->getWeight(z);
        a2 = lf2->getWeight(v);
        if (a2>0 && (c2->getUb()-zval*b2)/a2 < u1) {
          u1 = (c2->getUb()-zval*b2)/a2;
        }
        if (a2<0 && (c2->getUb()-zval*b2)/a2 > l1) {
          l1 = (c2->getUb()-zval*b2)/a2;
        }
        if (a2>0 && (c2->getLb()-zval*b2)/a2 > l1) {
          l1 = (c2->getLb()-zval*b2)/a2;
        }
        if (a2<0 && (c2->getLb()-zval*b2)/a2 < u1) {
          u1 = (c2->getLb()-zval*b2)/a2;
        }
      }
    }
    if (l1>v->getLb()) {
      m = (VarBoundModPtr) new VarBoundMod(v, Lower, l1);
      m->applyToProblem(p_);
      mods.push(m);
    }
    if (u1<v->getUb()) {
      m = (VarBoundModPtr) new VarBoundMod(v, Upper, u1);
      m->applyToProblem(p_);
      mods.push(m);
    }
  }
  c->getFunction()->getLinearFunction()->computeBounds(&ll, &uu);
  c->getFunction()->getNonlinearFunction()->computeBounds(&l1, &u1, &error);
  ll += l1;
  uu += u1;
  *lb = ll;
  *ub = uu;
  while (!mods.empty()) {
    m2 = mods.top();
    mods.pop();
    m2->undoToProblem(p_);
  }
}


std::string NlPresHandler::getName() const
{
  return "NlPresHandler (presolving nonlinear constraints).";
}


void NlPresHandler::perspRef_(ProblemPtr p, PreModQ *, bool *changed)
{
  ConstraintPtr c, c2;
  FunctionPtr f;
  LinearFunctionPtr lf;
  NonlinearFunctionPtr nlf, pnlf;  
  double alpha, beta;
  double lb, ub;
  double dtmp;
  int err = 0;

  VarSet indvars0;
  VarSet nlvars;
  VarSet indvars;
  VarSet candvars;
  VariablePtr v, z;

  p->calculateSize();

  // first find all binary variables
  for (VariableConstIterator vit=p->varsBegin(); vit!=p->varsEnd(); ++vit) {
    v = *vit;
    if (Binary==v->getType() || ImplBin==v->getType()) {
      indvars0.insert(v);
    }
  }

  // constraints
  for (ConstraintConstIterator cit=p->consBegin(); cit!=p->consEnd();
       ++cit) {
    c = *cit;
    f = c->getFunction(); 
    if (!f) {
      continue;
    }
    if (Quadratic!=f->getType() && Nonlinear!=f->getType()) {
      continue;
    }

    nlvars.clear();
    candvars.clear();
    indvars = indvars0;

    nlf = f->getNonlinearFunction();
    nlvars.insert(nlf->varsBegin(), nlf->varsEnd());

    // search for an indicator variable: We really need an implication graph
    // here.
    for (VarSetIter it=nlvars.begin(); it!=nlvars.end(); ++it) {
      v = *it;
      candvars = indvars;
      indvars.clear();
      for (ConstrSet::iterator cit2=v->consBegin(); cit2!=v->consEnd();
           ++cit2) {
        c2 = *cit2;
        if (c2->getFunctionType()!=Linear) {
          continue;
        }
        lf = c2->getLinearFunction();
        if (lf->getNumTerms()==2) {
          alpha = beta = 0.0;
          alpha = lf->getWeight(v);
          for (VariableGroupConstIterator it2 = lf->termsBegin();
               it2!=lf->termsEnd(); ++it2) {
            z = (it2)->first;
            if (z->getId() != v->getId()) {
              beta = (it2)->second;
              break;
            }
          } 
          if (candvars.end() == candvars.find(z)) {
            continue;
          }
          // check if z turns v off.
          // The constraint is of the form lb <= alpha.v + beta.z <= ub
          lb = c2->getLb();
          ub = c2->getUb();
          if (alpha>0) {
            ub = ub/alpha;
            lb = lb/alpha;
          } else {
            dtmp = lb;
            lb = ub/alpha;
            ub = dtmp/alpha;
          }
          lb = std::max(lb, v->getLb());
          ub = std::min(ub, v->getUb());
          // lb and ub are now bounds on v;
          if (fabs(lb)<zTol_ && fabs(ub)<zTol_) {
            indvars.insert(z);
          }
        }
      }
    }
    if (false == indvars.empty()) {
      // std::cout << "Found it!\n";
      // c->write(std::cout);
      // std::cout << "Indicator Variables:\n";

      // for (VarSet::const_iterator vit=indvars.begin(); vit!=indvars.end();
      //      ++vit) {
      //   (*vit)->write(std::cout);
      // }
      // doing the reformulation now.
      z = *(indvars.begin());
      nlf = c->getFunction()->getNonlinearFunction();
      if (false==nlf->hasVar(z)) {
        pnlf = nlf->getPersp(z, 1e-6, &err); assert(err==0);
        //pnlf->write(std::cout);
        p->changeConstraint(c, pnlf);
        //std::cout << "New constraint:\n";
        //c->write(std::cout);
        *changed = true;
        ++stats_.pRefs;
      }
    } else {
      //std::cout << "Did not find it!\n";
      //c->write(std::cout);
    }
  }
}


SolveStatus NlPresHandler::presolve(PreModQ *mods, bool *changed0)
{
  bool changed = true;
  Timer *tim = env_->getNewTimer();
  SolveStatus status = Started;

  tim->start();
  while(changed==true && stats_.iters < 5) {
    changed = false;
    chkRed_(&changed);
    p_->delMarkedCons();
    status = varBndsFromCons_(&changed);
    if (SolvedInfeasible==status) {
      stats_.time += tim->query();
      delete tim;
      return SolvedInfeasible;
    }
    coeffImpr_(&changed);
    bin2Lin_(p_, mods, &changed);
    if (doPersp_) {
      perspRef_(p_, mods, &changed);
    }
    if (doQuadCone_) {
      quadConeRef_(p_, mods, &changed);
    }
    ++stats_.iters;
    if (changed) {
      *changed0 = true;
    }
  }
  stats_.time += tim->query();
  delete tim;
  return Finished;
}


bool NlPresHandler::presolveNode(RelaxationPtr rel, NodePtr, SolutionPoolPtr s_pool,
                                 ModVector &, ModVector &r_mods)
{
  FunctionPtr f = rel->getObjective()->getFunction();
  NonlinearFunctionPtr nlf;
  LinearFunctionPtr lf;
  double nlfl, nlfu;
  double lfl, lfu;
  double olb;
  int error = 0;
  double a0;
  VariablePtr z;
  double ub = (s_pool)?s_pool->getBestSolutionValue():INFINITY;
  VarBoundModPtr mod;

  if (f && ub<INFINITY) {
    lf = f->getLinearFunction();
    nlf = f->getNonlinearFunction();
    if (lf) {
      lfl = 0.0;
      lfu = 0.0;
      nlfl = 0.;
      nlfu = 0.;
      if (nlf) {
        nlf->computeBounds(&nlfl, &nlfu, &error);
        assert(error==0);
      }
      lf->computeBounds(&lfl, &lfu);
      olb = lfl+nlfl;
      if (olb<=-INFINITY) {
        return false;
      } else if (olb>ub) {
        return true;
      }
      //std::cout << "obj lb = " << lfl+nlfl << " obj ub = " <<
      //lfu+nlfu << " ub = " << ub << std::endl;

      for (VariableGroupConstIterator it2=lf->termsBegin();
           it2 != lf->termsEnd(); ++it2) {
        z = it2->first;
        if ((z->getType()==Binary || z->getType()==ImplBin)
            && (z->getUb()-z->getLb()) > eTol_
            && !(rel->isMarkedDel(z))) {
          a0 = it2->second;
          if (a0>0 && olb+a0>ub) {
            mod = (VarBoundModPtr) new VarBoundMod(z, Upper, 0.0);
            mod->applyToProblem(rel);
            r_mods.push_back(mod);
            ++(stats_.nMods);
          } else if (a0<0 && olb-a0>ub) {
            mod = (VarBoundModPtr) new VarBoundMod(z, Lower, 1.0);
            r_mods.push_back(mod);
            mod->applyToProblem(rel);
            ++(stats_.nMods);
          }
        }
      }
    }
  }
  return false;
}


//
// x_1^2 + x_2^2 + x_3^2 <= K + Mz
// if (M>0), we rewrite it as
// sqrt(x_1^2 + x_2^2 + x_3^2 + eps) <= (1-z).sqrt(K+eps) + z.sqrt(K+M+eps)
// if (M<0), we rewrite it as
// x_1^2 + x_2^2 + x_3^2 <= (K+M) + (1-z).(-M), and then again we get the same 
// sqrt(x_1^2 + x_2^2 + x_3^2 + eps) <= (1-z).sqrt(K+eps) + z.sqrt(K+M+eps)
//
// TODO: implement for qf and K<0
// TODO: implement isSumOfSquares
//
void NlPresHandler::quadConeRef_(ProblemPtr p, PreModQ *, bool *changed)
{
  ConstraintPtr c;
  FunctionPtr f;
  QuadraticFunctionPtr qf;
  LinearFunctionPtr lf, lf2;
  NonlinearFunctionPtr nlf, nlf2;
  VariablePtr z;
  const double eps = 1e-4;
  int err = 0;
  double M, K;

  for (ConstraintConstIterator cit=p->consBegin(); cit!=p->consEnd();
       ++cit) {
    c = *cit;
    K = c->getUb();
    if (K==INFINITY || K<0.0) {
      continue;
    }

    f = c->getFunction(); 
    if (!f) {
      continue;
    }
    if (Quadratic!=f->getType()) {
      continue;
    }
    qf = f->getQuadraticFunction();
    lf = f->getLinearFunction();
    nlf = f->getNonlinearFunction();
    if (!lf || lf->getNumTerms()!=1) {
      continue;
    }

    z = lf->termsBegin()->first;
    if (z->getType()!=Binary && z->getType()!=ImplBin) {
      continue;
    }
    M = -1.0*lf->termsBegin()->second;

    if(nlf && qf) {
      continue;
    } else if (qf) {
      assert(!"implement me!");
    } else if (nlf) {
      if (nlf->isSumOfSquares()) {
        nlf2 = nlf->cloneWithVars(p->varsBegin(), &err);
        if (err) { continue; }
        nlf2->addConst(eps, err); if (err) { continue; }
        nlf2->sqrRoot(err);  if (err) { continue; }

        //c->write(std::cout);
        // nlf2 <= (1-z).sqrt(K+eps) + z.sqrt(K+M+eps)
        lf2 = (LinearFunctionPtr) new LinearFunction();
        lf2->addTerm(z, sqrt(K+eps)-sqrt(K+M+eps));
        f = (FunctionPtr) new Function(lf2, nlf2);

        p->changeConstraint(c, nlf2);
        p->changeConstraint(c, lf2, -INFINITY, sqrt(K+eps));
        ++(stats_.qCone);
        //c->write(std::cout);
        *changed = true;
      }
    }
  }
}


SolveStatus NlPresHandler::varBndsFromCons_(bool *changed)
{
  ConstraintPtr c;
  LinearFunctionPtr lf;
  NonlinearFunctionPtr nlf;
  double lfu, lfl, ub, lb;
  VarBoundModVector mods;
  SolveStatus status = Started;

  for (ConstraintConstIterator cit=p_->consBegin(); cit!=p_->consEnd();
       ++cit) {
    c = *cit;
    if (c->getFunctionType()!=Linear && c->getFunctionType()!=Constant) {
      lf = c->getFunction()->getLinearFunction();
      nlf = c->getFunction()->getNonlinearFunction();
      lfu = lfl = 0.0;
      if (lf) {
        lf->computeBounds(&lfl, &lfu);
      }
      ub = c->getUb()-lfl;
      lb = c->getLb()-lfu;
      nlf->varBoundMods(lb, ub, mods, &status);
      if (SolvedInfeasible == status) {
        mods.clear();
        break;
      } else if (SolveError == status) {
        mods.clear();
        break;
      }
      if (false==mods.empty()) {
        for (VarBoundModVector::iterator it=mods.begin(); it!=mods.end(); ++it) {
          (*it)->applyToProblem(p_);
          ++stats_.vBnd;
#if SPEW
          logger_->msgStream(LogDebug2) << me_ << " ";
          (*it)->write(logger_->msgStream(LogDebug2));
#endif
        }
        mods.clear();
        *changed = true;
      }
    }
  }
  return status;
}


void NlPresHandler::writePreStats(std::ostream &out) const
{
  out << me_ << "Statistics for presolve by NlPresHandler:"        << std::endl
    << me_ << "Number of iterations           = " << stats_.iters  << std::endl
    << me_ << "Time taken in initial presolve = " << stats_.time   << std::endl
    << me_ << "Number of variables deleted    = " << stats_.varDel << std::endl
    << me_ << "Number of constraints deleted  = " << stats_.conDel << std::endl
    << me_ << "Number of perspective reform.  = " << stats_.pRefs  << std::endl
    << me_ << "Times variables tightened      = " << stats_.vBnd   << std::endl
    << me_ << "Times constraints tightened    = " << stats_.cBnd   << std::endl
    << me_ << "Times coefficients improved    = " << stats_.cImp   << std::endl
    << me_ << "Times quad. changed to conic   = " << stats_.qCone  << std::endl
    << me_ << "Changes in nodes               = " << stats_.nMods  << std::endl
    ;
}


void NlPresHandler::writeStats(std::ostream &out) const
{
  writePreStats(out);
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
