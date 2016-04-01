//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file LinearHandler.cpp
 * \brief Handle linear contraints, including  simple bound constraints on
 * variables. Implements methods for relaxing, presolving and separating.
 * Should not be used for checking integrality or branching on variables.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iomanip>
#include <iostream>

#include "MinotaurConfig.h"
#include "Branch.h"
#include "BrCand.h"
#include "Constraint.h"
#include "Environment.h"
#include "Function.h"
#include "LinearFunction.h"
#include "Logger.h"
#include "Node.h"
#include "NonlinearFunction.h"
#include "Objective.h"
#include "Option.h"
#include "PreDelVars.h"
#include "PreSubstVars.h"
#include "Relaxation.h"
#include "Solution.h"
#include "SolutionPool.h"
#include "Timer.h"
#include "LinearHandler.h"
#include "VarBoundMod.h"
#include "Variable.h"

//#define SPEW 1
//TODO:
// remove nintmods

using namespace Minotaur;
const std::string LinearHandler::me_ = "LinearHandler: ";

LinearHandler::LinearHandler()
  : env_(EnvPtr()),
    problem_(ProblemPtr()),
    logger_(LoggerPtr()),
    intTol_(1e-6),
    eTol_(1e-8),
    infty_(1e20),
    pStats_(0),
    pOpts_(0)
{
  linVars_.clear();
}


LinearHandler::LinearHandler(EnvPtr env, ProblemPtr problem)
  : env_(env),
    problem_(problem),
    intTol_(1e-6),
    eTol_(1e-8),
    infty_(1e20),
    pStats_(0),
    pOpts_(0)
{
  logger_ = (LoggerPtr) new Logger((LogLevel)(env->getOptions()->
      findInt("handler_log_level")->getValue()));
  pStats_ = new LinPresolveStats();
  pOpts_  = new LinPresolveOpts();
  pOpts_->doPresolve = env->getOptions()->findBool("lin_presolve")->getValue();
  pOpts_->showStats  = env->getOptions()->findBool("lin_show_stats")
    ->getValue();
  pOpts_->maxIters    = 15;
  pOpts_->purgeVars   = true;
  pOpts_->purgeCons   = true;
  pOpts_->dualFix     = true;
  pOpts_->coeffImp    = true;

  pStats_->iters = 0;
  pStats_->varDel = 0;
  pStats_->conDel = 0;
  pStats_->var2Bin = 0;
  pStats_->var2Int = 0;
  pStats_->vBnd = 0;
  pStats_->cBnd = 0;
  pStats_->cImp = 0;
  pStats_->time = 0.;
  pStats_->timeN = 0.;
  pStats_->nMods = 0;

}


LinearHandler::~LinearHandler()
{
  delete pStats_;
  delete pOpts_;
  problem_.reset();
  env_.reset();
  linVars_.clear();
}


void LinearHandler::copyBndsFromRel_(RelaxationPtr rel, ModVector &p_mods)
{
  VarBoundModPtr mod;
  VariablePtr xp;
  VariablePtr xr;

  for (VariableConstIterator it=problem_->varsBegin(); it!=problem_->varsEnd();
       ++it) {
    xp = *it;
    xr = rel->getRelaxationVar(xp);
    if (!xr) {
      continue;
    }
    if (xr->getLb() > xp->getLb()+eTol_) {
      mod = (VarBoundModPtr) new  VarBoundMod(xp, Lower, xr->getLb());
      mod->applyToProblem(problem_);
      p_mods.push_back(mod);
    }
    if (xr->getUb() < xp->getUb()-eTol_) {
      mod = (VarBoundModPtr) new  VarBoundMod(xp, Upper, xr->getUb());
      mod->applyToProblem(problem_);
      p_mods.push_back(mod);
    }
  }
}


void LinearHandler::findLinVars_()
{
  VariablePtr v;
  bool is_lin;
  FunctionPtr of = problem_->getObjective()->getFunction();

  linVars_.clear();
  for (VariableConstIterator vit=problem_->varsBegin(); 
      vit!=problem_->varsEnd(); ++vit) {
    v      = *vit;
    is_lin = true;

    if (of && false==of->isLinearIn(v)) {
      continue;
    }
    
    for (ConstrSet::iterator cit=v->consBegin(); cit!=v->consEnd(); ++cit) {
      if (Linear != (*cit)->getFunction()->getType() && 
          Constant != (*cit)->getFunction()->getType()) {
        is_lin = false;
        break;
      }
    }
    if (true == is_lin) {
      linVars_.push_back(v);
    }
  }
}


void LinearHandler::findAllBinCons_()
{
  ConstraintConstIterator c_iter;
  ConstraintPtr c;
  LinearFunctionPtr lf;
  UInt nbins, nones, npos;

  for (c_iter=problem_->consBegin(); c_iter!=problem_->consEnd(); ++c_iter) {
    c = *c_iter;
    if (Linear!=c->getFunctionType()) {
      continue;
    }
    lf = c->getFunction()->getLinearFunction();
    nbins = 0;
    nones = 0;
    npos  = 0;
    for (VariableGroupConstIterator it = lf->termsBegin(); 
        it != lf->termsEnd(); ++it) {
      if (Binary == it->first->getType()) {
        ++nbins;
      }
      if (fabs(fabs(it->second)-1.0)<eTol_) {
        ++nones;
      }
      if (it->second>eTol_) {
        ++npos;
      }
      if (2 == nbins && 2 == nones && 2 == lf->getNumTerms() 
          && c->getUb()-c->getLb()<eTol_) {
        ++(pStats_->bImpl);
        problem_->setVarType(lf->termsBegin()->first, ImplBin);
        //c->write(std::cout);
      }
    }
  }
}


void LinearHandler::fixToCont_()
{
  VariablePtr v;
  VariableConstIterator v_iter;
  for (v_iter=problem_->varsBegin(); v_iter!=problem_->varsEnd(); ++v_iter) {
    v = (*v_iter);
    if ((v->getType()==Integer || v->getType() == Binary) && 
        fabs(v->getUb() - v->getLb()) < eTol_) {
      problem_->setVarType(v, Continuous);
    }
  }
}


void LinearHandler::separate(ConstSolutionPtr, NodePtr , RelaxationPtr , 
                             CutManager *, SolutionPoolPtr , bool *,
                             SeparationStatus *)
{
  // do nothing.
}


SolveStatus LinearHandler::presolve(PreModQ *pre_mods, bool *changed0)
{
  SolveStatus status = Started;
  bool changed = true;
  ModQ *dmods = 0; // NULL
  Timer *timer = env_->getNewTimer();
  UInt itemp=0;

  timer->start();
  if (false == pOpts_->doPresolve) {
    delete timer;
    return Finished;
  }

  while(changed==true && pStats_->iters < pOpts_->maxIters) {
    changed = false;
#if SPEW
    logger_->msgStream(LogDebug) << me_ << "presolve iteration " 
                                 << pStats_->iters << std::endl;
#endif
    if (true == pOpts_->purgeVars) delFixedVars_(&changed);
    status = checkBounds_(problem_);
    if (status == SolvedInfeasible) {
      delete timer;
      return SolvedInfeasible;
    }
    tightenInts_(problem_, true, &changed, dmods);
    status = varBndsFromCons_(problem_, true, &changed, dmods, &itemp);
    if (status == SolvedInfeasible) {
      delete timer;
      return SolvedInfeasible;
    }
    if (true==pOpts_->dualFix) dualFix_(&changed);
    if (true == pOpts_->purgeVars) {
      delFixedVars_(&changed);
      purgeVars_(pre_mods);
    }
    if (true == pOpts_->purgeCons) problem_->delMarkedCons();
    if (true == pOpts_->purgeCons && true == pOpts_->purgeVars) {
      substVars_(&changed, pre_mods);
    }
    if (true == pOpts_->purgeCons) problem_->delMarkedCons();
    if (true == pOpts_->purgeVars) purgeVars_(pre_mods);
    if (true == pOpts_->purgeVars && pStats_->iters+1 < pOpts_->maxIters) { 
      chkSing_(&changed);
      purgeVars_(pre_mods);
    }
    if (true == pOpts_->purgeCons) {
      dupRows_(&changed);
      problem_->delMarkedCons();
    }
    if (true == pOpts_->coeffImp) coeffImp_(&changed);
    ++(pStats_->iters);
    if (changed) {
      *changed0 = true;
    }
  }
  findAllBinCons_();
  fixToCont_();

  pStats_->time += timer->query();
  delete timer;
  return Finished;
}


SolveStatus LinearHandler::checkBounds_(ProblemPtr p)
{
  VariablePtr v;
  ConstraintPtr c;
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "checking bounds." << std::endl; 
#endif

  for (VariableConstIterator it=p->varsBegin(); it!=p->varsEnd(); 
      ++it) {
    v = *it;
    if (v->getLb() > v->getUb()+eTol_) {
#if SPEW
      logger_->msgStream(LogDebug) << me_ << "infeasible bounds: "; 
      v->write(logger_->msgStream(LogDebug));
#endif
      return SolvedInfeasible;
    }
  }

  for (ConstraintConstIterator it=p->consBegin(); it!=p->consEnd(); ++it) {
    c = *it;
    if (c->getLb() > c->getUb()+eTol_) {
#if SPEW
      logger_->msgStream(LogDebug) << me_ << "infeasible bounds: "; 
      c->write(logger_->msgStream(LogDebug));
#endif
      return SolvedInfeasible;
    }
  }
  return Finished;
}


void LinearHandler::chkSing_(bool *changed)
{
  ConstraintPtr c;
  LinearFunctionPtr lf;
  FunctionPtr of = problem_->getObjective()->getFunction();
  VariablePtr v;
  double coeff;
  bool del_var;

  findLinVars_();
  for (VarQueueConstIter vit=linVars_.begin(); vit!=linVars_.end(); ++vit) {
    v = *vit;
    if (1==v->getNumCons() && !(of->hasVar(v)) && v->getType()==Continuous) {
      c = *(v->consBegin());
      lf = c->getFunction()->getLinearFunction();
      coeff = lf->getWeight(v);
      del_var = false;
      if (coeff>0) {
        if (c->getLb()>-INFINITY && c->getUb()>=INFINITY) {
          problem_->changeBound(c, Lower, c->getLb()-coeff*v->getUb());
          del_var = true;
        } else if (c->getUb()<INFINITY && c->getLb()<=-INFINITY) {
          problem_->changeBound(c, Upper, c->getUb()-coeff*v->getLb());
          del_var = true;
        }
      } else {
        if (c->getLb()>-INFINITY && c->getUb()>=INFINITY) {
          problem_->changeBound(c, Lower, c->getLb()-coeff*v->getLb());
          del_var = true;
        } else  if (c->getUb()<INFINITY && c->getLb()<=-INFINITY) {
          problem_->changeBound(c, Upper, c->getUb()-coeff*v->getUb());
          del_var = true;
        }
      }
      if (del_var) {
        problem_->changeBound(v, 0.0, 0.0);
        problem_->markDelete(v);
        ++(pStats_->varDel);
        *changed = true;
#if SPEW
        logger_->msgStream(LogDebug) << me_ << "variable " << v->getName() 
                                     << " is a singleton. Fixed." << std::endl;
#endif
      }
    } 
  }
}


void LinearHandler::tightenInts_(ProblemPtr p, bool apply_to_prob, 
                                 bool *changed, ModQ *mods)
{
  VariablePtr v;
  double lb, ub;
  VarBoundModPtr mod;
  LinearFunctionPtr lf;
  ConstraintPtr c;

#if SPEW
  logger_->msgStream(LogDebug) << me_ << "tightening bounds." << std::endl; 
#endif

  // variable bounds
  for (VariableConstIterator it=p->varsBegin(); it!=p->varsEnd(); 
      ++it) {
    v = *it;
    if (apply_to_prob) {
      chkIntToBin_(v);
    }
    if (v->getType()==Integer || v->getType() == Binary) {
      lb = v->getLb();
      ub = v->getUb();
      if (lb > -infty_ && fabs(lb - floor(lb+0.5))>intTol_) {
        mod = (VarBoundModPtr) new VarBoundMod(v, Lower, ceil(lb));
        mod->applyToProblem(p);
        if (apply_to_prob) {
          chkIntToBin_(v);
        } else {
          mods->push_back(mod);
#if SPEW
          mod->write(logger_->msgStream(LogDebug));
#endif 
        }
        *changed = true;
      }
      if (ub < infty_ && fabs(ub - floor(ub+0.5))>intTol_) {
        mod = (VarBoundModPtr) new VarBoundMod(v, Upper, floor(ub));
        mod->applyToProblem(p);
        if (apply_to_prob) {
          chkIntToBin_(v);
        } else {
          mods->push_back(mod);
#if SPEW
          mod->write(logger_->msgStream(LogDebug));
#endif 
        }
        *changed = true;
      }
    }
  }

  /// TODO:
  // constraint bounds
  // for (ConstraintConstIterator it=problem_->consBegin();
  //     it!=problem_->consEnd(); ++it) {
  //   if ((*it)->getFunctionType() == Linear) {
  //     c  = (*it);
  //     g = c->getGran();
  //     if (g > eTol_) {
  //       lb = c->getLb();
  //       ub = c->getUb();
  //       newb = ceil(lb/g)*g;
  //       if ( newb > lb + eTol_) {
  //         mod = (VarBoundModPtr) new ConBoundMod(c, Lower, newb);
  //       }
  //       newb = floor(ub/g)*g;
  //       if ( newb < ub - eTol_) {
  //         mod = (VarBoundModPtr) new ConBoundMod(c, Upper, newb);
  //       }
  //     }
  //   }
  // }
}


SolveStatus LinearHandler::varBndsFromCons_(ProblemPtr p, bool apply_to_prob, 
                                            bool *changed, ModQ *mods,
                                            UInt *nintmods)
{
  ConstraintPtr c_ptr;
  bool t_changed;
  SolveStatus status = Started;

#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "bounds from constraints." 
                               << std::endl; 
#endif

  for (ConstraintConstIterator c_iter=p->consBegin(); 
       c_iter!=p->consEnd(); ++c_iter) {
    c_ptr = *c_iter;
    if (c_ptr->getFunctionType() == Linear && DeletedCons!=c_ptr->getState()) {
      t_changed = true;
      while (true == t_changed) {
        status = linBndTighten_(p, apply_to_prob, c_ptr, &t_changed, mods,
                                nintmods);
        if (SolvedInfeasible == status) {
          return SolvedInfeasible;
        }
        if (true == t_changed) {
          *changed = true;
          if (false==apply_to_prob) {
            t_changed = false;
          }
        }
      }
    } else if (c_ptr->getFunctionType() == Constant
               && DeletedCons!=c_ptr->getState()
               && apply_to_prob 
               && pOpts_->purgeCons) {
#if SPEW
      logger_->msgStream(LogDebug) << "constraint " << c_ptr->getName() 
                                   << " is redundant\n";
#endif
      //assert(!"check activity!");
      p->markDelete(c_ptr);
      ++(pStats_->conDel);
    }
  }
  return status;
}


SolveStatus LinearHandler::varBndsFromObj_(ProblemPtr p, double ub, bool apply_to_prob, 
                                           bool *changed, ModQ *mods)
{
  ObjectivePtr o_ptr;
  bool t_changed = true;
  SolveStatus status = Started;
  LinearFunctionPtr lf = LinearFunctionPtr(); // NULL
  double ll, uu, sing_ll, sing_uu;
  UInt nintmods = 0; // unused

#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "bounds from objective" 
                                << std::endl; 
#endif

  o_ptr = p->getObjective();
  if (o_ptr) {
    lf = o_ptr->getLinearFunction();
  }
  if (lf && o_ptr->getFunctionType() == Linear) {
    while (true == t_changed) {
      t_changed = false;
      sing_ll = INFINITY, sing_uu = INFINITY;
      getLfBnds_(lf, &ll, &uu);
      if (ll < -infty_ || uu > infty_) {
        getSingLfBnds_(lf, &sing_ll, &sing_uu);
      }

      if (ll > ub+eTol_) {
        logger_->msgStream(LogDebug) << "Problem infeasible." << std::endl 
                                     << "objective " << o_ptr->getName()
                                     << std::fixed    << std::setprecision(8) 
                                     << " lower = "   << ll 
                                     << " upper = "   << ub 
                                     << std::endl;
        return SolvedInfeasible;
      }
      if (ll > -infty_) {
        updateLfBoundsFromUb_(p, apply_to_prob, lf, ub, ll, false, &t_changed, 
                              mods, &nintmods);
      } else if (sing_ll > -infty_) {
        updateLfBoundsFromUb_(p, apply_to_prob, lf, ub, sing_ll, true, &t_changed, 
                              mods, &nintmods);
      }
      if (true == t_changed) {
        *changed = true;
      }
    }
  }
  return status;
}


void LinearHandler::coeffImp_(bool *changed)
{
  LinearFunctionPtr lf;
  VariablePtr v;
  double ll, uu; // lower and upper bounds implied by bounds on vars.
  double lb, ub; // lower and upper bounds on a constraint.
  double a0;
  ConstraintPtr c;
  bool implic = true;
  double coeftol = 1e-4; // improve coeffs only if improvement is more than
                         // this number.
  double bslack = 1e-4;  // relax bounds by this amount. Otherwise constraints
                         // may get numerically difficult to solve, e.g.
                         // Syn15H, Syn15M02H, Syn40M04H
                         // etc in CMU lib.

  for (ConstraintConstIterator it= problem_->consBegin(); 
      it!=problem_->consEnd(); ++it) {
    c = *it;
    if (!problem_->isMarkedDel(c) && (Linear==(*it)->getFunctionType() 
          && (c->getLb()<=-INFINITY || c->getUb()>=INFINITY))) {
      lf = c->getFunction()->getLinearFunction();
      if (lf->getNumTerms()<2) {
        continue;
      }
      if (lf->getNumTerms()<50) {
        implic = true;
      } else {
        implic = false;
      }
      lb = c->getLb();
      ub = c->getUb();
      getLfBnds_(lf, &ll, &uu);
      for (VariableGroupConstIterator it2=lf->termsBegin();
           it2 != lf->termsEnd(); ++it2) {
        v = it2->first;
        a0 = it2->second;
        if ((Binary == v->getType() || ImplBin == v->getType()) 
            && !problem_->isMarkedDel(v) 
            && (v->getUb()>v->getLb()+0.5)) {
          if (true==implic) {
            ll = uu = 0;
            computeImpBounds_(c, v, 1.0, &ll, &uu);
            ll -= bslack;
            uu += bslack;
            if (a0>0) {
              ll -= a0;
            } else {
              uu -= a0;
            }
          }
          if (uu+a0 < ub-coeftol && uu >= ub) {
            // constraint is redundant when x0=1
            assert(a0 < 0);
            lf->incTerm(v, ub-uu-a0);
            *changed = true;
            ++(pStats_->cImp);
            break;
          } else if (ll+a0 > lb+coeftol && ll <= lb) {
            // constraint is redundant when x0=1
            assert(a0 > 0);
            lf->incTerm(v, lb-ll-a0);
            *changed = true;
            ++(pStats_->cImp);
            break;
          } 
          if (true==implic) {
            ll = uu = 0;
            computeImpBounds_(c, v, 0.0, &ll, &uu);
            if (a0>0) {
              uu += a0;
            } else {
              ll += a0;
            }
          }
          if (uu-a0 < ub-coeftol && uu >= ub) {
            // constraint is redundant when x0=0
            assert(a0 > 0);
            lf->incTerm(v, uu-a0-ub);
            problem_->changeBound(c, Upper, uu-a0);
            *changed = true;
            ++(pStats_->cImp);
            break;
          } else if (ll-a0 > lb+coeftol && ll <= lb) {
            // constraint is redundant when x0=0
            assert(a0 < 0);
            lf->incTerm(v, ll-a0-lb);
            problem_->changeBound(c, Lower, ll-a0);
            *changed = true;
            ++(pStats_->cImp);
            break;
          }
        }
      }
    }
  }
}


void  LinearHandler::computeImpBounds_(ConstraintPtr c, VariablePtr z,
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

  if (zval<0.5) {
    m = (VarBoundModPtr) new VarBoundMod(z, Upper, 0.0);
  } else {
    m = (VarBoundModPtr) new VarBoundMod(z, Lower, 1.0);
  }
  m->applyToProblem(problem_);
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
      m->applyToProblem(problem_);
      mods.push(m);
    }
    if (u1<v->getUb()) {
      m = (VarBoundModPtr) new VarBoundMod(v, Upper, u1);
      m->applyToProblem(problem_);
      mods.push(m);
    }
  }
  c->getFunction()->getLinearFunction()->computeBounds(&ll, &uu);
  *lb = ll;
  *ub = uu;
  while (!mods.empty()) {
    m2 = mods.top();
    mods.pop();
    m2->undoToProblem(problem_);
  }
}


void LinearHandler::dualFix_(bool *changed)
{
  VariablePtr v;
  ConstraintPtr c;
  double w;
  int dir, cdir;
  ModificationPtr mod = ModificationPtr(); //null
  LinearFunctionPtr olf = problem_->getObjective()->getLinearFunction();

  findLinVars_();
  for (VarQueueConstIter vit=linVars_.begin(); vit!=linVars_.end(); ++vit) {
    v = *vit;
    dir = 0;
    cdir = 0;

    if (v->getUb() - v->getLb() < eTol_) {
      continue;
    }

    if (olf) {
      w = olf->getWeight(v);
      if (w > eTol_) {
        dir = -1;
      } else if (w < -eTol_) {
        dir = 1;
      }
    }
    for (ConstrSet::iterator cit=v->consBegin(); cit!=v->consEnd(); ++cit) {
      c = *cit;
      if (DeletedCons==c->getState()) {
        continue;
      }
      w = c->getLinearFunction()->getWeight(v);
      if (c->getLb() > -INFINITY && c->getUb() < INFINITY) {
        dir = 2;
        break;
      }
      if ((w > 0 && c->getLb() > -INFINITY) || 
          (w < 0 && c->getUb() < INFINITY)) {
        cdir = 1; // want to fix to ub
      } else {
        cdir = -1; // want to fix to lb
      }
      if (0 == dir) {
        dir = cdir;
      } else if (cdir != dir) {
        dir = 2;
        break;
      }
    }
    if (1 == dir && v->getUb() < INFINITY) {
      mod = (VarBoundModPtr) new VarBoundMod(v, Lower, v->getUb());
      mod->applyToProblem(problem_);
      *changed = true;
      ++(pStats_->vBnd);
    } else if (-1 == dir && v->getLb() > -INFINITY) {
      mod = (VarBoundModPtr) new VarBoundMod(v, Upper, v->getLb());
      mod->applyToProblem(problem_);
      *changed = true;
      ++(pStats_->vBnd);
    } else if (0 == dir) { // variable never occurs anywhere.
      if (v->getLb()>-INFINITY) {
        mod = (VarBoundModPtr) new VarBoundMod(v, Upper, v->getLb());
        mod->applyToProblem(problem_);
        *changed = true;
        ++(pStats_->vBnd);
      } else if (v->getUb()<INFINITY) {
        mod = (VarBoundModPtr) new VarBoundMod(v, Lower, v->getUb());
        mod->applyToProblem(problem_);
        *changed = true;
        ++(pStats_->vBnd);
      } else {
        mod = (VarBoundMod2Ptr) new VarBoundMod2(v, 0.0, 0.0);
        mod->applyToProblem(problem_);
        *changed = true;
        ++(pStats_->vBnd);
      }
    }
#if SPEW
    if (mod) {
      logger_->msgStream(LogDebug) << me_ << "variable " << v->getName() 
                                   << " fixed by dual fixing" << std::endl;
      mod->write(logger_->msgStream(LogDebug));
      mod.reset();
    }
#endif
  }
}


void LinearHandler::dupRows_(bool *changed)
{
  const UInt n = problem_->getNumVars();
  const UInt m = problem_->getNumCons();
  UInt i,j;
  int err = 0;
  DoubleVector r1, r2;
  DoubleVector h1;
  DoubleVector h2;
  ConstraintPtr c1, c2;
  bool is_deleted;

  r1.reserve(n);
  r2.reserve(n);
  h1.reserve(m);
  h2.reserve(m);

  for (i=0; i<n; ++i) {
    r1.push_back((double) rand()/(RAND_MAX)*10.0);
    r2.push_back((double) rand()/(RAND_MAX)*10.0);
  }

  i=0;
  for (ConstraintConstIterator it=problem_->consBegin();
       it!=problem_->consEnd(); ++it, ++i) {
    c1 = *it;
    if (c1->getFunctionType()==Linear) {
      h1[i] = c1->getActivity(&(r1[0]), &err);
      h2[i] = c1->getActivity(&(r2[0]), &err);
    } else {
      h1[i] = h2[i] = 1e30;
    }
  }

  for (i=0; i<m; ++i) {
    if (h1[i]<1e29) {
      for (j=i+1; j<m; ++j) {
        if (fabs(h1[j]-h1[i])<1e-10 ||
            fabs(h1[j]+h1[i])<1e-10) {
          c1 = problem_->getConstraint(i);
          c2 = problem_->getConstraint(j);
          is_deleted = treatDupRows_(c1, c2, 1.0, changed);
          //c1->write(std::cout);
          //c2->write(std::cout);
          //std::cout << h1[i] << " xxx " << h1[j] << "\n";
          if (is_deleted) {
            h1[j] = h2[j] = 1e30;
          }
        } else if (h1[j] < 1e29 && fabs(h1[i]/h1[j]-h2[i]/h2[j])<1e-10) {
          c1 = problem_->getConstraint(i);
          c2 = problem_->getConstraint(j);
          is_deleted = treatDupRows_(c1, c2, h1[i]/h1[j], changed);
          //c1->write(std::cout);
          //c2->write(std::cout);
          //std::cout << h1[i] << " *** " << h1[j] << "\n";
          if (is_deleted) {
            h1[j] = h2[j] = 1e30;
          }
        }
      }
    }
  }
}


SolveStatus LinearHandler::linBndTighten_(ProblemPtr p, bool apply_to_prob, 
                                          ConstraintPtr c_ptr, bool *changed,
                                          ModQ *mods, UInt *nintmods)
{
  LinearFunctionPtr lf = c_ptr->getLinearFunction();
  double lb = c_ptr->getLb();
  double ub = c_ptr->getUb();
  double ll, uu;
  double sing_ll = -INFINITY, sing_uu = INFINITY;

  *changed = false;

  if (p->isMarkedDel(c_ptr)) {
    return Started;
  }

  assert(lf);
  getLfBnds_(lf, &ll, &uu);
  if (ll < -infty_ || uu > infty_) {
    getSingLfBnds_(lf, &sing_ll, &sing_uu);
  }

  if (apply_to_prob && ll >= lb - eTol_ && uu <= ub + eTol_) {
#if SPEW
    logger_->msgStream(LogDebug) << "constraint " << c_ptr->getName() 
                                 << " is redundant\n";
#endif
    if (pOpts_->purgeCons == true) {
      p->markDelete(c_ptr);
      ++(pStats_->conDel);
    }
    return Started;
  }

  if (ll > lb+eTol_) {
    // change constraint bounds?
  } 
  if (uu < ub-eTol_) {
    // change constraint bounds?
  }

  if (ll > ub+eTol_) {
    logger_->msgStream(LogDebug) << "Problem infeasible." << std::endl 
                                << "constraint " << c_ptr->getName()
                                << std::fixed    << std::setprecision(8) 
                                << " lower = "   << ll 
                                << " upper = "   << ub 
                                << std::endl;
    return SolvedInfeasible;
  }
  if (uu < lb-eTol_) {
    logger_->msgStream(LogDebug) << "Problem infeasible." << std::endl 
                                << "constraint " << c_ptr->getName()
                                << std::fixed    << std::setprecision(8) 
                                << " lower = "   << lb 
                                << " upper = "   << uu 
                                << std::endl;
    return SolvedInfeasible;
  }

  if (lb > -infty_) {
    if (uu < infty_) {
      updateLfBoundsFromLb_(p, apply_to_prob, lf, lb, uu, false, changed, 
                            mods, nintmods);
    } else if (sing_uu < infty_) {
      updateLfBoundsFromLb_(p, apply_to_prob, lf, lb, sing_uu, true, changed, 
                            mods, nintmods);
    }
  }

  if (true == *changed) {
    getLfBnds_(lf, &ll, &uu);
    if (ll < -infty_ || uu > infty_) {
      getSingLfBnds_(lf, &sing_ll, &sing_uu);
    }
  }

  // c_ptr->write(std::cout);
  if (ub < infty_) {
    if (ll > -infty_) {
      updateLfBoundsFromUb_(p, apply_to_prob, lf, ub, ll, false, changed, 
                            mods, nintmods);
    } else if (sing_ll > -infty_) {
      updateLfBoundsFromUb_(p, apply_to_prob, lf, ub, sing_ll, true, changed, 
                            mods, nintmods);
    }
  }
  return Started;
}


void LinearHandler::updateLfBoundsFromLb_(ProblemPtr p, bool apply_to_prob, 
                                          LinearFunctionPtr lf, double lb,
                                          double uu, bool is_sing,
                                          bool *changed, 
                                          ModQ* mods, UInt *nintmods)
{
  ConstVariablePtr cvar;
  VariablePtr var;
  double coef, vlb, vub, nlb, nub;
  VarBoundModPtr mod;

  for (VariableGroupConstIterator it=lf->termsBegin(); it != lf->termsEnd(); 
      ++it) {
    cvar = it->first;
    coef = it->second;
    vlb  = cvar->getLb();
    vub  = cvar->getUb();
    if (coef > eTol_ && (is_sing==false || vub >= infty_)) {
      if (vub >= infty_) {
        vub = 0.; // should be zero when we have a singleton infinity.
      }
      nlb = (lb - uu)/coef + vub;
      if (nlb > vlb + eTol_) {
        // TODO: remove this hack of converting a const var to var.
        var  = p->getVariable(cvar->getIndex());
        if (nlb > var->getUb()-eTol_) {
          nlb = var->getUb();
        }
        mod = (VarBoundModPtr) new VarBoundMod(var, Lower, nlb);
        mod->applyToProblem(p);
#if SPEW
        logger_->msgStream(LogDebug2) << "mod 1: ";
        mod->write(logger_->msgStream(LogDebug2));
        logger_->msgStream(LogDebug2) << std::endl;
#endif
        if (apply_to_prob) {
          //mod->write(std::cout);
          //lf->write(std::cout);
          //std::cout << std::endl;
          chkIntToBin_(var);
          ++(pStats_->vBnd);
        } else {
          mods->push_back(mod);
        }
        if (var->getType()==Binary||var->getType()==Integer) {
          ++(*nintmods);
        }
        *changed = true;
      } 
      nlb = INFINITY; // for safety only.
    } else if (coef < -eTol_ && (is_sing==false || vlb <= -infty_)) {
      if (vlb <= -infty_) {
        vlb = 0.; // should be zero when we have a singleton infinity.
      }
      nub = (lb - uu)/coef + vlb;
      if (nub < vub - eTol_) {
        // TODO: remove this hack of converting a const var to var.
        var  = p->getVariable(cvar->getIndex());
        if (nub < var->getLb()+eTol_) {
          nub = var->getLb();
        }

        mod = (VarBoundModPtr) new VarBoundMod(var, Upper, nub);
        mod->applyToProblem(p);
#if SPEW
        logger_->msgStream(LogDebug2) << "mod 2: ";
        mod->write(logger_->msgStream(LogDebug2));
        logger_->msgStream(LogDebug2) << std::endl;
#endif
        if (apply_to_prob) {
          //mod->write(std::cout);
          //lf->write(std::cout);
          //std::cout << std::endl;
          ++(pStats_->vBnd);
          chkIntToBin_(var);
        } else {
          mods->push_back(mod);
        }
        if (var->getType()==Binary||var->getType()==Integer) {
          ++(*nintmods);
        }
        *changed = true;
      }
      nub = -INFINITY; // for safety only.
    }
  }
}


void LinearHandler::updateLfBoundsFromUb_(ProblemPtr p, bool apply_to_prob, 
                                          LinearFunctionPtr lf, double ub,
                                          double ll, bool is_sing,
                                          bool *changed, ModQ *mods,
                                          UInt *nintmods)
{
  ConstVariablePtr cvar;
  VariablePtr var;
  double coef, vlb, vub, nlb, nub;
  VarBoundModPtr mod;

  for (VariableGroupConstIterator it=lf->termsBegin(); it!=lf->termsEnd(); 
      ++it) {
    cvar = it->first;
    coef = it->second;
    vlb  = cvar->getLb();
    vub  = cvar->getUb();
    if (coef > eTol_ && (is_sing == false || vlb <= -infty_)) {
      if (vlb <= -infty_) {
        vlb = 0.; // should be zero when we have a singleton infinity.
      }
      nub = (ub - ll)/coef + vlb;
      if (nub < vub - eTol_) {
        // TODO: remove this hack of converting a const var to var.
        var  = p->getVariable(cvar->getIndex());
        if (nub < var->getLb()+eTol_) {
          nub = var->getLb();
        }
        mod = (VarBoundModPtr) new VarBoundMod(var, Upper, nub);
        mod->applyToProblem(p);
#if SPEW
        logger_->msgStream(LogDebug2) << "mod 3: ";
        mod->write(logger_->msgStream(LogDebug2));
        logger_->msgStream(LogDebug2) << std::endl;
#endif
        if (apply_to_prob) {
          //mod->write(std::cout);
          //lf->write(std::cout);
          //std::cout << std::endl;
          ++(pStats_->vBnd);
          chkIntToBin_(var);
        } else {
          mods->push_back(mod);
        }
        if (var->getType()==Binary||var->getType()==Integer) {
          ++(*nintmods);
        }
        *changed = true;
      } 
      nub = -INFINITY; // for safety only.
    } else if (coef < -eTol_ && (is_sing == false || vub >= infty_)) {
      if (vub >= infty_) {
        vub = 0.; // should be zero when we have a singleton infinity.
      }
      nlb = (ub - ll)/coef + vub;
      if (nlb > vlb + eTol_) {
        // TODO: remove this hack of converting a const var to var.
        var  = p->getVariable(cvar->getIndex());
        if (nlb > var->getUb()-eTol_) {
          nlb = var->getUb();
        }
        mod = (VarBoundModPtr) new VarBoundMod(var, Lower, nlb);
        mod->applyToProblem(p);
#if SPEW
        logger_->msgStream(LogDebug2) << "mod 4: ";
        mod->write(logger_->msgStream(LogDebug2));
        logger_->msgStream(LogDebug2) << std::endl;
#endif
        if (apply_to_prob) {
          //mod->write(std::cout);
          //lf->write(std::cout);
          //std::cout << std::endl;
          ++(pStats_->vBnd);
          chkIntToBin_(var);
        } else {
          mods->push_back(mod);
        }
        if (var->getType()==Binary||var->getType()==Integer) {
          ++(*nintmods);
        }
        *changed = true;
      }
      nlb = INFINITY; // for safety only.
    }
  }
}


void LinearHandler::getLfBnds_(LinearFunctionPtr lf, double *lo, double *up)
{
  double lb = 0;
  double ub = 0;
  double coef, vlb, vub;

  for (VariableGroupConstIterator it = lf->termsBegin(); it != lf->termsEnd(); 
      ++it) {
    coef = it->second;
    vlb = it->first->getLb();
    vub = it->first->getUb();
    if (coef>0) {
      lb += coef*vlb;
      ub += coef*vub;
    } else {
      lb += coef*vub;
      ub += coef*vlb;
    }
  }
  *lo = lb;
  *up = ub;
}


void LinearHandler::getSingLfBnds_(LinearFunctionPtr lf, double *lo, 
                                   double *up)
{
  double lb = 0;
  double ub = 0;
  double coef, vlb, vub;
  bool lo_is_sing = false, up_is_sing = false;
  bool lo_is_finite = true, up_is_finite = true;

  for (VariableGroupConstIterator it = lf->termsBegin(); it != lf->termsEnd(); 
      ++it) {
    coef = it->second;
    vlb = it->first->getLb();
    vub = it->first->getUb();
    if (coef>eTol_) {
      if (vub < infty_ && true == up_is_finite) {
        ub += coef*vub;
      } else if (true == up_is_sing) {
        up_is_sing = false;
        ub = INFINITY;
        up_is_finite = false;
      } else {
        up_is_sing = true;
      }

      if (vlb > -infty_ && true == lo_is_finite) {
        lb += coef*vlb;
      } else if (true == lo_is_sing) {
        lo_is_sing = false;
        lb = -INFINITY;
        lo_is_finite = false;
      } else {
        lo_is_sing = true;
      }
    } else if (coef<-eTol_) {
      if (vub < infty_ && true==lo_is_finite) {
        lb += coef*vub;
      } else if (true == lo_is_sing) {
        lo_is_sing = false;
        lb = -INFINITY;
        lo_is_finite = false;
      } else {
        lo_is_sing = true;
      }

      if (vlb > -infty_ && true==up_is_finite) {
        ub += coef*vlb;
      } else if (true == up_is_sing) {
        up_is_sing = false;
        ub = INFINITY;
        up_is_finite = false;
      } else {
        up_is_sing = true;
      }
    }
  }
  *lo = lb;
  *up = ub;
}


void LinearHandler::delFixedVars_(bool *changed)
{
  VariablePtr v;

#if SPEW
  logger_->msgStream(LogDebug) << me_ << "finding fixed variables." 
                               << std::endl; 
#endif
  for (VariableConstIterator it=problem_->varsBegin(); it!=problem_->varsEnd();
      ++it) {
    v = *it;
    if (fabs(v->getUb()-v->getLb()) < eTol_) {
      problem_->markDelete(v);
      ++(pStats_->varDel);
      *changed = true;
#if SPEW
      logger_->msgStream(LogDebug) << me_ << "fixed variable " 
                                   << v->getName() << " at value "
                                   << v->getLb() << std::endl; 
#endif
    }
  }
}


void LinearHandler::relax_(ProblemPtr p, RelaxationPtr rel, bool *is_inf)
{
  ConstraintConstIterator c_iter;
  ObjectivePtr oPtr;
  VariableConstIterator v_iter;
  VariablePtr v, v2;
  FunctionPtr f, newf;
  int err = 0;

  *is_inf = false;

  // clone all the variables one by one.
  for (v_iter=p->varsBegin(); v_iter!=p->varsEnd(); ++v_iter) {
    v2 = (*v_iter);
    v = rel->newVariable(v2->getLb(), v2->getUb(), v2->getType(), 
                         v2->getName(), v2->getSrcType());
  }

  // If the objective function is linear, add it.
  oPtr = p->getObjective();
  if (oPtr) {
    f = oPtr->getFunction();
    if (f && Linear==f->getType()) {
      // create a clone of this linear function.
      newf = f->cloneWithVars(rel->varsBegin(), &err);
      rel->newObjective(newf, oPtr->getConstant(), Minimize);
    } else if (!f) {
      newf.reset();
      rel->newObjective(newf, oPtr->getConstant(), Minimize);
    } else { 
      // NULL
    }
  } else {
    newf = (FunctionPtr) new Function(); 
    rel->newObjective(newf, 0.0, Minimize);
  }

  // clone all constraints that have linear functions only. Do not add
  // constraints that have some linear functions in them. e.g. if p has
  // a constraint:
  // x_0^2 + x_0 + x_1 + x_2 \leq 3,
  // We don't add anything here
  for (c_iter=p->consBegin(); c_iter!=p->consEnd(); ++c_iter) {
    f = (*c_iter)->getFunction();
    if (Linear == f->getType()) {
      // create a clone of this linear function.
      newf = f->cloneWithVars(rel->varsBegin(), &err);
      rel->newConstraint(newf, (*c_iter)->getLb(), (*c_iter)->getUb());
    }
  }
  //std::cout << "IntVarHandler gave this rel:\n";
  //rel->write(std::cout);
  //std::cout << "\nEnd of relaxation from IntVarHandler.\n";
}


void LinearHandler::relaxInitFull(RelaxationPtr rel, bool *is_inf)
{
  relax_(problem_, rel, is_inf);
}


void LinearHandler::relaxInitInc(RelaxationPtr rel, bool *is_inf)
{
  relax_(problem_, rel, is_inf);
}


void LinearHandler::relaxNodeFull(NodePtr, RelaxationPtr rel, bool *is_inf) 
{
  relax_(problem_, rel, is_inf);
}


void LinearHandler::relaxNodeInc(NodePtr, RelaxationPtr, bool *)
{
  // Do nothing.
}


void LinearHandler::substVars_(bool *, PreModQ *mods)
{
  ConstraintPtr c;
  LinearFunctionPtr lf;
  VariablePtr v1, v2, out, in;
  VariableType intype, outtype;
  double a1, a2;
  VariableGroupConstIterator git;
  PreSubstVarsPtr smod = (PreSubstVarsPtr) new PreSubstVars();
  VarBoundModPtr mod;

#if SPEW
  logger_->msgStream(LogDebug) << me_ << "substituting variables."
                               << std::endl; 
#endif

  for (ConstraintConstIterator it= problem_->consBegin(); 
      it!=problem_->consEnd(); ++it) {
    c = *it;
    if (c->getFunctionType()==Linear && 
        fabs(c->getUb()) < eTol_ && fabs(c->getLb()) < eTol_ &&
        c->getLinearFunction()->getNumTerms() == 2) {
      lf = c->getLinearFunction();
      git = lf->termsBegin();
      v1 = (git)->first;
      a1 = (git)->second;
      ++git;
      v2 = (git)->first;
      a2 = (git)->second;

      if (fabs(a1 + a2) < eTol_) {
        if (v2->getFunType()!=Linear && v1->getFunType()==Linear) {
          in = v2;
          out = v1;
        } else {
          in = v1;
          out = v2;
        }
        // the 'in' variable should get the right type.
        intype = in->getType();
        outtype = out->getType();
        if (intype==Continuous) {
          problem_->setVarType(in, outtype);
        } else if (intype==ImplInt && 
            (outtype==ImplBin || outtype==Integer || outtype==Binary)) {
          problem_->setVarType(in, outtype);
        } else if (intype==ImplBin && (outtype==Binary)) {
          problem_->setVarType(in, outtype);
        } else if (intype==Integer && (outtype==Binary || outtype==ImplBin)) {
          problem_->setVarType(in, Binary);
        }
        // the 'in' variable should get the right bounds.
        if (in->getLb() < out->getLb()) {
          mod = (VarBoundModPtr) new VarBoundMod(in, Lower, out->getLb());
          mod->applyToProblem(problem_);
        }
        if (in->getUb() > out->getUb()) {
          mod = (VarBoundModPtr) new VarBoundMod(in, Upper, out->getUb());
          mod->applyToProblem(problem_);
        }
#if SPEW
        logger_->msgStream(LogDebug) << me_ << "substituting " 
                                     << out->getName() << " in constraint " 
                                     << c->getName() << " by " << in->getName() 
                                     << std::endl;
#endif
        problem_->subst(out, in);
        problem_->markDelete(c);
        problem_->markDelete(out);
        ++(pStats_->varDel);
        ++(pStats_->conDel);
        smod->insert(out, in);
      } else if (v1->getType() == Continuous && v2->getType() == Continuous) {
        double rat = 1.0;
        if (v1->getNumCons()<v2->getNumCons()) {
          in  = v1;
          out = v2;
          rat = -a1/a2;
        } else {
          in  = v2;
          out = v1;
          rat = -a2/a1;
        }
        if (rat>0) {
          a1 = out->getLb()/rat;
          a2 = out->getUb()/rat;
        } else {
          a1 = out->getUb()/rat;
          a2 = out->getLb()/rat;
        }

        if (in->getLb() < a1) {
          mod = (VarBoundModPtr) new VarBoundMod(in, Lower, a1);
          mod->applyToProblem(problem_);
        }
        if (in->getUb() > a2) {
          mod = (VarBoundModPtr) new VarBoundMod(in, Upper, a2);
          mod->applyToProblem(problem_);
        }
        problem_->subst(out, in, rat);
        problem_->markDelete(c);
        problem_->markDelete(out);
        ++(pStats_->varDel);
        ++(pStats_->conDel);
        smod->insert(out, in, rat);
#if SPEW
        logger_->msgStream(LogDebug) << me_ << "substituting " 
                                     << out->getName() << " in constraint " 
                                     << c->getName() << " by " << in->getName() 
                                     << std::endl;
#endif
      }  
    }
  }
  if (smod->getSize()>0) {
    mods->push_front(smod);
  }
}


void LinearHandler::purgeVars_(PreModQ *pre_mods)
{
  VariablePtr v = VariablePtr(); // NULL

#if SPEW
  logger_->msgStream(LogDebug) << me_ << "removing variables."  
                               << std::endl; 
#endif
  // Aways push in a NULL pointer first. Helps identify the set of variables
  // that were deleted together.
  //preDelVars_.push_front(v);

  if (problem_->getNumDVars()>0) {
    PreDelVarsPtr dmod = (PreDelVarsPtr) new PreDelVars();
    for (VariableConstIterator it=problem_->varsBegin(); 
        it!=problem_->varsEnd(); ++it) {
      v = *it;
      if (problem_->isMarkedDel(v)) {
        dmod->insert(v);
        //preDelVars_.push_front(v);
      }
    }
    //std::cout << "**** " << //preDelVars_.size() << " ****\n";
    problem_->delMarkedVars();
    pre_mods->push_front(dmod);
  }
}


bool LinearHandler::presolveNode(RelaxationPtr rel, NodePtr node,
                                 SolutionPoolPtr spool, ModVector &p_mods,
                                 ModVector &r_mods)
{
  SolveStatus status = Started;
  if (false == pOpts_->doPresolve || 
      (node->getDepth()>3 && node->getId()%5!=0)) {
    return false;
  }

  simplePresolve(rel, spool, r_mods, status);
  if (true==modProb_) {
    copyBndsFromRel_(rel, p_mods);
  }
  return (status==SolvedInfeasible);
}


void LinearHandler::simplePresolve(ProblemPtr p, SolutionPoolPtr spool,
                                   ModVector &t_mods, SolveStatus &status) 
{
  bool changed = true;
  ModQ mods;
  UInt max_iters = 10;
  UInt min_iters = 2;
  UInt iters = 1;
  UInt nintmods;
  Timer *timer = 0;

  timer = env_->getNewTimer();
  timer->start();
  while (true==changed && iters<=max_iters && 
         (iters<=min_iters || nintmods>0) && 
         status!=SolvedInfeasible) {
    nintmods = 0;
    changed = false;
    ++iters;
    varBndsFromCons_(p, false, &changed, &mods, &nintmods);
    if (spool && spool->getNumSols()>0) {
      varBndsFromObj_(p, spool->getBestSolutionValue(), false, &changed, &mods);
    }
    tightenInts_(p, false, &changed, &mods); 
    status = checkBounds_(p);
  }
  
  for (ModQ::const_iterator it=mods.begin(); it!=mods.end(); ++it) {
    t_mods.push_back(*it);
  }
  pStats_->nMods += mods.size();
  pStats_->timeN += timer->query();

  delete timer;
}


void LinearHandler::chkIntToBin_(VariablePtr v)
{
  double lb = v->getLb();
  double ub = v->getUb();

  if (v->getType()==Integer && lb > -eTol_ && ub < 1+eTol_) {
    problem_->setVarType(v, Binary);
  }
}


const LinPresolveOpts* LinearHandler::getOpts() const
{
  return pOpts_;
}


bool LinearHandler::treatDupRows_(ConstraintPtr c1, ConstraintPtr c2,
                                  double mult, bool *changed)
{
  LinearFunctionPtr lf1 = c1->getFunction()->getLinearFunction();
  LinearFunctionPtr lf2 = c2->getFunction()->getLinearFunction();
  VariableGroupConstIterator it1, it2;
  bool b1=true; // if the two lf are alike.
  bool b2=true; // if one lf is negative of other.
  bool b3=true; // if one lf is some multiple of other.
  double d1, d2;
  double lb, ub;

  if (lf1->getNumTerms()!=lf2->getNumTerms()) {
    return false;
  }
  for (it1=lf1->termsBegin(), it2=lf2->termsBegin(); it1!=lf1->termsEnd();
       ++it1, ++it2) {
    if (it1->first!=it2->first) {
      return false;
    }
    d1 = it1->second;
    d2 = it2->second;
    if (fabs(d1-d2)>1e-12) {
      b1=false;
    }
    if (fabs(d1+d2)>1e-12) {
      b2=false;
    }
    if (fabs(d1/d2-mult)>1e-12) {
      b3=false;
    }
  }

  if (true==b1 || true==b2 || true==b3) {
    if (true==b1) {
      mult = 1.0;
    } else if (true==b2) {
      mult = -1.0;
    //} else {
    //  c1->write(std::cout);
    //  c2->write(std::cout);
    }
    if (mult>0) {
      lb = mult*c2->getLb();
      ub = mult*c2->getUb();
    } else {
      lb = mult*c2->getUb();
      ub = mult*c2->getLb();
    }
    lb = (c1->getLb()<lb)?lb:c1->getLb();
    ub = (c1->getUb()<ub)?c1->getUb():ub;
    problem_->changeBound(c1, lb, ub);
    problem_->markDelete(c2);
    ++(pStats_->conDel);
    *changed = true;
    return true;
  }
  return false;
}


void LinearHandler::writePreStats(std::ostream &out) const
{
  out << me_ << "Statistics for presolve by Linear Handler:"        << std::endl
    << me_ << "Number of iterations           = "<< pStats_->iters  << std::endl
    << me_ << "Time taken in initial presolve = "<< pStats_->time   << std::endl
    << me_ << "Time taken in node presolves   = "<< pStats_->timeN  << std::endl
    << me_ << "Number of variables deleted    = "<< pStats_->varDel << std::endl
    << me_ << "Number of constraints deleted  = "<< pStats_->conDel << std::endl
    << me_ << "Number of vars set to binary   = "<< pStats_->var2Bin<< std::endl
    << me_ << "Number of vars set to integer  = "<< pStats_->var2Int<< std::endl
    << me_ << "Times variables tightened      = "<< pStats_->vBnd   << std::endl
    << me_ << "Times constraints tightened    = "<< pStats_->cBnd   << std::endl
    << me_ << "Times coefficients improved    = "<< pStats_->cImp   << std::endl
    << me_ << "Times binary variable relaxed  = "<< pStats_->bImpl  << std::endl
    << me_ << "Changes in nodes               = "<< pStats_->nMods  << std::endl
    ;
}


void LinearHandler::writeStats(std::ostream &out) const
{
  writePreStats(out);
}


std::string LinearHandler::getName() const
{
   return "LinearHandler (Handling linear constraints).";
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
