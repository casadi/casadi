//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file TransPoly.cpp
 * \brief Define class for reformulating a polynomial problem.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iostream>

#include "MinotaurConfig.h"
#include "CGraph.h"
#include "CNode.h"
#include "Constraint.h"
#include "CxUnivarHandler.h"
#include "Function.h"
#include "LinearFunction.h"
#include "LinearHandler.h"
#include "MultilinearTermsHandler.h"
#include "Objective.h"
#include "PolynomialFunction.h"
#include "Problem.h"
#include "QuadHandler.h"
#include "Solution.h"
#include "TransPoly.h"
#include "YEqMonomial.h"
#include "YEqLFs.h"
#include "YEqUCGs.h"
#include "YEqVars.h"

using namespace Minotaur;

#undef DEBUG_TRANSPOLY

TransPoly::TransPoly(EnvPtr env, ConstProblemPtr p)
  : Transformer(env, p),
    mHandler_(MultilinearTermsHandlerPtr()),
    yMonoms_(0)
{
  env_ = env;
}


TransPoly::~TransPoly() 
{
  if (yMonoms_) {
    delete yMonoms_;
  }
}


void TransPoly::assignMH_(ConstraintPtr con, VariablePtr v, MonomialFunPtr mf)
{
  std::set<ConstVariablePtr> vars;
  for (VarIntMapConstIterator it=mf->termsBegin(); it!=mf->termsEnd(); ++it) {
    vars.insert(it->first);
  }
  mHandler_->addConstraint(con, v, vars);
}


std::string TransPoly::getName() const
{
  return "TransPoly";
}


SolutionPtr TransPoly::getSolOrig(ConstSolutionPtr , int &err)
{
  err = 1;
  return SolutionPtr();
}


SolutionPtr TransPoly::getSolTrans(ConstSolutionPtr , int &err)
{
  err = 1;
  return SolutionPtr();
}


MonomialFunPtr TransPoly::monomToMl_(MonomialFunPtr mf)
{
  MonomialFunPtr mf2 = (MonomialFunPtr) new MonomialFunction(1.0);
  VariablePtr v, v2;
  CGraphPtr cg;
  CNode *n1, *n2;
  LinearFunctionPtr lf;
  ConstraintPtr cnew;
  FunctionPtr f;

  for (VarIntMapConstIterator it=mf->termsBegin(); it!=mf->termsEnd(); ++it) {
    if (1==it->second) {
      v = it->first;
    } else {
      // create a nonlinear function x^k.
      cg = (CGraphPtr) new CGraph();
      n1 = cg->newNode(it->first);
      if (2==it->second) {
        n2 = 0;
        n1 = cg->newNode(OpSqr, n1, n2);
      } else {
        n2 = cg->newNode((double) it->second);
        n1 = cg->newNode(OpPowK, n1, n2);
      }
      cg->setOut(n1);
      cg->finalize();

      // find y = x^k or create it.
      v = yUniExprs_->findY(cg);
      if (!v) {

        // Jeff adding implied bounds on new variable.
        double lb = it->first->getLb();
        double ub = it->first->getUb();
        //XXX Fix this later
        if (lb < 0) {
          std::cerr << "Only solve polynomial problems over positive domain:";
          assert(0);
        }

        lb = pow(lb, double(it->second));
        ub = pow(ub, double(it->second));

        v2 = newp_->newVariable(lb, ub, Continuous); 
 
#if defined(DEBUG_TRANSPOLY)
        std::cout << "Transpoly -- Creating a new variable: ";
        v2->write(std::cout);
        std::cout << "For original variable: ";
        it->first->write(std::cout);
        std::cout << "Final lb: " << lb << " ub: " << ub << std::endl;        
#endif
         
        yUniExprs_->insert(v2, cg);
        lf = (LinearFunctionPtr) new LinearFunction();
        lf->addTerm(v2, -1.0);
        f = (FunctionPtr) new Function(lf, cg);
        cnew = newp_->newConstraint(f, 0.0, 0.0);

        //assignHandler_(cg, cnew);
        // XXX JTL.  Currently for TransPoly, we will give the quadratics to
        // uHandler.
        
        VariablePtr iv;
        iv = *(f->getNonlinearFunction()->varsBegin());
        uHandler_->addConstraint(cnew, iv, v2, 'E');
        v = v2;
      }
    }
    mf2->multiply(1.0, v, 1);
  }
  return mf2;
}


VariablePtr TransPoly::newPolyVar_(MonomialFunPtr mf)
{
  ConstraintPtr newcon;
  CGraphPtr cg;
  LinearFunctionPtr lf;
  FunctionPtr f;
  CNode *n0 = 0;
  MonomialFunPtr mf2 = monomToMl_(mf);
  VariablePtr v;

  if (1<mf2->getDegree()) {
    v = yMonoms_->findY(mf2);
  } else {
    v = mf2->termsBegin()->first;
  }

  assert(mf2);
  if (!v) {
    v = newp_->newVariable();
    yMonoms_->insert(v, mf2);
    cg = (CGraphPtr) new CGraph();
    n0 = mf2->fillCG(cg);
    cg->setOut(n0);
    cg->finalize();
    lf = (LinearFunctionPtr) new LinearFunction();
    lf->addTerm(v, -1.0);
    f = (FunctionPtr) new Function(lf, cg);
    newcon = newp_->newConstraint(f, 0.0, 0.0);
    assignMH_(newcon, v, mf2);
  }
  return v;
}


VariablePtr TransPoly::newPolyVar_(const CNode* cnode, MonomialFunPtr mf,
                                   LinearFunctionPtr lf, VariablePtr v,
                                   double d, double k)
{
  VariablePtr y;
  CGraphPtr cg;
  CNode *n1;
  CNode *n2 = 0;

  if (!v) {
    if (mf) {
      v = newPolyVar_(mf);
    } else if (lf) {
      v = newVar_(lf, d, newp_);
    } else if (v && fabs(d)>zTol_) {
      v = newVar_(v, d, newp_);
    } else if (v && fabs(k-1.0)>zTol_) {
      lf = (LinearFunctionPtr) new LinearFunction();
      lf->addTerm(v, k);
      v = newVar_(lf, 0.0, newp_);
    } else {
      assert(!"constant nonlinear function!!");
    }
    cg = (CGraphPtr) new CGraph();
    n1 = cg->newNode(v);
    n1 = cg->newNode(cnode->getOp(), n1, n2);
    cg->setOut(n1);
    cg->finalize();
  }
  y = yUniExprs_->findY(cg);
  if (!y) {
    cg = (CGraphPtr) new CGraph();
    ConstraintPtr cnew;
    FunctionPtr f;
    n1 = 0;
    n2 = 0;
    LinearFunctionPtr lf2 = (LinearFunctionPtr) new LinearFunction();

    y = newp_->newVariable();
    lf2->addTerm(y, -1.0);

    n1 = cg->newNode(v);
    n1 = cg->newNode(cnode->getOp(), n1, n2);
    cg->setOut(n1);
    cg->finalize();
    f = (FunctionPtr) new Function(lf2, cg);
    cnew = newp_->newConstraint(f, 0.0, 0.0);

    yUniExprs_->insert(y, cg);
    assignHandler_(cg, cnew);
  }
  return y;
}


// Returns one of the following 4:
// #1 mf,
// #2 lf + d, 
// #3 v + d, or
// #4 kv 
// Either k is 0 or d is 0 or both.
void TransPoly::recursPolyRef_(const CNode *node, 
                               MonomialFunPtr &mf, LinearFunctionPtr &lf,
                               VariablePtr &v, double &d, double &k)
{
  double dl, dr, kl, kr;
  LinearFunctionPtr lfl, lfr;
  VariablePtr vl, vr;
  MonomialFunPtr mfl, mfr;
  CNode *n1 = 0;

  lf = LinearFunctionPtr(); // NULL
  v = VariablePtr();
  d = 0.0;
  k = 0.0;

  switch (node->getOp()) {
  case (OpAbs):
  case (OpAcos):
  case (OpAcosh):
  case (OpAsin):
  case (OpAsinh):
  case (OpAtan):
  case (OpAtanh):
  case (OpCeil):
  case (OpCos):
  case (OpCosh):
  case (OpExp):
  case (OpFloor):
  case (OpInt):
  case (OpIntDiv):
  case (OpLog):
  case (OpLog10):
  case (OpRound):
  case (OpSin):
  case (OpSinh):
  case (OpSqrt):
  case (OpTan):
  case (OpTanh):
    recursPolyRef_(node->getL(), mfl, lfl, vl, dl, kl);
    refUnivarOpPoly_(node, mfl, lfl, vl, dl, kl, v, d, k);
    break;
  case (OpCPow):
  case (OpDiv):
    assert(!"not implemented!");
    break;
  case (OpMinus):
    recursPolyRef_(node->getL(), mfl, lfl, vl, dl, kl);
    recursPolyRef_(node->getR(), mfr, lfr, vr, dr, kr);
    refMinus_(mfl, mfr, lfl, lfr, vl, vr, dl, dr, kl, kr,
              mf, lf, v, d, k);
    break;
  case (OpMult):
    recursPolyRef_(node->getL(), mfl, lfl, vl, dl, kl);
    recursPolyRef_(node->getR(), mfr, lfr, vr, dr, kr);
    refMult_(mfl, mfr, lfl, lfr, vl, vr, dl, dr, kl, kr,
             mf, lf, v, d, k);
    break;
  case (OpNone):
    break;
  case (OpNum):
    d = node->getVal();
    break;
  case (OpPlus):
    recursPolyRef_(node->getL(), mfl, lfl, vl, dl, kl);
    recursPolyRef_(node->getR(), mfr, lfr, vr, dr, kr);
    refPlus_(mfl, mfr, lfl, lfr, vl, vr, dl, dr, kl, kr,
             mf, lf, v, d, k);
    break;
  case (OpPow):
    assert(!"OpPow not implemented!");
    break;
  case (OpPowK):
    recursPolyRef_(node->getL(), mfl, lfl, vl, dl, kl);
    recursPolyRef_(node->getR(), mfr, lfr, vr, dr, kr);
    assert(!mfr && !lfr && !vr);
    if (fabs(dr)<zTol_) {
      d = 1.0;
    } else if (fabs(dr-1.0)<zTol_) {
      mf = mfl;
      lf = lfl;
       v =  vl;
       d =  dl;
       k =  kl;
    } else if (dr>0 && fabs(floor(dr+0.5)-dr)<zTol_) { // natural number
      if (mfl) {
        mf = mfl;
        mf->toPower((int) dr);
      } else if (lfl) {
         vl = newVar_(lfl, dl, newp_);
         mf = (MonomialFunPtr) new MonomialFunction(1.0);
         mf->multiply(1.0, vl, (int) dr);
      } else if (vl) {
        if (fabs(dl)>zTol_ && fabs(kl-1.0)>zTol_) {
          lfl = (LinearFunctionPtr) new LinearFunction();
          lfl->addTerm(vl, kl);
          vl = newVar_(lfl, dl, newp_);
        } else if (fabs(dl)>zTol_) {
          vl = newVar_(vl, dl, newp_);
        } 
        mf = (MonomialFunPtr) new MonomialFunction(1.0);
        mf->multiply(1.0, vl, (int) dr);
        mf->multiply(pow(kl, dr));
      } else {
        d = pow(dl, dr);
      }
    } else {
      assert(!"OpPowK not implemented for fractional or negative powers!");
    }
    break;
  case (OpSqr):
    recursPolyRef_(node->getL(), mfl, lfl, vl, dl, kl);
    if (mfl) {
      mf = mfl;
      mf->toPower(2);
    } else if (lfl) {
      vl = newVar_(lfl, dl, newp_);
      mf = (MonomialFunPtr) new MonomialFunction(1.0);
      mf->multiply(1.0, vl, 2);
    } else if (vl) {
      if (fabs(dl)>zTol_ && fabs(kl-1.0)>zTol_) {
        lfl = (LinearFunctionPtr) new LinearFunction();
        lfl->addTerm(vl, kl);
        vl = newVar_(lfl, dl, newp_);
      } else if (fabs(dl)>zTol_) {
        vl = newVar_(vl, dl, newp_);
      } 
      mf = (MonomialFunPtr) new MonomialFunction(1.0);
      mf->multiply(kl*kl, vl, 2);
    } else {
      d = dl*dl;
    }
    break;
  case (OpSumList):
    lf = (LinearFunctionPtr) new LinearFunction();
    for (CNode **it=node->getListL(); it!=node->getListR(); ++it) {
       n1 = *it;
       mfl.reset(); lfl.reset(); vl.reset(); dl = 0; kl = 0;
       recursPolyRef_(n1, mfl, lfl, vl, dl, kl);
       d += dl;
       if (mfl) {
         vl = newPolyVar_(mfl);
         lf->incTerm(vl, mfl->getCoeff());
       } else if (lfl) {
         lf->add(lfl);
       } else if (vl) {
         lf->incTerm(vl, kl);
       }
    }
    break;
  case (OpUMinus):
    recursPolyRef_(node->getL(), mfl, lfl, vl, dl, kl);
    d = -dl;
    if (mfl) {
      mf = mfl;
      mf->multiply(-1.0);
    } else if (lfl) {
      lf = lfl;
      lf->multiply(-1.0);
    } else if (vl) {
      k = -kl;
      v = vl;
    }
    break;
  case (OpVar):
    v = newp_->getVariable(node->getV()->getId());
    k = 1.0;
    break;
  default:
    assert(!"cannot evaluate!");
  }

  assert(!mf || !lf);
  assert(!mf || !v);
  assert(!lf || !v);

  if (lf && lf->getNumTerms()==1 && fabs(d-1.0)<zTol_) { // return return k.v
    k = lf->termsBegin()->second;
    v = lf->termsBegin()->first;
    d = 0.0;
  }
}


void TransPoly::refMinus_(MonomialFunPtr mfl, MonomialFunPtr mfr,
                         LinearFunctionPtr lfl, LinearFunctionPtr lfr,
                         VariablePtr vl, VariablePtr vr,
                         double dl, double dr,
                         double kl, double kr,
                         MonomialFunPtr &mf, LinearFunctionPtr &lf,
                         VariablePtr &v, double &d, double &k) 
{
  d = dl - dr;
  if (!mfr && !lfr && !vr) {
    v = vl;
    lf = lfl;
    mf = mfl;
    k =  kl;
  } else if (!mfl && !lfl && !vl) {
    if (mfr) {
      mf = mfr;
      mf->multiply(-1.0);
    } else if (lfr) {
      lf = lfr;
      lf->multiply(-1.0);
    } else if (vr && fabs(dr)>zTol_) {
      if (fabs(d)<zTol_) {
        v = vr;
        k = -1.0;
        d = 0.0;
      } else {
        assert (vr);
        lf = (LinearFunctionPtr) new LinearFunction();
        lf->addTerm(vr, -1.0);
      }
    } else if (vr) {
      // dl - kr.vr
      if (fabs(dl)<zTol_) {
        v = vr;
        k = -kr;
      } else {
        lf = (LinearFunctionPtr) new LinearFunction();
        lf->addTerm(vr, -kr);
      }
    }
  } else {
    lf = (LinearFunctionPtr) new LinearFunction();
    if (mfl) {
      vl = newPolyVar_(mfl);
      lf->addTerm(vl, mfl->getCoeff());
    } else if (lfl) {
      lf->add(lfl);
    } else if (vl) {
      lf->incTerm(vl, kl);
    }
    if (mfr) {
      vr = newPolyVar_(mfr);
      lf->incTerm(vr, -1.0*mfr->getCoeff());
    } else if (lfr) {
      lfr->multiply(-1.0);
      lf->add(lfr);
    } else if (vr) {
      lf->incTerm(vr, -kr);
    }
    v.reset();
    if (lf->getNumTerms()==1 && fabs(d)<zTol_) {
      v = lf->termsBegin()->first;
      k = lf->termsBegin()->second;
      d = 0.0;
    } else if (lf->getNumTerms()==0) {
      lf.reset();
      k = 0;
    } 
  }
}


void TransPoly::refMult_(MonomialFunPtr mfl, MonomialFunPtr mfr,
                         LinearFunctionPtr lfl, LinearFunctionPtr lfr,
                         VariablePtr vl, VariablePtr vr,
                         double dl, double dr,
                         double kl, double kr,
                         MonomialFunPtr &mf, LinearFunctionPtr &lf,
                         VariablePtr &v, double &d, double &k) 
{
  d = 0.0;
  if (mfl) {
    mf = mfl;
    if (mfr) {
      mf->multiply(mfr);
    } else if (lfr) {
      vr = newVar_(lfr, dr, newp_);
      mf->multiply(1.0, vr, 1);
    } else if (vr) {
      if (fabs(dr)<zTol_) {
        mf->multiply(kr, vr, 1);
      } else {
        vr = newVar_(vr, dr, newp_);
        mf->multiply(1.0, vr, 1);
      }
    } else if (fabs(dr)>zTol_) {
      mf->multiply(dr);
    } else {
      mf.reset();
    }
  } else if (lfl) {
    if (mfr) {
      mf = mfr;
      vl = newVar_(lfl, dl, newp_);
      mf->multiply(1.0, vl, 1);
    } else if (lfr) {
      vl = newVar_(lfl, dl, newp_);
      vr = newVar_(lfr, dr, newp_);
      mf = (MonomialFunPtr) new MonomialFunction(1.0);
      mf->multiply(1.0, vl, 1);
      mf->multiply(1.0, vr, 1);
    } else if (vr && fabs(dr)<zTol_) {
      vl = newVar_(lfl, dl, newp_);
      mf = (MonomialFunPtr) new MonomialFunction(1.0);
      mf->multiply(1.0, vl, 1);
      mf->multiply(kr, vr, 1);
    } else if (vr) {
      vl = newVar_(lfl, dl, newp_);
      vr = newVar_(vr, dr, newp_);
      mf = (MonomialFunPtr) new MonomialFunction(1.0);
      mf->multiply(1.0, vl, 1);
      mf->multiply(1.0, vr, 1);
    } else if (fabs(dr)>zTol_) {
      lf = lfl;
      lf->multiply(dr);
    } 
  } else if (vl) {
    if (fabs(dl)>zTol_) {
      vl = newVar_(vl, dl, newp_);
      kl = 1.0;
    }
    if (mfr) {
      mf = mfr;
      mf->multiply(kl, vl, 1);
    } else if (lfr) {
      vr = newVar_(lfr, dr, newp_);
      mf = (MonomialFunPtr) new MonomialFunction(1.0);
      mf->multiply(kl, vl, 1);
      mf->multiply(1.0, vr, 1);
    } else if (vr && fabs(dr)<zTol_) {
      mf = (MonomialFunPtr) new MonomialFunction(1.0);
      mf->multiply(kl, vl, 1);
      mf->multiply(kr, vr, 1);
    } else if (vr) {
      mf = (MonomialFunPtr) new MonomialFunction(1.0);
      mf->multiply(kl, vl, 1);
      mf->multiply(1.0, vr, 1);
    } else if (fabs(dr)>zTol_) {
      v = vl;
      k = kl*dr;
    }
  } else if (fabs(dl)>zTol_) {
    if (mfr) {
      mf = mfr;
      mf->multiply(dl);
    } else if (lfr) {
      lf = lfr;
      lf->multiply(dl);
    } else if (vr && fabs(dr)>zTol_) {
      v = newVar_(vr, dr, newp_);
      k = dl;
    } else if (vr) {
      v = vr;
      k = kr*dl;
    } else {
      d = dr*dl;
    }
  }
}


void TransPoly::refPlus_(MonomialFunPtr mfl, MonomialFunPtr mfr,
                         LinearFunctionPtr lfl, LinearFunctionPtr lfr,
                         VariablePtr vl, VariablePtr vr,
                         double dl, double dr,
                         double kl, double kr, 
                         MonomialFunPtr &mf, LinearFunctionPtr &lf,
                         VariablePtr &v, double &d, double &k) 
{
  d = dl + dr;
  if (!mfr && !lfr && !vr) {
    mf = mfl;
    lf = lfl;
    v =  vl;
    k =  kl;
  } else if (!mfl && !lfl && !vl) {
    mf = mfr;
    lf = lfr;
    v =  vr;
    k =  kr;
  } else {
    lf = (LinearFunctionPtr) new LinearFunction();
    if (mfl) {
      vl = newPolyVar_(mfl);
      lf->addTerm(vl, mfl->getCoeff());
    } else if (lfl) {
      lf->add(lfl);
    } else {
      lf->incTerm(vl, kl);
    }
    if (mfr) {
      vr = newPolyVar_(mfr);
      lf->incTerm(vr, mfr->getCoeff());
    } else if (lfr) {
      lf->add(lfr);
    } else {
      lf->incTerm(vr, kr);
    }
    v.reset();
    if (lf->getNumTerms()==1 && fabs(d)<zTol_) {
      v = lf->termsBegin()->first;
      k = lf->termsBegin()->second;
      d = 0.0;
      lf.reset();
    } else if (lf->getNumTerms()==0) {
      lf.reset();
      k = 0;
    } 
  }
}


void TransPoly::refNonlinCons_() 
{
  ConstraintPtr c, cnew;
  FunctionPtr f;
  CGraphPtr cg;
  LinearFunctionPtr lf;
  MonomialFunPtr mf;
  VariablePtr v;
  double d;
  double k;

  assert (p_ && newp_);

  for (ConstraintConstIterator it=p_->consBegin(); it!=p_->consEnd();
       ++it) {
    c = *it;
    f = c->getFunction();
    if (f->getType()!=Constant && f->getType()!=Linear) {
      cg = boost::dynamic_pointer_cast <CGraph> (f->getNonlinearFunction());
      mf.reset(); lf.reset(); v.reset(); d = 0; k=0;
      recursPolyRef_(cg->getOut(), mf, lf, v, d, k);
      if (mf) {
        lf = (LinearFunctionPtr) new LinearFunction();
        v = newPolyVar_(mf);
        lf->addTerm(v, mf->getCoeff());
      } else if (v) {
        lf = (LinearFunctionPtr) new LinearFunction();
        lf->addTerm(v, k);
      } 
      if (lf) {
        lf->add(f->getLinearFunction());
        f = (FunctionPtr) new Function(lf);
        cnew = newp_->newConstraint(f, c->getLb()-d, c->getUb()-d);
        lHandler_->addConstraint(cnew);
      } else {
        assert(!"empty constraint?");
      }
    }
  }
}


void TransPoly::refNonlinObj_() 
{
  ObjectivePtr obj;
  FunctionPtr f;
  CGraphPtr cg;
  double d = 0;
  VariablePtr v = VariablePtr();
  LinearFunctionPtr lf;
  MonomialFunPtr mf;
  double k;

  assert(newp_);
  assert(p_);

  obj = p_->getObjective();
  if (!obj) {
    return;
  }
  
  f = obj->getFunction();
  if (f && f->getType()!=Constant && f->getType()!=Linear) {
    cg = boost::dynamic_pointer_cast <CGraph> (f->getNonlinearFunction());
    assert(cg);
    recursPolyRef_(cg->getOut(), mf, lf, v, d, k);
    if (mf) {
      lf = (LinearFunctionPtr) new LinearFunction();
      v = newPolyVar_(mf);
      lf->addTerm(v, mf->getCoeff());
    } else if (v) {
      lf = (LinearFunctionPtr) new LinearFunction();
      lf->addTerm(v, k);
    } 
    if (lf) {
      lf->add(f->getLinearFunction());
      f = (FunctionPtr) new Function(lf);
      newp_->newObjective(f, obj->getConstant()+d, Minimize);
    } else {
      assert(!"empty constraint?");
    }
  }
}


void TransPoly::reformulate(ProblemPtr &newp, HandlerVector &handlers,
                            int &status) 
{
  assert(p_);
  newp = (ProblemPtr) new Problem();
  newp_ = newp;
  yMonoms_ = new YEqMonomial(2*p_->getNumVars());
  yLfs_ = new YEqLFs(2*p_->getNumVars());
  yUniExprs_ = new YEqUCGs();
  lHandler_ = (LinearHandlerPtr) new LinearHandler(env_, newp_);
  handlers.push_back(lHandler_);
  qHandler_ = (QuadHandlerPtr) new QuadHandler(env_, newp_);
  handlers.push_back(qHandler_);
  uHandler_ = (CxUnivarHandlerPtr) new CxUnivarHandler(env_, newp_);
  handlers.push_back(uHandler_);
  mHandler_ = (MultilinearTermsHandlerPtr) new MultilinearTermsHandler(env_,
                                                                       newp_);
  handlers.push_back(mHandler_);

  copyVars_(p_, newp_);
  copyLinear_(p_, newp_);
  refNonlinCons_();
  refNonlinObj_();

  clearUnusedHandlers_(handlers);

  status = 0;
}


void TransPoly::refUnivarOpPoly_(const CNode *node, MonomialFunPtr mfl,
                                 LinearFunctionPtr lfl, VariablePtr vl,
                                 double dl, double kl, 
                                 VariablePtr &v, double &d, double &k)
{
  assert(node!=0);
  if (!mfl && !lfl && !vl) {
    std::cout << "constant function operation in transpoly!!\n";
    d = 0;
    k = 0;
  } else {
    v = newPolyVar_(node, mfl, lfl, vl, dl, kl);
    d = 0;
    k = 1.0;
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
