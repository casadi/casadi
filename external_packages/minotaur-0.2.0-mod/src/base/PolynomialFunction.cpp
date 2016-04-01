//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file PolynomialFunction.h
 * \brief Define the PolynomialFunction class for handling polynomial
 * constraints.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iostream>

#include "MinotaurConfig.h"
#include "CGraph.h"
#include "CNode.h"
#include "LinearFunction.h"
#include "PolynomialFunction.h"
#include "QuadraticFunction.h"
#include "Variable.h"

using namespace Minotaur;

MonomialFunction::MonomialFunction()
: coeff_(0),
  deg_(0),
  eTol_(1e-10)
{
  terms_.clear();
}


MonomialFunction::MonomialFunction(double c)
: coeff_(c),
  deg_(0),
  eTol_(1e-10)
{
  terms_.clear();
}


MonomialFunction::MonomialFunction(double c, ConstVariablePtr v, UInt p)
: coeff_(c),
  deg_(p),
  eTol_(1e-10)
{
  if (p>0) {
    terms_[v] = p;
  }
}


MonomialFunPtr MonomialFunction::clone() const
{
  MonomialFunPtr m2 = (MonomialFunPtr) new MonomialFunction(coeff_);
  m2->terms_ = terms_; // creates copy
  m2->deg_   = deg_;
  return m2;
}


NonlinearFunctionPtr 
MonomialFunction::cloneWithVars(VariableConstIterator vbeg, int *) const
{
  MonomialFunPtr m2 = (MonomialFunPtr) new MonomialFunction(coeff_);
  for (VarIntMapConstIterator it=terms_.begin(); it!=terms_.end(); ++it) {
    m2->terms_[*(vbeg+it->first->getIndex())] = it->second;
  }
  m2->deg_   = deg_;
  return m2;
}


MonomialFunction::~MonomialFunction()
{
  vars_.clear();
  terms_.clear();
}


CNode* MonomialFunction::fillCG(CGraphPtr cg)
{
  CNode *n0 = 0;
  CNode *n1 = 0;
  CNode *n2 = 0;
  for (VarIntMapConstIterator it=terms_.begin(); it!=terms_.end(); 
       ++it) {
    n1 = 0; n2 = 0;
    if (it->second>0) {
      n1 = cg->newNode(it->first);
      if (it->second>1) {
        if (it->second>2) {
          n2 = cg->newNode((int) it->second);
          n1 = cg->newNode(OpPowK, n1, n2);
        } else {
          n1 = cg->newNode(OpSqr, n1, 0);
        }
      } 
    }
    if (0==n0) {
      n0 = n1;
    } else if (0!=n1) {
      n0 = cg->newNode(OpMult, n0, n1);
    }
  }
  if (fabs(coeff_ + 1.0) < eTol_) {
    n0 = cg->newNode(OpUMinus, n0, 0);
  } else if (fabs(coeff_-1.0)>eTol_) {
    n1 = cg->newNode(coeff_);
    n0 = cg->newNode(OpMult, n0, n1);
  }
  return n0;
}


double MonomialFunction::getCoeff() const
{
  return coeff_;
}


int MonomialFunction::getDegree() const
{
  return deg_;
}


VarIntMapConstIterator MonomialFunction::termsBegin()
{
  return terms_.begin();
}


VarIntMapConstIterator MonomialFunction::termsEnd()
{
  return terms_.end();
}


const VarIntMap* MonomialFunction::getTerms() const
{
  return &terms_;
}


void MonomialFunction::multiply(double coeff, ConstVariablePtr v, int p)
{
  if (fabs(coeff) < eTol_) {
    terms_.clear();
    coeff_ = 0;
    deg_   = 0;
  } else if ((int)terms_[v] + p >= 0) {
    terms_[v] += p;
    deg_ += p;
    coeff_ *= coeff;
  } else {
    assert(!"can not have negative powers in monomial function.");
  }
}


void MonomialFunction::multiply(ConstMonomialFunPtr m2)
{
  if (fabs(coeff_) < eTol_) {
    return;
  } else if (!m2 || fabs(m2->coeff_) < eTol_) {
    terms_.clear();
    coeff_ = 0.;
    deg_ = 0;
  } else {
    for (VarIntMapConstIterator it=m2->terms_.begin(); it!=m2->terms_.end(); 
        ++it) {
      terms_[it->first] += it->second;
    }
    coeff_ *= m2->coeff_;
    deg_   += m2->deg_;
  }
}


void MonomialFunction::multiply(double c)
{
  coeff_ *= c;
  if (fabs(coeff_)<eTol_) {
    terms_.clear();
    coeff_ = 0;
    deg_ = 0;
  }
}


void MonomialFunction::operator*=(double c)
{
  multiply(c);
}


void MonomialFunction::operator*=(ConstMonomialFunPtr m2)
{
  multiply(m2);
}


double MonomialFunction::eval(const double *x, int *error)
{
  double prod = 1.;
  *error = 0;
  for (VarIntMapConstIterator it=terms_.begin(); it!=terms_.end(); ++it) {
    prod *= pow(x[it->first->getIndex()],it->second);
  }
  return coeff_*prod;
}


double MonomialFunction::eval(const DoubleVector &x, int *error)
{
  return eval(&x[0], error);
}


void MonomialFunction::evalGradient(const double *x, double *grad_f, 
                                    int *error) 
{
  double prod;
  VariablePtr v;

  *error = 0;
  for (VarIntMapConstIterator it=terms_.begin(); it!=terms_.end(); ++it) {
    prod = coeff_ * pow(x[it->first->getIndex()], it->second-1) * it->second;
    for (VarIntMapConstIterator it2=terms_.begin(); it2!=terms_.end(); ++it2) {
      if (it!=it2) {
        prod *= pow(x[it2->first->getIndex()],it2->second);
      }
    }
    grad_f[it->first->getIndex()] += prod;
  }
}


void MonomialFunction::toPower(int k)
{
  assert(k>=0);
  if (0==k) {
    coeff_ = 1;
    deg_ = 0;
    terms_.clear();
  } else {
    coeff_ = pow(coeff_, k);
    deg_ *= k;
    for (VarIntMapIterator it=terms_.begin(); it!=terms_.end(); ++it) {
      it->second *= k;
    }
  }
}


void MonomialFunction::write(std::ostream &out) const
{
  out << "(" << coeff_ << ")";
  for (VarIntMapConstIterator it=terms_.begin(); it!=terms_.end(); ++it) {
    out << "(" << it->first->getName() << ")^(" << it->second << ")";
  }
}


// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

PolynomialFunction::PolynomialFunction()
 : cb_(0),
   cg_(CGraphPtr()), // NULL
   eTol_(1e-10),
   terms_(0)
{
}


PolynomialFunction::PolynomialFunction(CGraphPtr cg)
 : cb_(0),
   cg_(cg), // NULL
   eTol_(1e-10),
   terms_(0)
{
  if (cg->getOut()) {
    recCG_(cg->getOut(), &cb_, &terms_);
  }
}


PolynomialFunction::~PolynomialFunction()
{
  terms_.clear();
  vars_.clear();
}


void PolynomialFunction::add(ConstMonomialFunPtr m)
{
  if (m) {
    MonomialFunPtr m2 = m->clone();
    // TODO: check for duplicacy.
    terms_.push_back(m2);
  }
}

void PolynomialFunction::clear_()
{
  cb_ = 0;
  terms_.clear();
}


PolyFunPtr PolynomialFunction::clone() const
{
  PolyFunPtr p = (PolyFunPtr) new PolynomialFunction();
  MonomialFunPtr m2;
  for (MonomialConstIter it = terms_.begin(); it!=terms_.end(); ++it) {
    m2 = (*it)->clone();
    p->terms_.push_back(m2);
  }
  p->cb_ = cb_;
  return p;
}


NonlinearFunctionPtr
PolynomialFunction::cloneWithVars(VariableConstIterator vbeg, int *err) const
{
  PolyFunPtr p = (PolyFunPtr) new PolynomialFunction();
  MonomialFunPtr m2;
  for (MonomialConstIter it = terms_.begin(); it!=terms_.end(); ++it) {
    m2 = boost::dynamic_pointer_cast <MonomialFunction> 
      ((*it)->cloneWithVars(vbeg, err));
    p->terms_.push_back(m2);
  }
  p->cb_   = cb_;
  if (cg_) {
    p->cg_ = boost::dynamic_pointer_cast<CGraph> (cg_->cloneWithVars(vbeg, err));
    // hack 
    p->cg_->getVars(&(p->vars_));
  }
  return p;
}


void PolynomialFunction::createCG()
{
  CNode **cnodes;
  CNode *n1;
  UInt size = terms_.size();
  UInt nz=0;

  if (fabs(cb_)<eTol_) {
    ++size;
  }
  cnodes = new CNode*[size];

  cg_ = (CGraphPtr) new CGraph();
  if (fabs(cb_)<eTol_) {
    cnodes[0] = cg_->newNode(cb_);
    ++nz;
  }

  for (MonomialConstIter it = terms_.begin(); it!=terms_.end(); ++it) {
    if (fabs((*it)->getCoeff())>eTol_) {
      cnodes[nz] = (*it)->fillCG(cg_);
      ++nz;
    }
  }
  if (nz>0) {
    n1 = cg_->newNode(OpSumList, cnodes, nz);
  } else {
    n1 = cg_->newNode(0.0);
  }
  cg_->setOut(n1);
  cg_->finalize();

  delete [] cnodes;
}


double PolynomialFunction::eval(const double *x, int *error)
{
  if (cg_) {
    return (cg_->eval(x, error));
  } else {
    double sum = cb_;
    *error = 0;
    for (MonomialConstIter it=terms_.begin(); it!=terms_.end(); ++it) {
      sum += (*it)->eval(x, error);
      if (0!=*error) break;
    }
    return sum;
  }
}


void PolynomialFunction::evalGradient(const double *x, double *grad_f, 
                                      int *error)
{
  *error = 0;
  for (MonomialConstIter it=terms_.begin(); it!=terms_.end(); ++it) {
    (*it)->evalGradient(x, grad_f, error);
    if (0!=*error) break;
  }
}


void PolynomialFunction::evalHessian(const double mult, const double *x, 
                                     const LTHessStor *stor, double *values, 
                                     int *error) 
{
  if (cg_) {
    cg_->evalHessian(mult, x, stor, values, error);
  } else {
    assert(!"Can not get derivatives in polynomialFunction without Cgraph!");
  }
}


void  PolynomialFunction::fillHessStor(LTHessStor *stor)
{
  if (cg_) {
    cg_->fillHessStor(stor);
  } else {
    assert(!"Can not get derivatives in polynomialFunction without Cgraph!");
  }
}


void PolynomialFunction::fillJac(const double *x, double *values, int *error)
{
  if (cg_) {
    cg_->fillJac(x, values, error);
  } else {
    assert(!"Can not get derivatives in polynomialFunction without Cgraph!");
  }
}


void PolynomialFunction::finalHessStor(const LTHessStor *stor)
{
  if (cg_) {
    cg_->finalHessStor(stor);
  } else {
    assert(!"Can not get derivatives in polynomialFunction without Cgraph!");
  }
}


FunctionType PolynomialFunction::getType() const
{
  return Polynomial;
}


void PolynomialFunction::getVars(VariableSet *vars)
{
  if (cg_) {
    cg_->getVars(vars);
  }
}


bool PolynomialFunction::isEmpty() const 
{
  return (fabs(cb_)<eTol_ && terms_.empty());
}


void PolynomialFunction::multiply(ConstLinearFunctionPtr lf, double c)
{
  MonomialFunPtr m;
  MonomialVector terms2 = terms_;
  terms_.clear();
  if (lf) {
    for (VariableGroupConstIterator it=lf->termsBegin(); it!=lf->termsEnd();
        ++it) {
      for (MonomialConstIter it2 = terms2.begin(); it2!=terms2.end(); ++it2) {
        m = (*it2)->clone();
        m->multiply(it->second, it->first, 1);
        terms_.push_back(m);
      }
      if (fabs(cb_)>eTol_) {
        m = (MonomialFunPtr) new MonomialFunction(cb_);
        m->multiply(it->second, it->first, 1);
        terms_.push_back(m);
      }
    }
  }

  if (fabs(c)>eTol_) {
    for (MonomialConstIter it2 = terms2.begin(); it2!=terms2.end(); ++it2) {
      m = (*it2)->clone();
      (*m) *= c;
      terms_.push_back(m);
    }
  }

  cb_ *= c;
  terms2.clear();
}


void PolynomialFunction::multiply(double c)
{
  if (fabs(c)>eTol_) {
    cb_ *= c;
    for (MonomialConstIter it = terms_.begin(); it!=terms_.end(); ++it) {
      **it *= c;
    }
  } else {
    clear_();
  }
}


void PolynomialFunction::prepJac(VarSetConstIter vb, VarSetConstIter ve)
{
  if (cg_) {
    cg_->prepJac(vb, ve);
  } else {
    assert(!"Can not get derivatives in polynomialFunction without Cgraph!");
  }
}


void PolynomialFunction::recCG_(const CNode* cnode, double *c,
                                MonomialVector *terms)
{
  MonomialVector t1, t2;
  double c1 = 0;
  double c2 = 0;
  MonomialFunPtr m;

  switch (cnode->getOp()) {
  case (OpDiv):
    recCG_(cnode->getL(), &c1, &t1);
    recCG_(cnode->getR(), &c2, &t2);
    *c = c1+c2;
    terms->insert(terms->end(), t1.begin(), t1.end());
    break;
  case (OpInt):
    *c = cnode->getVal();
    break;
  case (OpMinus):
    recCG_(cnode->getL(), &c1, &t1);
    recCG_(cnode->getR(), &c2, &t2);
    *c = c1-c2;
    terms->insert(terms->end(), t1.begin(), t1.end());
    for (MonomialConstIter it = t2.begin(); it!=t2.end(); ++it) {
      m = (*it)->clone();
      (*m) *= -1;
      terms->push_back(m);
    }
    break;
  case (OpMult):
    recCG_(cnode->getL(), &c1, &t1);
    recCG_(cnode->getR(), &c2, &t2);
    recCGMult_(&t1, &t2, c1, c2, terms, c);
    break;
  case (OpNum):
    *c = cnode->getVal();
    break;
  case (OpPlus):
    recCG_(cnode->getL(), &c1, &t1);
    recCG_(cnode->getR(), &c2, &t2);
    *c = c1+c2;
    terms->insert(terms->end(), t1.begin(), t1.end());
    terms->insert(terms->end(), t2.begin(), t2.end());
    break;
  case (OpPowK):
    c2 = 1.0;
    for (int i=0; i<cnode->getR()->getVal(); ++i) {
      recCG_(cnode->getL(), &c1, &t1);
      terms->clear();
      *c = 0.0;
      recCGMult_(&t1, &t2, c1, c2, terms, c);
      t2.clear();
      t2.insert(t2.end(), terms->begin(), terms->end());
      c2 = *c;
    }
    break;
  case (OpSqr):
    recCG_(cnode->getL(), &c1, &t1);
    recCGMult_(&t1, &t1, c1, c2, terms, c);
    break;
  case (OpSumList):
    {
      for (CNode **cp=cnode->getListL(); cp<cnode->getListR(); ++cp) {
        recCG_(*cp, &c1, &t1);
        *c += c1;
        terms->insert(terms->end(), t1.begin(), t1.end());
        t1.clear();
      }
    }
    break;
  case (OpUMinus):
    recCG_(cnode->getL(), &c1, &t1);
    *c = -c1;
    for (MonomialConstIter it = t2.begin(); it!=t2.end(); ++it) {
      m = (*it)->clone();
      (*m) *= -1;
      terms->push_back(m);
    }
    break;
  case (OpVar):
    m = (MonomialFunPtr) new MonomialFunction(1.0, cg_->getVar(cnode), 1);
    terms->push_back(m);
    break;
  default:
    assert(!"can't create polynomial from this opcode!");
  }
}


void PolynomialFunction::recCGMult_(MonomialVector *t1, MonomialVector *t2,
                                    double c1, double c2,
                                    MonomialVector *terms, double *c)
{
  MonomialFunPtr m;
  *c = c1*c2;
  for (MonomialConstIter it = t1->begin(); it!=t1->end(); ++it) {
    for (MonomialConstIter it2 = t2->begin(); it2!= t2->end(); ++it2) {
      m = (*it)->clone();
      (*m) *= (*it2);
      terms->push_back(m);
    }
  }
  if (fabs(c2)>eTol_) {
    for (MonomialConstIter it = t1->begin(); it!=t1->end(); ++it) {
      m = (*it)->clone();
      (*m) *= c2;
      terms->push_back(m);
    }
  }
  if (fabs(c1)>eTol_) {
    for (MonomialConstIter it = t1->begin(); it!=t1->end(); ++it) {
      m = (*it)->clone();
      (*m) *= c1;
      terms->push_back(m);
    }
  }
}


double PolynomialFunction::removeConstant()
{
  double c = cb_;
  cb_ = 0.;
  return c;
}


void PolynomialFunction::removeLinear(LinearFunctionPtr lf)
{
  for (MonomialIter it = terms_.begin(); it!=terms_.end(); ) {
    if ((*it)->getDegree() == 1) {
      VarIntMapConstIterator it2 = (*it)->termsBegin();
      lf->incTerm(it2->first, (*it)->getCoeff());
      it = terms_.erase(it);
    } else {
      ++it;
    }
  }
}


void PolynomialFunction::removeQuadratic(QuadraticFunctionPtr qf)
{
  for (MonomialIter it = terms_.begin(); it!=terms_.end(); ) {
    if ((*it)->getDegree() == 2) {
      VarIntMapConstIterator it2 = (*it)->termsBegin();
      if (it2->second == 2) {
        // square term
        qf->incTerm(it2->first, it2->first, (*it)->getCoeff());
      } else if (it2->second == 1) {
        // bilinear term
        ConstVariablePtr v1 = it2->first;
        ++it2;
        qf->incTerm(v1, (it2)->first, (*it)->getCoeff());
      }
      it = terms_.erase(it);
    } else {
      ++it;
    }
  }
}


MonomialConstIter PolynomialFunction::termsBegin()
{
  return terms_.begin();
}
  

MonomialConstIter PolynomialFunction::termsEnd()
{
  return terms_.end();
}
  

void PolynomialFunction::write(std::ostream &out) const
{
  for (MonomialConstIter it = terms_.begin(); it!=terms_.end(); ++it) {
    out << " + ";
    (*it)->write(out);
  }
  out << " + (" << cb_ << ")";
}


void PolynomialFunction::operator+=(ConstMonomialFunPtr m)
{
  add(m);
}


void PolynomialFunction::operator+=(ConstPolyFunPtr p)
{
  if (p) {
    MonomialFunPtr m2;
    for (MonomialConstIter it = p->terms_.begin(); it!=p->terms_.end(); ++it) {
      m2 = (*it)->clone();
      terms_.push_back(m2);
    }
    cb_ += p->cb_;
  }
}


void PolynomialFunction::operator*=(double c)
{
  multiply(c);
}


void PolynomialFunction::operator+=(double c)
{
  cb_ += c;
}


void PolynomialFunction::operator+=(ConstLinearFunctionPtr lf)
{
  if (lf) {
    MonomialFunPtr m;
    for (VariableGroupConstIterator it=lf->termsBegin(); it!=lf->termsEnd();
        ++it) {
      m = (MonomialFunPtr) new MonomialFunction(1.);
      m->multiply(it->second, it->first, 1);
      terms_.push_back(m);
    }
  }
}


void PolynomialFunction::operator-=(ConstLinearFunctionPtr lf)
{
  if (lf) {
    MonomialFunPtr m;
    for (VariableGroupConstIterator it=lf->termsBegin(); it!=lf->termsEnd();
        ++it) {
      m = (MonomialFunPtr) new MonomialFunction(1.);
      m->multiply(-1*(it->second), it->first, 1);
      terms_.push_back(m);
    }
  }
}


void  PolynomialFunction::operator +=(ConstQuadraticFunctionPtr qf)
{
  if (qf) {
    MonomialFunPtr m;
    for (VariablePairGroupConstIterator it=qf->begin(); it!=qf->end(); ++it) {
      m = (MonomialFunPtr) new MonomialFunction(it->second);
      m->multiply(1, it->first.first, 1);
      m->multiply(1, it->first.second, 1);
      terms_.push_back(m);
    }
  }
}


void PolynomialFunction::operator*=(ConstLinearFunctionPtr lf)
{
  if (lf) {
    MonomialFunPtr m;
    MonomialVector terms2 = terms_;
    terms_.clear();
    for (VariableGroupConstIterator it=lf->termsBegin(); it!=lf->termsEnd();
        ++it) {
      for (MonomialConstIter it2 = terms2.begin(); it2!=terms2.end(); ++it2) {
        m = (*it2)->clone();
        m->multiply(it->second, it->first, 1);
        terms_.push_back(m);
      }
      if (fabs(cb_)>eTol_) {
        m = (MonomialFunPtr) new MonomialFunction(cb_);
        m->multiply(it->second, it->first, 1);
        terms_.push_back(m);
      }
    }
    terms2.clear();
  } else {
    clear_();
  }
  cb_ = 0;
}


void PolynomialFunction::operator*=(ConstQuadraticFunctionPtr qf)
{
  if (qf) {
    MonomialFunPtr m;
    MonomialVector terms2 = terms_;
    terms_.clear();
    for (VariablePairGroupConstIterator it = qf->begin(); 
        it != qf->end(); ++it) {
      for (MonomialConstIter it2 = terms2.begin(); it2!=terms2.end(); ++it2) {
        m = (*it2)->clone();
        m->multiply(it->second, it->first.first, 1);
        m->multiply(1., it->first.second, 1);
        terms_.push_back(m);
      }
      if (fabs(cb_)>eTol_) {
        m = (MonomialFunPtr) new MonomialFunction(cb_);
        m->multiply(it->second, it->first.first, 1);
        m->multiply(1., it->first.second, 1);
        terms_.push_back(m);
      }
    }
  }
  cb_ = 0;
}


void PolynomialFunction::operator*=(ConstPolyFunPtr p2)
{
  if (p2) {
    double c2 = p2->cb_;
    MonomialVector oldterms = terms_;
    MonomialFunPtr m;

    terms_.clear();
    for (MonomialConstIter it = oldterms.begin(); it!=oldterms.end(); ++it) {
      for (MonomialConstIter it2 = p2->terms_.begin(); it2!=p2->terms_.end(); 
          ++it2) {
        m = (*it)->clone();
        (*m) *= (*it2);
        terms_.push_back(m);
      }
    }

    if (fabs(c2)>eTol_) {
      for (MonomialConstIter it = oldterms.begin(); it!=oldterms.end(); ++it) {
        m = (*it)->clone();
        (*m) *= c2;
        terms_.push_back(m);
      }
    }

    if (fabs(cb_)>eTol_) {
      for (MonomialConstIter it = p2->terms_.begin(); it!=p2->terms_.end();
          ++it) {
        m = (*it)->clone();
        (*m) *= cb_;
        terms_.push_back(m);
      }
    }

    cb_ *= p2->cb_;
  } else {
    clear_();
  }
}


namespace Minotaur {
PolyFunPtr operator+(ConstPolyFunPtr p1, ConstPolyFunPtr p2)
{
  PolyFunPtr p = PolyFunPtr();  //NULL
  if (!p1 && !p2) {
    // do nothing.
  } else if (!p1) {
    p = p2->clone();
  } else if (!p2) {
    p = p1->clone();
  } else {
    p = p1->clone();
    (*p) += p2;
  }
  return p;
}


PolyFunPtr operator-(ConstPolyFunPtr p1, ConstPolyFunPtr p2)
{
  PolyFunPtr p = PolyFunPtr();  //NULL
  if (!p1 && !p2) {
    // do nothing.
  } else if (!p1) {
    p = p2->clone();
    (*p) *= -1;
  } else if (!p2) {
    p = p1->clone();
  } else {
    p = p2->clone();
    (*p) *= -1;
    (*p) += p1;
  }
  return p;
}


PolyFunPtr operator*(double c, ConstPolyFunPtr p2)
{
  PolyFunPtr p = PolyFunPtr(); // NULL
  if (p2 && fabs(c)>p2->eTol_) {
    LinearFunctionPtr lf = LinearFunctionPtr(); // NULL
    p = p2->clone();
    p->multiply(lf, c);
  }
  return p;
}


PolyFunPtr operator*(ConstQuadraticFunctionPtr q2, ConstLinearFunctionPtr l1)
{
  PolyFunPtr p = PolyFunPtr(); // NULL
  if (l1 && q2) {
    double c=0.;
    p = (PolyFunPtr) new PolynomialFunction();
    (*p) += 1.;
    (*p) *= q2;
    p->multiply(l1, c);
  } 
  return p;
}


PolyFunPtr operator*(ConstPolyFunPtr p2, ConstLinearFunctionPtr l1)
{
  PolyFunPtr p = PolyFunPtr(); // NULL
  if (l1 && p2) {
    p = p2->clone();
    (*p) *= l1;
  }
  return p;
}


PolyFunPtr operator*(ConstPolyFunPtr p1, ConstQuadraticFunctionPtr q2)
{
  PolyFunPtr p = PolyFunPtr(); // NULL
  if (p1 && q2) {
    p = p1->clone();
    (*p) *= q2;
  }
  return p;
}


PolyFunPtr operator*(ConstQuadraticFunctionPtr q1, ConstQuadraticFunctionPtr q2)
{
  PolyFunPtr p = PolyFunPtr(); // NULL
  if (q1 && q2) {
    p = (PolyFunPtr) new PolynomialFunction();
    (*p) += 1.;
    (*p) *= q1;
    (*p) *= q2;
  }
  return p;
}


PolyFunPtr operator*(ConstPolyFunPtr p1, ConstPolyFunPtr p2)
{
  PolyFunPtr p = PolyFunPtr(); // NULL
  if (p1 && p2) {
    p = p1->clone();
    (*p) *= p2;
  }
  return p;
}


} // namespace


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
