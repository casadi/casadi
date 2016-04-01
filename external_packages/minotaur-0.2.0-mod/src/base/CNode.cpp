//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//


/**
 * \file CGraph.cpp
 * \brief Define class CNode for representing a node of computational graph 
 * of a nonlinear function.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cerrno>
#include <cmath>
#include <iostream>

#include "MinotaurConfig.h"
#include "CNode.h"
#include "Operations.h"
#include "Variable.h"

#define PI 3.141592653589793
#define MINFTY 1e25


using namespace Minotaur;

CNode::CNode()
  : b_(false),
    child_(0),
    d_(0),
    fType_(UnknownFunction),
    g_(0),
    gi_(0),
    h_(0),
    i_(0),
    id_(0),
    l_(0),
    lb_(-INFINITY),
    numChild_(0),
    numPar_(0),
    op_(OpNone),
    parB_(0),
    parE_(0),
    r_(0),
    ti_(0),
    ub_(INFINITY),
    uPar_(0),
    v_(0),
    val_(0)
{
}


CNode::CNode(OpCode op, CNode *lchild, CNode *rchild)
  : b_(false),
    child_(0),
    d_(0),
    fType_(UnknownFunction),
    g_(0),
    gi_(0),
    h_(0),
    i_(0),
    id_(0),
    l_(lchild),
    lb_(-INFINITY),
    numChild_(2),
    numPar_(0),
    op_(op),
    parB_(0),
    parE_(0),
    r_(rchild),
    ti_(0),
    ub_(INFINITY),
    uPar_(0),
    v_(0),
    val_(0)
{
  if (rchild) {
    rchild->addPar(this);
  } else {
    --numChild_;
  }
  if (lchild) {
    lchild->addPar(this);
  } else {
    --numChild_;
  }
  findFType();
}


CNode::CNode(OpCode op, CNode **children, UInt num_child)
  : child_(0),
    d_(0),
    fType_(UnknownFunction),
    g_(0),
    gi_(0),
    h_(0),
    i_(0),
    id_(0),
    l_(0),
    lb_(-INFINITY),
    numChild_(num_child),
    numPar_(0),
    op_(op),
    parB_(0),
    parE_(0),
    r_(0),
    ti_(0),
    ub_(-INFINITY),
    uPar_(0),
    v_(0),
    val_(0)
{
  if (0<num_child) {
    l_ = *children;
    if (1<num_child) {
      r_ = *(children+num_child-1);
    }
    if (2<num_child || OpSumList==op_) {
      child_ = new CNode*[num_child+1];
      for (UInt i=0; i<num_child; ++i) {
        child_[i] = *(children+i);
        child_[i]->addPar(this);
      }
      child_[numChild_] = 0;
    } else {
      if (l_) {
        l_->addPar(this);
      }
      if (r_) {
        r_->addPar(this);
      }
    }
    findFType();
  } else if (OpNum==op_) { 
    fType_ = Constant;
  } else if (OpInt==op_) {
    fType_ = Constant;
  } else if (OpVar==op_) {
    fType_ = Linear;
  } else {
    assert(!"unknown function. Can't proceed further!");
  }
}


CNode::~CNode()
{
  CQIter2 *it;
  if (child_) {
    delete [] child_;
    child_ = 0;
  }
  l_ = 0;
  while (parB_) {
    it = parB_;
    parB_ = parB_->next;
    delete it;
  }
  r_ = 0;
  v_ = 0;
}


void CNode::addPar(CNode *par)
{
  if (0==numPar_) {
    numPar_ = 1;
    uPar_ = par;
  } else if (1==numPar_) {
    parB_ = new CQIter2();
    parB_->prev = 0;
    parB_->node = uPar_;
    uPar_ = 0;

    parE_       = new CQIter2();
    parE_->prev = parB_;
    parE_->next = 0;
    parE_->node = par;
    parB_->next = parE_;
    numPar_     = 2;
  } else {
    CQIter2 *it  = new CQIter2();
    parE_->next  = it;
    it->prev     = parE_;
    it->next     = 0;
    it->node     = par;
    parE_        = it;
    ++numPar_;
  }
}


CNode* CNode::clone() const
{
  CNode *node = 0;

  node = new CNode();
  node->b_ = b_;
  node->d_ = d_;
  node->fType_ = fType_;
  node->g_ = g_;
  node->gi_ = gi_;
  node->h_ = h_;
  node->i_ = i_;
  node->id_ = id_;
  node->lb_ = lb_;
  node->numChild_ = numChild_;
  node->numPar_ = numPar_;
  node->op_ = op_;
  node->ti_ = ti_;
  node->ub_ = ub_;
  node->val_ = val_;
  // do not copy par_, child_, l_, r_ etc.
  return node;
}


void CNode::changeChild(CNode *out, CNode *in)
{
  CNode **c;
  UInt i;

  switch (numChild_) {
  case 1:
    assert(l_==out);
    l_ = in;
    break;
  case 2:
    assert(l_==out || r_==out);
    if (l_==out) {
      l_ = in;
    } else if (r_==out) {
      r_ = in;
    }
    break;
  default:
    c = child_;
    for (i=0; i<numChild_; ++i, ++c) {
      if ((*c)==out) {
        *c = in;
        break;
      }
    }
    assert(i<numChild_);
  }
}


void CNode::copyParChild(CNode *out,
                         std::map<const CNode*, CNode*> *nmap) const
{
  std::map<const CNode*, CNode*>::iterator mit;
  out->numPar_ = numPar_;
  if (uPar_) {
    out->uPar_   = nmap->find(uPar_)->second;
  }
  if (parB_) {
    CQIter2 *it = 0;
    CQIter2 *it2 = 0;
    out->parB_ = new CQIter2();
    out->parB_->prev = 0;
    out->parB_->node = nmap->find(parB_->node)->second;
    out->parB_->next = 0;
    out->parE_ = out->parB_;
    it = parB_->next;
    while (it) {
      it2 = out->parE_;
      out->parE_ = new CQIter2();
      out->parE_->prev = it2;
      out->parE_->node = nmap->find(it->node)->second;
      out->parE_->next = 0;
      it2->next = out->parE_;
      it = it->next;
    }
  }

  if (child_) {
    out->child_ = new CNode*[numChild_+1];
    for (UInt i=0; i<numChild_; ++i) {
      mit = nmap->find(child_[i]);
      assert (mit!=nmap->end());
      out->child_[i] = mit->second;
    }
    out->child_[numChild_] = 0;
  }

  if (l_) {
    mit = nmap->find(l_);
    assert (mit!=nmap->end());
    out->l_ = mit->second;
  }

  if (r_) {
    mit = nmap->find(r_);
    assert (mit!=nmap->end());
    out->r_ = mit->second;
  }
}


double CNode::eval(double x, int *error) const
{
  errno = 0; //declared in cerrno
  double val = 0;
  switch (op_) {
  case (OpAbs):
    val = fabs(x);
    break;
  case (OpAcos):
    val = acos(x);
    break;
  case (OpAcosh):
    val = acosh(x);
    break;
  case (OpAsin):
    val = asin(x);
    break;
  case (OpAsinh):
    val = asinh(x);
    break;
  case (OpAtan):
    val = atan(x);
    break;
  case (OpAtanh):
    val = atanh(x);
    break;
  case (OpCeil):
    val = ceil(x);
    break;
  case (OpCos):
    val = cos(x);
    break;
  case (OpCosh):
    val = cosh(x);
    break;
  case (OpCPow):
    val = pow(l_->val_, x);
    break;
  case (OpExp):
    val = exp(x);
    break;
  case (OpFloor):
    val = floor(x);
    break;
  case (OpInt):
    break;
  case (OpLog):
    val = log(x);
    break;
  case (OpLog10):
    val = log10(x);
    break;
  case (OpNone):
    break;
  case (OpNum):
    break;
  case (OpPowK):
    val = pow(x, r_->val_);
    break;
  case (OpRound):
    val = floor(x+0.5);
    break;
  case (OpSin):
    val = sin(x);
    break;
  case (OpSinh):
    val = sinh(x);
    break;
  case (OpSqr):
    val = x*x;
    break;
  case (OpSqrt):
    val = sqrt(x);
    break;
  case (OpTan):
    val = tan(x);
    break;
  case (OpTanh):
    val = tanh(x);
    break;
  case (OpUMinus):
    val = -(x);
    break;
  default:
    *error = 1;
  }
  if (errno!=0) {
    *error = errno;
  }
  return val;
  //std::cout << "value = " << val_ << std::endl;
}


void CNode::eval(const double *x, int *error)
{
  errno = 0; //declared in cerrno
  //writeSubExp(std::cout);
  //std::cout << "\n";
  switch (op_) {
  case (OpAbs):
    val_ = fabs(l_->val_);
    break;
  case (OpAcos):
    val_ = acos(l_->val_);
    break;
  case (OpAcosh):
    val_ = acosh(l_->val_);
    break;
  case (OpAsin):
    val_ = asin(l_->val_);
    break;
  case (OpAsinh):
    val_ = asinh(l_->val_);
    break;
  case (OpAtan):
    val_ = atan(l_->val_);
    break;
  case (OpAtanh):
    val_ = atanh(l_->val_);
    break;
  case (OpCeil):
    val_ = ceil(l_->val_);
    break;
  case (OpCos):
    val_ = cos(l_->val_);
    break;
  case (OpCosh):
    val_ = cosh(l_->val_);
    break;
  case (OpCPow):
    val_ = pow(l_->val_, r_->val_);
    break;
  case (OpDiv):
    val_ = l_->val_/r_->val_;
    break;
  case (OpExp):
    val_ = exp(l_->val_);
    break;
  case (OpFloor):
    val_ = floor(l_->val_);
    break;
  case (OpInt):
    break;
  case (OpIntDiv):
    // always round towards zero
    val_ = l_->val_/r_->val_;
    if (val_>0) {
      val_ = floor(val_); 
    } else {
      val_ = ceil(val_);
    }
    break;
  case (OpLog):
    val_ = log(l_->val_);
    break;
  case (OpLog10):
    val_ = log10(l_->val_);
    break;
  case (OpMinus):
    val_ = l_->val_ - r_->val_;
    break;
  case (OpMult):
    val_ = l_->val_ * r_->val_;
    break;
  case (OpNone):
    break;
  case (OpNum):
    break;
  case (OpPlus):
    val_ = l_->val_ + r_->val_;
    break;
  case (OpPow):
    val_ = pow(l_->val_, r_->val_);
    break;
  case (OpPowK):
    val_ = pow(l_->val_, r_->val_);
    break;
  case (OpRound):
    val_ = floor(l_->val_+0.5);
    break;
  case (OpSin):
    val_ = sin(l_->val_);
    break;
  case (OpSinh):
    val_ = sinh(l_->val_);
    break;
  case (OpSqr):
    val_ = l_->val_*l_->val_;
    break;
  case (OpSqrt):
    val_ = sqrt(l_->val_);
    break;
  case (OpSumList):
    {
    CNode **c = child_;
    val_ = 0;
    for (UInt i=0; i<numChild_; ++i, ++c) {
      val_ += (*c)->val_;
    }
    }
    break;
  case (OpTan):
    val_ = tan(l_->val_);
    break;
  case (OpTanh):
    val_ = tanh(l_->val_);
    break;
  case (OpUMinus):
    val_ = -(l_->val_);
    break;
  case (OpVar):
    val_ = x[v_->getIndex()];
    break;
  default:
    assert(!"cannot evaluate!");
  }
  if (errno!=0) {
    *error = errno;
  }
  //std::cout << "value = " << val_ << std::endl;
}


FunctionType CNode::findFType()
{
  errno = 0; // declared in cerrno
  FunctionType f0;
  fType_ = UnknownFunction;
  switch (op_) {
  case (OpAbs):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = fabs(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpAcos):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = acos(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpAcosh):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = acosh(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpAsin):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = asin(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpAsinh):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = asinh(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpAtan):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = atan(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpAtanh):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = atanh(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpCeil):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = ceil(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpCos):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = cos(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpCosh):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = cosh(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpCPow):
    if (Constant==r_->fType_) {
      fType_ = Constant;
      val_ = pow(l_->val_, r_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpDiv):
    if (Constant==l_->fType_ && Constant==r_->fType_) {
      fType_ = Constant;
      val_ = l_->val_/r_->val_;
    } else if (Constant==r_->fType_) {
      fType_ = l_->fType_;
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpExp):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = exp(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpFloor):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = floor(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpInt):
    fType_ = Constant;
    break;
  case (OpIntDiv):
    if (Constant==l_->fType_ && Constant==r_->fType_) {
      fType_ = Constant;
      val_ = (l_->val_/r_->val_);
      if (val_ < 0) {
        val_ = ceil(val_);
      } else {
        val_ = floor(val_);
      }
    } else if (Constant==r_->fType_) {
      fType_ = l_->fType_;
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpLog):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = log(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpLog10):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = log10(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpMinus):
    fType_ = funcTypesAdd(l_->fType_, r_->fType_);
    if (Constant==fType_) {
      val_ = l_->val_ - r_->val_;
    }
    break;
  case (OpMult):
    fType_ = funcTypesMult(l_->fType_, r_->fType_);
    if (Constant==fType_) {
      val_ = l_->val_ * r_->val_;
    }
    break;
  case (OpNone):
  case (OpNum):
    fType_ = Constant;
    break;
  case (OpPlus):
    fType_ = funcTypesAdd(l_->fType_, r_->fType_);
    if (Constant==fType_) {
      val_ = l_->val_ + r_->val_;
    }
    break;
  case (OpPow):
    fType_ = Nonlinear;
    break;
  case (OpPowK):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = pow(l_->val_,r_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpRound):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = floor(l_->val_+0.5);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpSumList):
    f0 = Constant;
    for (UInt i=0; i<numChild_; ++i) {
      f0 = funcTypesAdd(f0, child_[i]->fType_);
    }
    fType_ = f0;
    if (Constant==fType_) {
      val_ = 0.0;
      for (UInt i=0; i<numChild_; ++i) {
        val_ += child_[i]->val_;
      }
    }
    break;
  case (OpSqr):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = l_->val_ * l_->val_;
    } else if (Linear==l_->fType_) {
      fType_ = Quadratic;
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpSqrt):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = sqrt(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpSin):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = sin(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpSinh):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = sinh(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpTan):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = tan(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpTanh):
    if (Constant==l_->fType_) {
      fType_ = Constant;
      val_ = tanh(l_->val_);
    } else {
      fType_ = Nonlinear;
    }
    break;
  case (OpUMinus):
    fType_ = l_->fType_;
    if (Constant==fType_) {
      val_ = -(l_->val_);
    }
    break;
  case (OpVar):
    fType_ = Linear;
    break;
  default:
    fType_ = UnknownFunction;
  }
  assert(errno==0);
  return fType_;
}


void CNode::fwdGrad()
{
  switch (op_) {
  case (OpAbs):
    if (l_->val_>1e-10) {
      gi_ += 1.0;
    } else if (l_->val_<-1e-10) {
      gi_ -= 1.0;
    } else {
      gi_ += 0.0;
    }
    break;
  case (OpAcos):
    gi_ -= l_->gi_/sqrt(1-l_->val_*l_->val_); // -1/sqrt(1-x^2)
    break;
  case (OpAcosh):
    gi_ += l_->gi_/sqrt(l_->val_*l_->val_ - 1.0); // 1/sqrt(x^2-1)
    break;
  case (OpAsin):
    gi_ += l_->gi_/sqrt(1-l_->val_*l_->val_); // 1/sqrt(1-x^2)
    break;
  case (OpAsinh):
    gi_ += l_->gi_/sqrt(l_->val_*l_->val_ + 1.0); // 1/sqrt(x^2+1)
    break;
  case (OpAtan):
    gi_ += l_->gi_/(1+l_->val_*l_->val_); // 1/(1+x^2)
    break;
  case (OpAtanh):
    gi_ += l_->gi_/(1-l_->val_*l_->val_); // 1/(1-x^2)
    break;
  case (OpCeil):
    if (fabs(l_->val_ - floor(0.5+l_->val_))<1e-12) {
      gi_ += l_->gi_;
    } else {
      gi_ += 0.0;
    }
    break;
  case (OpCos):
    gi_ -= l_->gi_*sin(l_->val_);
    break;
  case (OpCosh):
    gi_ += l_->gi_*sinh(l_->val_);
    break;
  case (OpCPow):
    gi_ += r_->gi_*log(l_->val_)*val_; // val_ = a^(r_->val_)
    break;
  case (OpDiv):
    gi_ += l_->gi_/r_->val_;
    gi_ -= r_->gi_*l_->val_/(r_->val_*r_->val_);
    break;
  case (OpExp):
    gi_ += l_->gi_*val_; // val_ = e^(l_val_)
    break;
  case (OpFloor):
    gi_ += l_->gi_; // val_ = e^(l_val_)
    break;
  case (OpInt):
    break;
  case (OpIntDiv):
    assert(!"derivative of OpIntDiv not implemented!");
    break;
  case (OpLog):
    gi_ += l_->gi_/l_->val_;
    break;
  case (OpLog10):
    gi_ += l_->gi_/l_->val_/log(10.0);
    break;
  case (OpMinus):
    gi_ += l_->gi_;
    gi_ -= r_->gi_;
    break;
  case (OpMult):
    gi_ += l_->gi_*r_->val_;
    gi_ += r_->gi_*l_->val_;
    break;
  case (OpNone):
    break;
  case (OpNum):
    break;
  case (OpPlus):
    gi_ += l_->gi_;
    gi_ += r_->gi_;
    break;
  case (OpPow):
    assert(!"derivative of OpPow not implemented!");
    break;
  case (OpPowK):
    gi_ += l_->gi_*r_->val_*pow(l_->val_,r_->val_-1.0);
    break;
  case (OpRound):
    assert(!"derivative of OpRound not implemented!");
    break;
  case (OpSin):
    gi_ += l_->gi_*cos(l_->val_);
    break;
  case (OpSinh):
    gi_ += l_->gi_*cosh(l_->val_);
    break;
  case (OpSqr):
    gi_ += 2.0*l_->gi_*l_->val_; 
    break;
  case (OpSqrt):
    gi_ += l_->gi_*0.5/val_; // since val_ = sqrt(l_->val_).
    break;
  case (OpSumList):
    {
    CNode **c = child_;
    for (UInt i=0; i<numChild_; ++i, ++c) {
      gi_ += (*c)->gi_;
    }
    }
    break;
  case (OpTan):
    {
      double r = cos(l_->val_);
      gi_ += l_->gi_/(r*r);
    }
    break;
  case (OpTanh):
    { 
      double r = cosh(l_->val_);
      gi_ += l_->gi_/(r*r);
    }
    break;
  case (OpUMinus):
    gi_ -= l_->gi_;
    break;
  case (OpVar):
    break;
  default:
    break;
  }
}


void CNode::grad(int *error)
{
  errno = 0; // declared in cerrno
  switch (op_) {
  case (OpAbs):
    if (l_->val_>1e-10) {
      l_->g_ += g_;
    } else if (l_->val_<-1e-10) {
      l_->g_ -= g_;
    } else {
      l_->g_ += 0.0;
    }
    break;
  case (OpAcos):
    l_->g_ -= g_/sqrt(1-l_->val_*l_->val_); // -1/sqrt(1-x^2)
    break;
  case (OpAcosh):
    l_->g_ += g_/sqrt(l_->val_*l_->val_ - 1.0); // 1/sqrt(x^2-1)
    break;
  case (OpAsin):
    l_->g_ += g_/sqrt(1-l_->val_*l_->val_); // 1/sqrt(1-x^2)
    break;
  case (OpAsinh):
    l_->g_ += g_/sqrt(l_->val_*l_->val_ + 1.0); // 1/sqrt(x^2+1)
    break;
  case (OpAtan):
    l_->g_ += g_/(1+l_->val_*l_->val_); // 1/(1+x^2)
    break;
  case (OpAtanh):
    l_->g_ += g_/(1-l_->val_*l_->val_); // 1/(1-x^2)
    break;
  case (OpCeil):
    if (fabs(l_->val_ - floor(0.5+l_->val_))<1e-12) {
      l_->gi_ += gi_;
    } else {
      l_->gi_ += 0.0;
    }
    break;
  case (OpCos):
    l_->g_ -= g_*sin(l_->val_);
    break;
  case (OpCosh):
    l_->g_ += g_*sinh(l_->val_);
    break;
  case (OpCPow):
    r_->g_ += g_*log(l_->val_)*val_;
    break;
  case (OpDiv):
    l_->g_ += g_/r_->val_;
    r_->g_ -= g_*l_->val_/(r_->val_*r_->val_);
    break;
  case (OpExp):
    l_->g_ += g_*val_; // val_ = e^(l_val_)
    break;
  case (OpFloor):
    l_->g_ += g_; // assuming that gradient is 1.
    break;
  case (OpInt):
    break;
  case (OpIntDiv):
    assert(!"derivative of OpIntDiv not implemented!");
    break;
  case (OpLog):
    l_->g_ += g_/l_->val_;
    break;
  case (OpLog10):
    l_->g_ += g_/l_->val_/log(10.0);
    break;
  case (OpMinus):
    l_->g_ += g_;
    r_->g_ -= g_;
    break;
  case (OpMult):
    l_->g_ += g_*r_->val_;
    r_->g_ += g_*l_->val_;
    break;
  case (OpNone):
    break;
  case (OpNum):
    break;
  case (OpPlus):
    l_->g_ += g_;
    r_->g_ += g_;
    break;
  case (OpPow):
    assert(!"derivative of OpPow not implemented!");
    break;
  case (OpPowK):
    l_->g_ += g_*r_->val_*pow(l_->val_, r_->val_-1.0);
    break;
  case (OpRound):
    assert(!"derivative of OpRound not implemented!");
    break;
  case (OpSin):
    l_->g_ += g_*cos(l_->val_);
    break;
  case (OpSinh):
    l_->g_ += g_*cosh(l_->val_);
    break;
  case (OpSqr):
    l_->g_ += 2.0*g_*l_->val_; 
    break;
  case (OpSqrt):
    l_->g_ += g_*0.5/val_; // same as l_->g += g_*0.5/sqrt(l_->val_).
    break;
  case (OpSumList):
    if (g_!=0.0) {
      CNode **c = child_;
      for (UInt i=0; i<numChild_; ++i, ++c) {
        (*c)->g_ += g_;
      }
    }
    break;
  case (OpTan):
    {
      double r = cos(l_->val_);
      l_->g_ += g_/(r*r);
    }
    break;
  case (OpTanh):
    { 
      double r = cosh(l_->val_);
      l_->g_ += g_/(r*r);
    }
    break;
  case (OpUMinus):
    l_->g_ -= g_;
    break;
  case (OpVar):
    break;
  default:
    break;
  }
  if (errno != 0) {
    *error = errno;
  }
}


void CNode::hess(int *error)
{
  errno = 0;
  switch (op_) {
  case (OpAbs):
    l_->h_ += 0.0;
    break;
  case (OpAcos):
    l_->h_ += -h_/sqrt(1-l_->val_*l_->val_) 
            - g_ * l_->gi_ * l_->val_/pow((1.0-l_->val_*l_->val_),1.5) ; 
    // -x/(1-x^2)^1.5
    break;
  case (OpAcosh):
    l_->h_ +=  h_/sqrt(l_->val_*l_->val_-1.0) 
            - g_ * l_->gi_ * l_->val_/pow((l_->val_*l_->val_-1.0),1.5) ; 
    break;
  case (OpAsin):
    l_->h_ +=  h_/sqrt(1-l_->val_*l_->val_) 
            + g_ * l_->gi_ * l_->val_/pow((1-l_->val_*l_->val_),1.5) ; 
    break;
  case (OpAsinh):
    l_->h_ +=  h_/sqrt(1+l_->val_*l_->val_) 
            - g_ * l_->gi_ * l_->val_/pow((1+l_->val_*l_->val_),1.5) ; 
    break;
  case (OpAtan):
    {
    double d = 1+l_->val_*l_->val_;
    l_->h_ +=  h_/d - 2.0 * g_ * l_->gi_ * l_->val_/(d*d); 
    }
    break;
  case (OpAtanh):
    {
    double d = (1.0 - l_->val_*l_->val_);
    l_->h_ +=  h_/d + 2.0 * g_ * l_->gi_ * l_->val_/(d*d);
    }
    break;
  case (OpCeil):
    l_->h_ +=  h_;
    break;
  case (OpCos):
    l_->h_ += -h_*sin(l_->val_) - g_ * l_->gi_ * val_; // since val_ = cos(). 
    break;
  case (OpCosh):
    l_->h_ += h_*sinh(l_->val_) + g_ * l_->gi_ * val_; // since val_ = cosh().
    break;
  case (OpCPow):
    r_->h_ += h_*log(l_->val_)*val_ 
            + g_*r_->gi_*log(l_->val_)*log(l_->val_)*val_;
    break;
  case (OpDiv):
    l_->h_ += h_/r_->val_ - g_ * r_->gi_ /(r_->val_*r_->val_);
    r_->h_ += -h_*l_->val_/(r_->val_*r_->val_) 
            - g_ * l_->gi_ /(r_->val_*r_->val_)
            + g_ * r_->gi_ * l_->val_ * 2.0 /(r_->val_*r_->val_*r_->val_);
    break;
  case (OpExp):
    l_->h_ += h_*val_ + g_ * l_->gi_ * val_;
    break;
  case (OpFloor):
    l_->h_ += h_;
    break;
  case (OpInt):
    break;
  case (OpIntDiv):
    assert(!"derivative of OpIntDiv not implemented!");
    break;
  case (OpLog):
    l_->h_ += h_/l_->val_ - g_ * l_->gi_  / (l_->val_ * l_->val_); // -1/x^2
    break;
  case (OpLog10):
    l_->h_ += h_/l_->val_/log(10) - g_ * l_->gi_ / (log(10)*l_->val_*l_->val_);
    break;
  case (OpMinus):
    l_->h_ += h_;
    r_->h_ -= h_;
    break;
  case (OpMult):
    l_->h_ += h_*r_->val_ + g_*r_->gi_;
    r_->h_ += h_*l_->val_ + g_*l_->gi_;
    break;
  case (OpNone):
    break;
  case (OpNum):
    break;
  case (OpPlus):
    l_->h_ += h_;
    r_->h_ += h_;
    break;
  case (OpPow):
    assert(!"derivative of OpPow not implemented!");
    break;
  case (OpPowK):
    l_->h_ += h_ * r_->val_ * pow(l_->val_,r_->val_-1.0) 
            + g_*l_->gi_*r_->val_*(r_->val_-1.0)*pow(l_->val_, r_->val_-2.0);
    break;
  case (OpRound):
    assert(!"derivative of OpRound not implemented!");
    break;
  case (OpSin):
    l_->h_ += h_*cos(l_->val_) - g_ * l_->gi_ * val_; // since val_ = sin(). 
    break;
  case (OpSinh):
    l_->h_ += h_*cosh(l_->val_) + g_ * l_->gi_ * val_; // since val_ = sinh().
    break;
  case (OpSqr):
    l_->h_ += 2.0*h_*l_->val_ + g_ * 2.0 * l_->gi_;
    break;
  case (OpSqrt):
    l_->h_ += h_*0.5/val_ - g_ * l_->gi_ * 0.25 /(val_ * l_->val_); 
    // -1/4/x^1.5
    break;
  case (OpSumList):
    if (h_!=0.0) {
      CNode **c = child_;
      for (UInt i=0; i<numChild_; ++i, ++c) {
        (*c)->h_ += h_;
      }
    }
    break;
  case (OpTan):
    {
    double d = cos(l_->val_);
           d *= d;
    l_->h_ += h_/d  + 2.0 * g_ * l_->gi_ * tan(l_->val_)/d;
    }
    break;
  case (OpTanh):
    { 
    double d = cosh(l_->val_);
           d *= d;
    l_->h_ += h_/d  - 2.0 * g_ * l_->gi_ * tanh(l_->val_)/d;
      
    }
    break;
  case (OpUMinus):
    l_->h_ -= h_;
    break;
  case (OpVar):
    break;
  default:
    break;
  }
  if (errno != 0) {
    *error = errno;
  }
  //writeSubExp(std::cout);
  //std::cout << "h_ = " << h_ << std::endl;
}


void CNode::hess2(CNodeRSet *nset, int *error)
{
  propHessSpa2(nset);
  hess(error);
}


UInt CNode::numChild() const
{
  return numChild_;
}


void CNode::propBounds(bool *is_inf, int *error)
{
  errno = 0; //declared in cerrno
  double lb = -INFINITY;
  double ub = INFINITY;
  switch (op_) {
  case (OpAbs):
    assert(lb_>-1e-12);
    lb = -ub_;
    ub = ub_;
    l_->propBounds_(lb, ub, is_inf);
    break;
  case (OpAcos):
    lb = -1.0;
    ub = 1.0;
    l_->propBounds_(lb, ub, is_inf);
    break;
  case (OpAcosh):
    // TODO: Implement me
    break;
  case (OpAsin):
    lb = -1.0;
    ub = 1.0;
    l_->propBounds_(lb, ub, is_inf);
    break;
  case (OpAsinh):
    // TODO: Implement me
    break;
  case (OpAtan):
    // TODO: Implement me
    break;
  case (OpAtanh):
    // TODO: Implement me
    break;
  case (OpCeil):
    lb = floor(lb_);
    ub = floor(ub_);
    l_->propBounds_(lb, ub, is_inf);
    break;
  case (OpCos):
    // TODO: Implement me
    break;
  case (OpCosh):
    // TODO: Implement me
    break;
  case (OpCPow):
    // TODO: Implement me
    break;
  case (OpDiv):
    BoundsOnProduct(r_->lb_, r_->ub_, lb_, ub_, lb, ub);
    l_->propBounds_(lb, ub, is_inf);
    BoundsOnDiv(l_->lb_, l_->ub_, lb_, ub_, lb, ub);
    r_->propBounds_(lb, ub, is_inf);
    break;
  case (OpExp):
    lb = log(lb_);
    ub = log(ub_);
    l_->propBounds_(lb, ub, is_inf);
    break;
  case (OpFloor):
    lb = ceil(lb_);
    ub = ceil(ub_);
    l_->propBounds_(lb, ub, is_inf);
    break;
  case (OpIntDiv):
    // TODO: Implement me
    break;
  case (OpLog):
    lb = exp(lb_);
    ub = exp(ub_);
    l_->propBounds_(lb, ub, is_inf);
    break;
  case (OpLog10):
    lb = pow(10.0, lb_);
    ub = pow(10.0, ub_);
    l_->propBounds_(lb, ub, is_inf);
    break;
  case (OpMinus):
    lb = lb_ + r_->lb_;
    ub = ub_ + r_->ub_;
    l_->propBounds_(lb, ub, is_inf);
    lb = l_->lb_ - ub_;
    ub = l_->ub_ - lb_;
    r_->propBounds_(lb, ub, is_inf);
    break;
  case (OpMult):
    BoundsOnDiv(lb_, ub_, r_->lb_, r_->ub_, lb, ub);
    l_->propBounds_(lb, ub, is_inf);
    BoundsOnDiv(lb_, ub_, l_->lb_, l_->ub_, lb, ub);
    r_->propBounds_(lb, ub, is_inf);
    break;
  case (OpNone):
    break;
  case (OpNum):
    break;
  case (OpPlus):
    lb = lb_-r_->ub_;
    ub = ub_-r_->lb_;
    l_->propBounds_(lb, ub, is_inf);
    lb = lb_-l_->ub_;
    ub = ub_-l_->lb_;
    r_->propBounds_(lb, ub, is_inf);
    break;
  case (OpPow):
    // TODO: Implement me
    break;
  case (OpPowK):
    if (r_->val_>0) {
      if (IsInt(r_->val_/2.0)) { 
        if (ub_<-1e-12) {
          *error = 3141;
        } else {
          ub =  pow(ub_, 1.0/r_->val_);
          lb = -ub;
          l_->propBounds_(lb, ub, is_inf);
          // std::cout << "new bounds = " << lb << " " << ub << std::endl;
        }
      } else if (IsInt((r_->val_+1)/2.0)) { 
        if (lb < 0) {
          lb = -pow(-lb_, 1.0/r_->val_);
        } else {
          lb = pow(lb_, 1.0/r_->val_);
        }
        if (ub < 0) {
          ub = -pow(-ub_, 1.0/r_->val_);
        } else {
          ub = pow(ub_, 1.0/r_->val_);
        }
        l_->propBounds_(lb, ub, is_inf);
      }
    }
    break;
  case (OpRound):
    break;
  case (OpSin):
    // TODO: Implement me
    break;
  case (OpSinh):
    // TODO: Implement me
    break;
  case (OpSqr):
    ub = sqrt(ub);
    lb = -ub;
    l_->propBounds_(lb, ub, is_inf);
    break;
  case (OpSqrt):
    if (ub_<0.0) {
      *is_inf = true;
    } else if (lb>=0.0) {
      l_->propBounds_(lb*lb, ub*ub, is_inf);
    } else {
      l_->propBounds_(0.0, ub*ub, is_inf);
    }
    break;
  case (OpSumList):
    {
      // TODO: move to a separate function
      CNode **c = child_;
      bool inf_lb = false;
      bool inf_ub = false;
      double tlb;
      double tub;

      lb = 0.0;
      for (UInt i=0; i<numChild_; ++i, ++c) {
        if ((*c)->lb_>-INFINITY) {
          lb += (*c)->lb_;
        } else if (true==inf_lb) {
          lb = -INFINITY;
          break;
        } else {
          inf_lb = true;
        }
      }

      c = child_;
      ub = 0.0;
      for (UInt i=0; i<numChild_; ++i, ++c) {
        if ((*c)->ub_ < INFINITY) {
          ub += (*c)->ub_;
        } else if (true==inf_ub) {
          ub = INFINITY;
          break;
        } else {
          inf_ub = true;
        }
      }

      if (lb > -INFINITY || ub < INFINITY) {
        c = child_;
        for (UInt i=0; i<numChild_; ++i, ++c) {
          tlb = -INFINITY;
          tub =  INFINITY;
          if (ub<INFINITY) {
            if (false == inf_ub) {
              tlb = lb_-(ub-(*c)->ub_);
            } else if ((*c)->ub_ < INFINITY) {
              tlb = -INFINITY;
            } else {
              tlb = lb_-ub;
            }
          } else {
            tlb = -INFINITY;
          }
          if (lb>-INFINITY) {
            if (false == inf_lb) {
              tub = ub_-(lb-(*c)->lb_);
            } else if ((*c)->lb_ > -INFINITY) {
              tub = INFINITY;
            } else {
              tub = ub_-lb;
            }
          } else {
            tub = -INFINITY;
          }
          (*c)->propBounds_(tlb, tub, is_inf);
          if (true==(*is_inf)) {
            break;
          }
        }
      }
    }
    break;
  case (OpTan):
    // TODO: Implement me
    break;
  case (OpTanh):
    // TODO: Implement me
    break;
  case (OpUMinus):
    lb =-ub_;
    ub =-lb_;
    l_->propBounds_(lb, ub, is_inf);
    break;
  case (OpVar):
    break;
  default:
    break;
  }
  if (errno!=0) {
    *error = errno;
  }
}


void CNode::propBounds_(double lb, double ub, bool *is_inf)
{
  double etol = 1e-7;
  assert(false==std::isnan(lb));
  assert(false==std::isnan(ub));

  if (lb<-MINFTY) {
    lb = -INFINITY;
  }
  if (ub>MINFTY) {
    ub = INFINITY;
  }
  if (lb > ub + etol || ub < lb_- etol || lb > ub_ + etol) {
    *is_inf = true;
  } else {
    if (lb > lb_) {
      lb_ = lb;
    }
    if (ub < ub_) {
      ub_ = ub;
    } 
  }
}


void CNode::propHessSpa2(CNodeRSet *nset)
{
  propHessSpa();
  switch(numChild_) {
  case(0):
    break;
  case(1):
    if ((l_->b_ || l_->ti_>0) && l_->op_!=OpNum && l_->op_!=OpInt) {
      nset->insert(l_);
    }
    break;
  case(2):
    if ((l_->b_ || l_->ti_>0) && l_->op_!=OpNum && l_->op_!=OpInt) {
      nset->insert(l_);
    }
    if ((r_->b_ || r_->ti_>0) && r_->op_!=OpNum && r_->op_!=OpInt) {
      nset->insert(r_);
    }
    break;
  default:
    {
    CNode **c = child_;
    for (UInt i=0; i<numChild_; ++i, ++c) {
      if ( ((*c)->b_  || (*c)->ti_>0) &&
           (*c)->op_!=OpNum && (*c)->op_!=OpInt) {
        nset->insert(*c);
      }
    }
    }
  }
}


void CNode::propHessSpa()
{
  if (numChild_<1) {
    //write(std::cout);
    return;
  } else if (b_) {
    switch(numChild_) {
    case(1):
      l_->b_ = true;
      break;
    case(2):
      l_->b_ = true;
      r_->b_ = true;
      break;
    default:
      {
      CNode **c = child_;
      for (UInt i=0; i<numChild_; ++i, ++c) {
        (*c)->b_ = true;
      }
      }
    }
  } else {
    // g_ is false
    switch(op_) {
    case (OpAbs):
    case (OpAcos):
    case (OpAcosh):
    case (OpAsin):
    case (OpAsinh):
    case (OpAtan):
    case (OpAtanh):
      if (l_->ti_>0) {
        l_->b_ = true;
      }
      break;
    case (OpCeil):
      break;
    case (OpCos):
    case (OpCosh):
      if (l_->ti_>0) {
        l_->b_ = true;
      }
      break;
    case (OpCPow):
      if (r_->ti_>0) {
        r_->b_ = true;
      }
      break;
    case (OpDiv):
      if (r_->ti_>0) {
        l_->b_ = true;
        r_->b_ = true;
      } else if (l_->ti_>0) {
        r_->b_ = true;
      }
      break;
    case(OpExp):
      if (l_->ti_>0) {
        l_->b_ = true;
      }
      break;
    case(OpFloor):
      break;
    case(OpInt):
      break;
    case (OpIntDiv):
      assert(!"not implemented yet!");
      break;
    case (OpLog):
    case (OpLog10):
      if (l_->ti_>0) {
        l_->b_ = true;
      }
      break;
    case (OpMinus):
      break;
    case (OpMult):
      if (r_->ti_>0) {
        l_->b_ = true;
      }
      if (l_->ti_>0) {
        r_->b_ = true;
      }
      break;
    case (OpNone):
    case (OpNum):
      break;
    case (OpPlus):
      break;
    case (OpPow):
      assert(!"not implemented yet!");
      break;
    case (OpPowK):
    case (OpRound):
    case (OpSin):
    case (OpSinh):
    case (OpSqr):
    case (OpSqrt):
      if (l_->ti_>0) {
        l_->b_ = true;
      }
      break;
    case (OpSumList):
      break;
    case (OpTan):
    case (OpTanh):
      if (l_->ti_>0) {
        l_->b_ = true;
      }
    case (OpUMinus):
      // do nothing.
      break;
    case (OpVar):
      break;
    default:
      assert(!"not dealing with this opcode!");
    }
  }
  //writeSubExp(std::cout);
  //std::cout << " g_ = " << g_ << " ti_ " << ti_;
  //std::cout << "\n";
}


void CNode::setType(FunctionType t)
{
  fType_ = t;
}


void CNode::setVal(double v)
{
  val_ = v;
  if (OpNum==op_) {
    lb_ = ub_ = val_;
  }
}


void CNode::updateBnd(int *error)
{
  errno = 0; //declared in cerrno
  switch (op_) {
  case (OpAbs):
    if (l_->ub_<0) {
      lb_ = -l_->ub_;
      ub_ = -l_->lb_;
    } else if (l_->lb_<0) {
      if (-l_->lb_>l_->ub_) {
        lb_ =  l_->ub_;
        ub_ = -l_->lb_;
      } else {
        lb_ = -l_->ub_;
        ub_ =  l_->lb_;
      }
    } else {
      lb_ = l_->lb_;
      ub_ = l_->ub_;
    }
    break;
  case (OpAcos):
    lb_ = 0.0;
    ub_ = PI;
    break;
  case (OpAcosh):
    // TODO: Implement me
    lb_ = -INFINITY;
    ub_ = INFINITY;
    break;
  case (OpAsin):
    lb_ = -PI/2;
    ub_ =  PI/2;
    break;
  case (OpAsinh):
    // TODO: Implement me
    lb_ = -INFINITY;
    ub_ = INFINITY;
    break;
  case (OpAtan):
    lb_ = -PI/2;
    ub_ =  PI/2;
    break;
  case (OpAtanh):
    // TODO: Implement me
    lb_ = -INFINITY;
    ub_ = INFINITY;
    break;
  case (OpCeil):
    lb_ = ceil(l_->lb_);
    ub_ = ceil(l_->ub_);
    break;
  case (OpCos):
    lb_ = -1.0;
    ub_ = 1.0;
    break;
  case (OpCosh):
    // TODO: Implement me
    lb_ = -INFINITY;
    ub_ = INFINITY;
    break;
  case (OpCPow):
    // TODO: Implement me
    lb_ = -INFINITY;
    ub_ = INFINITY;
    break;
  case (OpDiv):
    BoundsOnDiv(l_->lb_, l_->ub_, r_->lb_, r_->ub_, lb_, ub_);
    break;
  case (OpExp):
    if (l_->lb_==-INFINITY) {
      lb_ = 0.0;
    } else {
      lb_ = exp(l_->lb_);
    }
    if (l_->ub_==INFINITY) {
      ub_ = INFINITY;
    } else {
      ub_ = exp(l_->ub_);
    }
    break;
  case (OpFloor):
    lb_ = floor(l_->lb_);
    ub_ = floor(l_->ub_);
    break;
  case (OpIntDiv):
    // TODO: Implement me
    lb_ = -INFINITY;
    ub_ = INFINITY;
    break;
  case (OpLog):
    if (l_->lb_<=0.0) {
      lb_ = -INFINITY;
    } else {
      lb_ = log(l_->lb_);
    }
    ub_ = log(l_->ub_);
    break;
  case (OpLog10):
    if (l_->lb_<=0.0) {
      lb_ = -INFINITY;
    } else {
      lb_ = log10(l_->lb_);
    }
    ub_ = log10(l_->ub_);
    break;
  case (OpMinus):
    lb_ = l_->lb_ - r_->ub_;
    ub_ = l_->ub_ - r_->lb_;
    break;
  case (OpMult):
    BoundsOnProduct(l_->lb_, l_->ub_, r_->lb_, r_->ub_, lb_, ub_);
    break;
  case (OpNone):
    break;
  case (OpNum):
    lb_ = d_;
    ub_ = d_;
    break;
  case (OpPlus):
    lb_ = l_->lb_+r_->lb_;
    ub_ = l_->ub_+r_->ub_;
    break;
  case (OpPow):
    // TODO: Implement me
    lb_ = -INFINITY;
    ub_ = INFINITY;
    break;
  case (OpPowK):
    // TODO: Implement me
    lb_ = -INFINITY;
    ub_ = INFINITY;
    break;
  case (OpRound):
    // TODO: Implement me
    lb_ = -INFINITY;
    ub_ = INFINITY;
    break;
  case (OpSin):
    lb_ = -1.0;
    ub_ =  1.0;
    break;
  case (OpSinh):
    // TODO: Implement me
    lb_ = -INFINITY;
    ub_ = INFINITY;
    break;
  case (OpSqr):
    BoundsOnSquare(l_->lb_, l_->ub_, lb_, ub_);
    break;
  case (OpSqrt):
    if (l_->lb_<1e-12) {
      lb_ = 0.0;
    } else {
      lb_ = sqrt(l_->lb_);
    }
    ub_ = sqrt(l_->ub_);
    break;
  case (OpSumList):
    {
      lb_ = ub_ = 0.0;
      CNode **c = child_;
      for (UInt i=0; i<numChild_; ++i, ++c) {
        lb_ += (*c)->lb_;
        ub_ += (*c)->ub_;
      }
    }
    break;
  case (OpTan):
    // TODO: Implement me
    lb_ = -INFINITY;
    ub_ = INFINITY;
    break;
  case (OpTanh):
    // TODO: Implement me
    lb_ = -INFINITY;
    ub_ = INFINITY;
    break;
  case (OpUMinus):
    lb_ =-l_->ub_;
    ub_ =-l_->lb_;
    break;
  case (OpVar):
    lb_ = v_->getLb();
    ub_ = v_->getUb();
    break;
  default:
    break;
  }
  if (errno!=0) {
    *error = errno;
  }
  if (lb_<-MINFTY) {
    lb_ = -INFINITY;
  }
  if (ub_>MINFTY) {
    ub_ = INFINITY;
  }
}


void CNode::write(std::ostream &out) const
{
  out 
    << "\n"
    << "d_ = " << d_ << std::endl
    << "function type = " << fType_ << std::endl
    << "g_ = " << g_ << std::endl
    << "number of children = " << numChild_ << std::endl
    << "number of parents = " << numPar_ << std::endl
    << "opcode = " << op_ << std::endl
    << "value = " << val_ << std::endl
    << "\n";
}


void CNode::writeSubExp(std::ostream &out) const
{
  switch (op_) {
  case (OpAbs):
    out << "abs(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpAcos):
    out << "acos(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpAcosh):
    out << "acosh(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpAsin):
    out << "asin(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpAsinh):
    out << "asinh(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpAtan):
    out << "atan(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpAtanh):
    out << "atanh(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpCeil):
    out << "ceil(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpCos):
    out << "cos(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpCosh):
    out << "cosh(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpCPow):
    out << l_->val_ << "^(";
    r_->writeSubExp(out);
    out << ")";
    break;
  case (OpDiv):
    out << "(";
    l_->writeSubExp(out);
    out << "/";
    r_->writeSubExp(out);
    out << ")";
    break;
  case (OpExp):
    out << "exp(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpFloor):
    out << "floor(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpIntDiv):
    out << "intdiv(";
    l_->writeSubExp(out);
    r_->writeSubExp(out);
    out << ")";
    break;
  case (OpLog):
    out << "log(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpLog10):
    out << "log10(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpMinus):
    out << "(";
    l_->writeSubExp(out);
    out << " - ";
    r_->writeSubExp(out);
    out << ")";
    break;
  case (OpMult):
    out << "(";
    l_->writeSubExp(out);
    out << " * ";
    r_->writeSubExp(out);
    out << ")";
    break;
  case (OpNone):
    out << "OpNone()";
    break;
  case (OpNum):
    out << d_;
    break;
  case (OpPlus):
    out << "(";
    l_->writeSubExp(out);
    out << " + ";
    r_->writeSubExp(out);
    out << ")";
    break;
  case (OpPow):
    out << "(";
    l_->writeSubExp(out);
    out << ")^(";
    r_->writeSubExp(out);
    out << ")";
    break;
  case (OpPowK):
    out << "(";
    l_->writeSubExp(out);
    out << ")^" << r_->val_;
    break;
  case (OpRound):
    out << "round(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpSin):
    out << "sin(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpSinh):
    out << "sinh(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpSqr):
    out << "(";
    l_->writeSubExp(out);
    out << ")^2";
    break;
  case (OpSqrt):
    out << "sqrt(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpSumList):
    out << "(";
    for (UInt i=0; i<numChild_; ++i) {
      child_[i]->writeSubExp(out);
      if (i<numChild_-1) {
        out << " + ";
      }
    }
    out << ")";
    break;
  case (OpTan):
    out << "tan(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpTanh):
    out << "tanh(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpUMinus):
    out << "-";
    out << "(";
    l_->writeSubExp(out);
    out << ")";
    break;
  case (OpVar):
    out << v_->getName();
    break;
  default:
    break;
  }
}


bool Minotaur::CompareCNodes::operator()(const CNode* n1, const CNode *n2) const
{
  // process lower ids first and lower variables first.
  if (n1->getId()==n2->getId() && n1->getId()==0) {
    return (n1->getV()->getIndex()<n2->getV()->getIndex());
  }
  return (n1->getId()<n2->getId());
}


bool Minotaur::CompareCNodesR::operator()(const CNode* n1, const CNode *n2) const
{
  // process higher ids first, but lower variables first.
  if (n1->getId()==n2->getId() && n1->getId()==0) {
    return (n1->getV()->getIndex()<n2->getV()->getIndex());
  }
  return (n1->getId()>n2->getId());
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
