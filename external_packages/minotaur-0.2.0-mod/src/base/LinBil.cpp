// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

/**
 * \file LinBil.cpp
 * \brief Implement routines to store and modify linear constraints obtained
 * by relexing bilinear constraints of the form
 * \f$ y = x_1x_2 \f$,
 * \author Ashutosh Mahajan, IIT Bombay
 */


#include <cmath>
#include <iostream>
#include <iomanip>

#include "MinotaurConfig.h"
#include "Constraint.h"
#include "LinBil.h"
#include "Variable.h"

using namespace Minotaur;

LinBil::LinBil(VariablePtr x0, VariablePtr x1, VariablePtr y)
 : aTol_(1e-5),
   rTol_(1e-4),
   y_(y)
{
  if (x0->getIndex()>x1->getIndex()) {
    x0_ = x1;
    x1_ = x0;
  } else {
    x0_ = x0;
    x1_ = x1;
  }
  c0_ = c1_ = c2_ = c3_ = ConstraintPtr(); // NULL
}


LinBil::~LinBil() 
{
  x0_.reset();
  x1_.reset();
  y_.reset();
  c0_.reset();
  c1_.reset();
  c2_.reset();
  c3_.reset();
}


bool Minotaur::CompareLinBil::operator()(LinBil* b0, LinBil* b1) const
{
  UInt b0x0 = b0->getX0()->getId();
  UInt b0x1 = b0->getX1()->getId();

  UInt b1x0 = b1->getX0()->getId();
  UInt b1x1 = b1->getX1()->getId();

  if (b0x0 == b1x0) {
    return (b0x1 < b1x1);
  }
  return (b0x0 < b1x0);
}


bool LinBil::isViolated(const double *x, double &vio) const
{

  double xval = x[x0_->getIndex()] * x[x1_->getIndex()];
  double yval = x[y_->getIndex()];

  vio = fabs(xval - yval);
  if (vio > aTol_ && vio > fabs(yval)*rTol_) {
    return true;
  }
  return false;
}


bool LinBil::isViolated(const double x0val, const double x1val, 
                        const double yval) const
{
  double xval = x1val*x0val;
  if (fabs(xval - yval) > aTol_ && fabs(xval - yval) > fabs(yval)*rTol_) {
    return true;
  }
  return false;
}


VariablePtr LinBil::getOtherX(ConstVariablePtr x) const
{
  if (x0_==x) {
   return x1_;
  } else if (x1_==x) {
    return x0_;
  } else {
    VariablePtr v = VariablePtr(); // NULL
    return v;
  }
}


void LinBil::setCons(ConstraintPtr c0, ConstraintPtr c1, ConstraintPtr c2,
                     ConstraintPtr c3)
{
  c0_ = c0;
  c1_ = c1;
  c2_ = c2;
  c3_ = c3;
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
