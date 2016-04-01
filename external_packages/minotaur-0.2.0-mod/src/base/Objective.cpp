//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file Objective.cpp
 * \brief Define the Objective class for storing and manipulating objective
 * function.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iostream>

#include "MinotaurConfig.h"
#include "Function.h"
#include "LinearFunction.h"
#include "NonlinearFunction.h"
#include "Objective.h"
#include "QuadraticFunction.h"

using namespace Minotaur;

Objective::Objective() 
{
}


Objective::Objective(double cb, ObjectiveType otyp)
  : cb_(cb), 
    otyp_(otyp), 
    state_(NormalObj) 
{
}


Objective::Objective(FunctionPtr fPtr, double cb, ObjectiveType otyp)
  : cb_(cb),
    otyp_(otyp),
    state_(NormalObj),
    f_(fPtr)
{
}


Objective::Objective(FunctionPtr fPtr, double cb, ObjectiveType otyp, 
    std::string name)
  : cb_(cb),
    otyp_(otyp),
    state_(NormalObj),
    f_(fPtr),
    name_(name)
{
}


const LinearFunctionPtr Objective::getLinearFunction() const 
{ 
  if (f_) {
    return f_->getLinearFunction(); 
  } else {
    return LinearFunctionPtr(); // NULL
  }
}


const QuadraticFunctionPtr Objective::getQuadraticFunction() const 
{ 
  if (f_) {
    return f_->getQuadraticFunction(); 
  } else {
    return QuadraticFunctionPtr(); // NULL
  }
}


const NonlinearFunctionPtr Objective::getNonlinearFunction() const 
{ 
  if (f_) {
    return f_->getNonlinearFunction(); 
  } else {
    return NonlinearFunctionPtr(); // NULL
  }
}


double Objective::eval(const double *x, int *err) const
{
  if (f_) {
    return (cb_ + f_->eval(x, err));
  } else {
    return cb_;
  }
}


void Objective::evalGradient(const double *x, double *grad_f, int *error)
{
  // first zero out everything
  //std::fill(grad_f, grad_f+indices.size(), 0);
  if (f_) {
    f_->evalGradient(x, grad_f, error);
  }

}


void Objective::negate_()
{
  cb_ = cb_*-1.0;
  if (f_) {
    (*f_) *= -1.0;
  }
  if (otyp_ == Minimize) {
    otyp_ = Maximize;
  } else {
    otyp_ = Minimize;
  }
}


void Objective::add_(ConstLinearFunctionPtr lPtr)
{
  if (f_) {
    f_->add(lPtr);
  }
}


void Objective::add_(double cb)
{
  cb_ += cb;
}


QuadraticFunctionPtr Objective::removeQuadratic_()
{
  QuadraticFunctionPtr qf = QuadraticFunctionPtr(); // NULL
  if (f_) {
    qf = f_->removeQuadratic();
  }
  return qf;
}


void Objective::delFixedVar_(VariablePtr v, double val)
{
  if (f_ && f_->hasVar(v)) {
    double offset = f_->getFixVarOffset(v, val);
    // v->write(std::cout);
    cb_ += offset;
    f_->removeVar(v, val);
  }
}


void Objective::subst_(VariablePtr out, VariablePtr in, double rat)
{
  if (f_ && f_->hasVar(out)) {
    f_->subst(out, in, rat);
  }
}


ObjectiveType Objective::getObjectiveType() const
{
  return otyp_;
}


FunctionType Objective::getFunctionType() const 
{ 
  if (f_) {
    return f_->getType(); 
  } else {
    return Constant;
  }
}


const std::string Objective::getName() const 
{
  return name_;
}


void Objective::write(std::ostream &out) const 
{

  switch (otyp_) {
   case Minimize:
     out << "minimize ";
     break;
   case Maximize:
     out << "maximize ";
     break;
   default:
     out << "obj-sense-undefined ";
     break;
  }
  out << name_ << ": ";

  if (f_) {
    f_->write(out);
  }

  if (cb_ > 0) {
    out << " + " << cb_;
  } else {
    out << " - " << fabs(cb_);
  }
  out << ";" << std::endl;
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
