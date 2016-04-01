//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file Constraint.cpp
 * \brief Define the Constraint class.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

#include "MinotaurConfig.h"
#include "Constraint.h"
#include "Function.h"
#include "LinearFunction.h"
#include "NonlinearFunction.h"
#include "QuadraticFunction.h"
#include "Variable.h"

using namespace Minotaur;


Constraint::Constraint() 
  : f_(FunctionPtr()),
    id_(0),
    index_(0),
    lb_(-INFINITY),
    name_(""),
    state_(NormalCons),
    ub_(INFINITY)
{
}


Constraint::Constraint(UInt id, UInt index, FunctionPtr f, double lb, double
                       ub, std::string name)
  :f_(f),
   id_(id),
   index_(index),
   lb_(lb),
   name_(name),
   state_(NormalCons),
   ub_(ub)
{
}


Constraint::~Constraint()
{
  f_.reset();
}


void Constraint::add_(double cb)
{
  lb_-=cb;
  ub_-=cb;
}


void Constraint::changeLf_(LinearFunctionPtr lf)
{
  f_->changeLf(lf);
}


void Constraint::changeNlf_(NonlinearFunctionPtr nlf)
{
  f_->changeNlf(nlf);
}


void Constraint::delFixedVar_(VariablePtr v, double val)
{
  if (f_ && f_->hasVar(v)) {
    double offset = f_->getFixVarOffset(v, val);
    lb_ -= offset;
    ub_ -= offset;
    f_->removeVar(v, val);
  }
}


double Constraint::getActivity(const double *x, int *error) const
{
  return f_->eval(x, error);
}


FunctionType Constraint::getFunctionType() const
{ 
  return f_->getType();
}


const LinearFunctionPtr Constraint::getLinearFunction() const
{ 
  return f_->getLinearFunction();
}


const std::string Constraint::getName() const 
{
  return name_;
}


const NonlinearFunctionPtr Constraint::getNonlinearFunction() const 
{ 
  return f_->getNonlinearFunction();
}


const QuadraticFunctionPtr Constraint::getQuadraticFunction() const 
{ 
  return f_->getQuadraticFunction();
}


void Constraint::reverseSense_()
{
  double tmp = lb_;
  (*f_) *= -1;
  lb_    = -1*ub_;
  ub_    = -1*tmp;
}


void Constraint::setName_(std::string name)
{
  name_ = name;
}


void Constraint::subst_(VariablePtr out, VariablePtr in, double rat,
                        bool *instay)
{
  *instay = false;
  if (f_ && f_->hasVar(out)) {
    f_->subst(out, in, rat);
    *instay = f_->hasVar(in);
  }
}


void Constraint::write(std::ostream &out) const
{

  out << "subject to " << name_ << ": ";

  if (f_) {
    if (lb_ > -INFINITY) {
      out << lb_ << " <= ";
    } 
    f_->write(out);
    if (ub_ < INFINITY) {
      out << " <= " << ub_;
    }

    if (lb_<=-INFINITY && ub_<=-INFINITY) {
      out << lb_ << " <= ";
      f_->write(out);
      out << " <= " << ub_;
    }
  } else {
    out << "no function in this constraint" << std::endl;
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
