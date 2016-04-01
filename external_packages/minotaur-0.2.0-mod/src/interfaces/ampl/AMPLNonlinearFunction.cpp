//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file AMPLNonlinearFunction.cpp
 * \brief Define the AMPLNonlinearFunction class for calling evaluation and
 * gradient routines of ASL for nonlinear functions.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <iostream>

#include "MinotaurConfig.h"
#include "AMPLNonlinearFunction.h"

using namespace MINOTAUR_AMPL;

AMPLNonlinearFunction::AMPLNonlinearFunction()
  : nVars_(0), 
    amplIndex_(0), 
    myAsl_(0), 
    neg_(false),
    isInObj_(false)
{
}


AMPLNonlinearFunction::AMPLNonlinearFunction(Minotaur::UInt i, 
                                             Minotaur::UInt n, ASL* my_asl, 
                                             bool is_in_obj)
  : nVars_(n), 
    amplIndex_(i), 
    myAsl_(my_asl), 
    neg_(false),
    isInObj_(is_in_obj)
{

}


Minotaur::NonlinearFunctionPtr 
AMPLNonlinearFunction::cloneWithVars(Minotaur::VariableConstIterator, 
                                     int *err) const 
{ 
  *err = 1;
  return Minotaur::NonlinearFunctionPtr(); 
}


double AMPLNonlinearFunction::eval(const double *x, int *error) 
{
  double r;
  double *xx = const_cast<double *>(x);
  fint ferror = 0;

  if (true == isInObj_) {
    r = (*((ASL*)myAsl_)->p.Objval)((ASL*)myAsl_, amplIndex_, xx, &ferror);
  } else {
    // call the corresponding eval function from ampl.
    r = (*((ASL*)myAsl_)->p.Conival)((ASL*)myAsl_, amplIndex_, xx, &ferror);
  }
  if (neg_) {
    r *= -1.;
  }
  *error = (int) ferror;
  return r;
}


void AMPLNonlinearFunction::evalGradient (const double *x, double *grad_f,
                                          int *error) 
{
  double *xx = const_cast<double *>(x);
  fint ferror = 0;

  std::fill(grad_f, grad_f+nVars_, 0);
  if (true == isInObj_) {
    (*(myAsl_)->p.Objgrd)(myAsl_, 0, xx, grad_f, &ferror);
  } else {
    (*(myAsl_)->p.Congrd)(myAsl_, amplIndex_, xx, grad_f, &ferror);
  }
  if (neg_) {
    for (Minotaur::UInt i=0; i<nVars_; ++i) {
      grad_f[i] *= -1.;
    }
  }
  *error = (int) ferror;
}


void AMPLNonlinearFunction::evalHessian(const double, const double *, 
                                        const Minotaur::LTHessStor *,
                                        double *, int *)
{
  assert(!"can't fill hessian in a AMPL nonlinear function.");
}


void AMPLNonlinearFunction::fillHessStor(Minotaur::LTHessStor *)
{
  assert(!"can't fill hessian in a AMPL nonlinear function.");
}


void AMPLNonlinearFunction::finalHessStor(const Minotaur::LTHessStor *) 
{
  assert(!"can't fill hessian in a AMPL nonlinear function.");
}


void AMPLNonlinearFunction::fillJac(const double *, double *, int *) 
{
  assert(!"can't fill jacobian in a ampl nonlinear function.");
}


void AMPLNonlinearFunction::getVars(Minotaur::VariableSet *vars)
{
  vars->insert(vars_.begin(), vars_.end());
}


#ifdef NDEBUG
void AMPLNonlinearFunction::multiply(double )
#else
void AMPLNonlinearFunction::multiply(double c)
#endif
{
  assert(fabs(c+1.0)<1e-9); // only allowed to multiply by -1.0
  neg_ = !neg_;
}


void AMPLNonlinearFunction::prepJac(Minotaur::VarSetConstIter, 
                                    Minotaur::VarSetConstIter)
{
  assert(!"can't fill jacobian in a ampl nonlinear function.");
}


void AMPLNonlinearFunction::setVars(Minotaur::VarSetConstIterator vb, 
                                    Minotaur::VarSetConstIterator ve)
{
  vars_.insert(vb, ve);
}


void AMPLNonlinearFunction::write(std::ostream &out) const
{
  out << "nonlinear function";
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
