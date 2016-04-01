//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file AMPLJacobian.cpp
 * \brief Define the AMPLJacobian class for extracting Jacobian from ASL.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <sstream>
#include <stdint.h>
#include <iostream>

#include "MinotaurConfig.h"
#include "AMPLInterface.h"
#include "AMPLJacobian.h"

using namespace MINOTAUR_AMPL;

AMPLJacobian::AMPLJacobian(AMPLInterfacePtr iface)
  : myAsl_(iface->getAsl()),
    tmp_(0),
    tmpSize_(0)
{
  nNz_ = myAsl_->i.nzc_;
}


AMPLJacobian::~AMPLJacobian()
{
  if (tmp_) {
    delete [] tmp_;
  }
}


Minotaur::UInt AMPLJacobian::getNumNz() 
{
  return nNz_;
}


void AMPLJacobian::fillColRowIndices(Minotaur::UInt *jcol, Minotaur::UInt *irow)
{
  // iRow and jCol are each of dimension nNz_. we can visit each constraint
  // and obtain the value of cg->varno and cg->goff. 
  //
  // AMPL saves the jacobian in a 1-D array. in order to get the value at row
  // i, col j of the constraints, we need to get the value of jac(cg->goff),
  // where cg corresponds to row i and col j.
  cgrad *cg; 
  Minotaur::UInt cnt = 0;

  for (int i=0; i<myAsl_->i.n_con_; ++i) {
    for (cg = myAsl_->i.Cgrad_[i]; cg; cg = cg->next) {
      irow[cg->goff] = i;
      jcol[cg->goff] = cg->varno;
      ++cnt;
    }
  }

  assert (cnt==nNz_);
}


void AMPLJacobian::fillColRowValues(const double *x, 
                                    double *values, int *error)
{
  double *xx = const_cast<double *>(x);
  fint ferror = 0;
  (*(myAsl_)->p.Jacval)(myAsl_, xx, values, &ferror);
  if (ferror!=0) {
    *error = (int) ferror;
  }
}


void AMPLJacobian::fillRowColIndices(Minotaur::UInt *irow, Minotaur::UInt *jcol)
{

  Minotaur::UInt cnt = 0;
  cgrad *cg; 
  for (int i=0; i<myAsl_->i.n_con_; ++i) {
    for (cg = myAsl_->i.Cgrad_[i]; cg; cg = cg->next) {
      irow[cnt] = i;
      jcol[cnt] = cg->varno;
      ++cnt;
    }
  }
  assert (cnt==nNz_);
  if (tmpSize_!=nNz_) {
    delete [] tmp_;
    tmp_ = new double[nNz_];
    tmpSize_ = nNz_;
  }
}


void AMPLJacobian::fillRowColValues(const double *x, 
                                    double *values, int *error)
{
  double *xx = const_cast<double *>(x);
  cgrad *cg; 
  Minotaur::UInt cnt = 0;
  fint ferror = 0;

  (*(myAsl_)->p.Jacval)(myAsl_, xx, tmp_, &ferror);
  for (int i=0; i<myAsl_->i.n_con_; ++i) {
    for (cg = myAsl_->i.Cgrad_[i]; cg; cg = cg->next) {
      values[cnt] = tmp_[cg->goff];
      ++cnt;
    }
  }
  if (ferror!=0) {
    *error = (int) ferror;
  }
  assert (cnt==nNz_);
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
