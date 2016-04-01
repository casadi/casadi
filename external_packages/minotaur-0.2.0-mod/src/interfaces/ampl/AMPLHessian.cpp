//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 
/**
 * \file AMPLHessian.cpp
 * \brief Define the AMPLHessian class for extracting Hessian of Lagrangian 
 * from ASL.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#include <sstream>
#include <stdint.h>
#include <iostream>

#include "MinotaurConfig.h"
#include "AMPLHessian.h"
#include "AMPLInterface.h"

using namespace MINOTAUR_AMPL;


AMPLHessian::AMPLHessian(AMPLInterfacePtr iface)
  : myAsl_(iface->getAsl()),
    negObj_(false)
{
  //
  // nNz_ = sphsetup(int nobj, int ow, int y, int uptri);
  // the "-1" below means that ampl will assume there is only one objective
  // and use the multiplier value for it. Even though it may seem that '0' is
  // right value, it is not so :-(. The Sphes function must also have a '-1'. 
  //
  // returns the number of nonzeros in the sparse Hessian W of the Lagrangian
  // (if uptri = 0) or its upper triangle (if uptri = 1), and stores in
  // fields sputinfo->hrownos and sputinfo->hcolstarts a description of the
  // sparsity of W
  if (iface->getReaderType() == PFGHReader) {
    nNz_ = (*(myAsl_)->p.Sphset)(myAsl_, 0, -1, 1, 1, 1);
  }
}


AMPLHessian::~AMPLHessian()
{
}


void AMPLHessian::fillRowColIndices(Minotaur::UInt *iRow, Minotaur::UInt *jCol)
{
  int i,j,cnt;
  if (nNz_>0) {
    assert(myAsl_->i.sputinfo_);
    cnt = 0;
    for (i=0; i<myAsl_->i.n_var_; ++i) {
      for (j = myAsl_->i.sputinfo_->hcolstarts[i]; 
          j<myAsl_->i.sputinfo_->hcolstarts[i+1]; ++j) {
        // if we want the sparsity pattern of the upper-triangular. do this:
        // iRow[j] = myAsl_->i.sputinfo_->hrownos[j];
        // jCol[j] = i;
        //
        // but we want the lower triangular for ipopt. so we will exchange the
        // above value for now.
        iRow[j] = i;
        jCol[j] = myAsl_->i.sputinfo_->hrownos[j];
        ++cnt;
      }
    }
    assert (cnt== (int) nNz_);
  }
  
}


void AMPLHessian::fillRowColValues(const double *x, double obj_mult,
                                   const double *con_mult, double *values,
                                   int *error)
{
  // We don't need err_jmp flags, because error is always a valid pointer. Not
  // a null pointer.
  double *xx = const_cast<double *>(x);
  double *mm = const_cast<double *>(con_mult);
  assert(error);

  if (negObj_) {
    obj_mult*= -1.;
  }

  // recompute values and gradients before calling hessian
  (*(myAsl_)->p.Xknown)(myAsl_,xx,0);

  //sphes(real *H, int nobj, real *OW, real *Y) // no error flag?
  (*(myAsl_)->p.Sphes)(myAsl_, 0, values, -1, &obj_mult, mm);
  *error = 0;

  // recompute values and gradients before calling hessian
  myAsl_->i.x_known = 0;
}


Minotaur::UInt AMPLHessian::getNumNz() const
{
  return nNz_;
}


void AMPLHessian::negateObj()
{
  negObj_ = !negObj_;
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

