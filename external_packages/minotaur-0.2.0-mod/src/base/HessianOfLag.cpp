// 
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

// /**
// \file HessianOfLag.cpp
// \brief Implement methods for obtaining the Hessian of the Lagrangean for a
// given problem
// \author Ashutosh Mahajan, Argonne National Laboratory
//
// */

#include <cmath>
#include <iostream>


#include "MinotaurConfig.h"
#include "Constraint.h"
#include "Function.h"
#include "HessianOfLag.h"
#include "Objective.h"
#include "Problem.h"
#include "Variable.h"


using namespace Minotaur;

HessianOfLag::HessianOfLag()
: etol_(1e-12),
  obj_(FunctionPtr()),
  p_(0)  // NULL
{
  stor_.nz = 0;
  stor_.nlVars = 0;
  stor_.rows = 0;
  stor_.colQs = 0;
  stor_.cols = 0;
  stor_.starts = 0;
}


HessianOfLag::HessianOfLag(Problem *p)
: etol_(1e-12),
  obj_(FunctionPtr()),
  p_(p) // NULL
{
  if (p_->getObjective()) {
    obj_ = p_->getObjective()->getFunction();
  }
  stor_.nz = 0;
  stor_.nlVars = 0;
  stor_.rows = 0;
  stor_.colQs = 0;
  stor_.cols = 0;
  stor_.starts = 0;
  setupRowCol();
}


HessianOfLag::~HessianOfLag()
{
  if (stor_.cols) {
    delete [] stor_.cols;
    stor_.cols = 0;
    delete [] stor_.rows;
    stor_.rows = 0;
    delete [] stor_.starts;
    stor_.starts = 0;
  }
}


UInt HessianOfLag::getNumNz() const
{
  return stor_.nz;
}


void HessianOfLag::fillRowColIndices(UInt *irow, UInt *jcol)
{
  UInt vindex;
  UInt nz=0;

  for (UInt i=0; i<stor_.nlVars; ++i) {
    vindex = stor_.rows[i]->getIndex();
    for (UInt j=stor_.starts[i]; j<stor_.starts[i+1]; ++j, ++nz) {
      irow[nz] = vindex;
      jcol[nz] = stor_.cols[j];
    }
  }
}


void HessianOfLag::fillRowColValues(const double *x, double obj_mult, 
                                    const double *con_mult, double *values,
                                    int *error)
{
  UInt i=0;
  FunctionPtr f;

  std::fill(values, values+stor_.nz, 0);
  if (p_->getObjective()) {
    f = p_->getObjective()->getFunction();
    if (f) {
      if (fabs(obj_mult) > etol_) {
        f->evalHessian(obj_mult, x, &stor_, values, error);
      }
    }
  }

  for (ConstraintConstIterator c_iter=p_->consBegin(); c_iter!=p_->consEnd(); 
       ++c_iter, ++i) {
    f = (*c_iter)->getFunction();
    if (fabs(con_mult[i]) > etol_) {
      f->evalHessian(con_mult[i], x, &stor_, values, error);
    }
  }
}


void HessianOfLag::setupRowCol()
{
  UInt nz;
  VariablePtr v;
  std::deque<UInt> *indq;
  UInt *cols;
  UInt i;

  // remember, we need lower triangle.
  if (stor_.cols) {
    delete [] stor_.cols;
    stor_.cols = 0;
    delete [] stor_.rows;
    stor_.rows = 0;
    delete [] stor_.starts;
    stor_.starts = 0;
  }

  stor_.nlVars = 0;
  for(VariableConstIterator it=p_->varsBegin(); it!=p_->varsEnd(); ++it) {
    if (Linear==(*it)->getFunType() || Constant==(*it)->getFunType()) {
      continue;
    }
    ++(stor_.nlVars);
  }

  i = 0;
  stor_.rows   = new VariablePtr[stor_.nlVars];
  stor_.colQs  = new std::deque<UInt>[stor_.nlVars];
  stor_.starts = new UInt[stor_.nlVars+1];
  stor_.nz = 0;
  for(VariableConstIterator it=p_->varsBegin(); it!=p_->varsEnd(); ++it) {
    if (Linear==(*it)->getFunType() || Constant==(*it)->getFunType()) {
      continue;
    }
    stor_.rows[i] = (*it);
    ++i;
  }

  if (obj_) {
    obj_->fillHessStor(&stor_);
  }
  for (ConstraintConstIterator c_iter=p_->consBegin(); c_iter!=p_->consEnd(); 
       ++c_iter) {
    (*c_iter)->getFunction()->fillHessStor(&stor_);
  }


  nz = 0;
  for (i=0; i<stor_.nlVars; ++i) {
    stor_.starts[i] = nz;
    nz += (stor_.colQs+i)->size();
  }
  stor_.starts[i] = nz;
  stor_.nz = nz;


  stor_.cols = new UInt[nz];
  cols = stor_.cols;
  indq = stor_.colQs;
  nz = 0;
  for (i=0; i<stor_.nlVars; ++i, ++indq) {
    for (std::deque<UInt>::iterator it2=indq->begin(); it2!=indq->end(); 
         ++it2,++cols) {
      *cols = *it2;
      ++nz;
    }
  }
  assert(nz == stor_.nz);

  if (obj_) {
    obj_->finalHessStor(&stor_);
  }
  for (ConstraintConstIterator c_iter=p_->consBegin(); c_iter!=p_->consEnd(); 
       ++c_iter) {
    (*c_iter)->getFunction()->finalHessStor(&stor_);
  }
  delete [] stor_.colQs;
  stor_.colQs = 0;

}


void HessianOfLag::write(std::ostream &out) const
{
  std::string me = "HessOfLag: ";
  out << me << "nz = " << stor_.nz << std::endl
      << me << "nlvars = " << stor_.nlVars << std::endl;
  for (UInt i=0; i<stor_.nlVars; ++i) {
    out << me << " matrix row " << stor_.rows[i]->getName() << " ";
    for (UInt j=stor_.starts[i]; j<stor_.starts[i+1]; ++j) {
      out << stor_.cols[j] << " ";
    }
    out << std::endl;
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
