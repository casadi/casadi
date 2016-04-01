//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

#include <iostream>

#include "MinotaurConfig.h"
#include "Constraint.h"
#include "Function.h"
#include "Jacobian.h"
#include "Variable.h"

using namespace Minotaur;


Jacobian::Jacobian()
  : cons_(0),
    nz_(0)
{
}


Jacobian::Jacobian(const std::vector<ConstraintPtr> & cons, const UInt)
{
  ConstraintConstIterator c_iter;
  FunctionPtr f;

  nz_ = 0;
  cons_ = &cons;
  for (c_iter=cons_->begin(); c_iter!=cons_->end(); ++c_iter) {
    nz_ += (*c_iter)->getFunction()->getNumVars();
    (*c_iter)->getFunction()->prepJac();
  }
}


Jacobian::~Jacobian()
{
  cons_ = 0; // do not free.
}


UInt Jacobian::getNumNz()
{
  return nz_;
}


void Jacobian::fillRowColIndices(UInt *iRow, UInt *jCol)
{
  ConstraintConstIterator c_iter;
  FunctionPtr f;
  UInt r_cnt=0;
  UInt nz_cnt=0;

  for (c_iter=cons_->begin(); c_iter!=cons_->end(); ++c_iter) {
    f = (*c_iter)->getFunction();
    for (VarSetConstIter it=f->varsBegin(); it!=f->varsEnd(); ++it, ++nz_cnt) {
      iRow[nz_cnt] = r_cnt;
      jCol[nz_cnt] = (*it)->getIndex();
      //std::cout << "irow["<<nz_cnt<<"] = "<<iRow[nz_cnt]<<" jCol["<<nz_cnt<<"] = "
      // << jCol[nz_cnt] << std::endl; 
    }
    ++r_cnt;
  }
  assert( (r_cnt==0) || (r_cnt==cons_->size()) );
  assert(nz_cnt==nz_);
}


void Jacobian::fillRowColValues(const double *x, double *values, int *error)
{
  ConstraintConstIterator c_iter;
  UInt nz_cnt = 0;
  FunctionPtr f;

  *error = 0;
  std::fill(values, values+nz_, 0.0);
  for (c_iter=cons_->begin(); c_iter!=cons_->end(); ++c_iter) {
    f = (*c_iter)->getFunction();
    f->fillJac(x, values+nz_cnt, error);
    nz_cnt += f->getNumVars();
    if (*error != 0) {
      return;
    }
  }
  assert(nz_cnt==nz_);
}


void Jacobian::write(std::ostream &out) const
{
  out << "nz_ = " << nz_ << std::endl;
  out << "number of constraints = " << cons_->size() << std::endl;
  //exit(0);
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
