// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file Solution.h
 * \brief Implement base class Solution.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 * 
 * Implement the base class Solution. 
 */

#include <cmath>
#include <iomanip>
#include <iostream>

#include "MinotaurConfig.h"
#include "Constraint.h"
#include "Problem.h"
#include "Solution.h"
#include "Variable.h"

using namespace Minotaur;

Solution::Solution()
: n_(0),
  m_(0),
  problem_(ProblemPtr()),
  x_(0),
  dualCons_(0),
  dualX_(0),
  consViol_(INFINITY),
  objValue_(INFINITY),
  comple_(INFINITY)
{
}


Solution::Solution(double obj_value, const double *x, ProblemPtr problem)
: n_(problem->getNumVars()),
  m_(problem->getNumCons()),
  problem_(problem),
  dualCons_(0),
  dualX_(0),
  consViol_(INFINITY),
  objValue_(obj_value)
{
  if (x) {
    x_ = new double[n_];
    std::copy(x, x+n_, x_);
  } else {
    x_ = 0;
  }
}


Solution::Solution(double objval, const DoubleVector &x, ProblemPtr p)
: n_(x.size()),
  m_(p->getNumCons()),
  problem_(p),
  dualCons_(0),
  dualX_(0),
  consViol_(INFINITY),
  objValue_(objval)
{
  if (x.size() > 0) {
    x_ = new double[n_];
    std::copy(x.begin(), x.end(), x_);
  } else {
    x_ = 0;
  }
}


Solution::~Solution()
{
  if (x_) {
    delete [] x_;
  }
  if (dualCons_) {
    delete [] dualCons_;
  }
  if (dualX_) {
    delete [] dualX_;
  }
}


Solution::Solution(ConstSolutionPtr sol)
{
  n_ = sol->n_;
  m_ = sol->m_;

  if (sol->x_) {
    x_ = new double[n_];
    std::copy(sol->x_, sol->x_+n_, x_);
  } else {
    x_ = 0;
  }

  if (sol->dualCons_) {
    dualCons_ = new double[m_];
    std::copy(sol->dualCons_, sol->dualCons_+m_, dualCons_);
  } else {
    dualCons_ = 0;
  }

  if (sol->dualX_) {
    dualX_ = new double[n_];
    std::copy(sol->dualX_, sol->dualX_+n_, dualX_);
  } else {
    dualX_ = 0;
  }
  consViol_ = sol->consViol_;
  objValue_ = sol->objValue_;
  comple_   = sol->comple_;
}


void Solution::setPrimal(const double *x)
{
  if (!x_) {
    x_ = new double[n_];
  }
  std::copy(x, x+n_, x_);
}

void Solution::setDualOfCons(const double *dualCons)
{
  if (!dualCons_) {
    dualCons_ = new double[m_];
  }
  std::copy(dualCons, dualCons+m_, dualCons_);
}


void Solution::setDualOfVars(const double *dualX)
{
  if (!dualX_) {
    dualX_ = new double[n_];
  }
  std::copy(dualX, dualX+n_, dualX_);
}


void Solution::write(std::ostream &out) const
{
  writePrimal(out);
  writeDual(out);
}


void Solution::writePrimal(std::ostream &out, const VarVector *v) const
{
  out << "primal values:" << std::endl;
  out << std::fixed << std::setprecision(9);
  if (x_) {
    if (v) {
      UInt i=0;
      assert(v->size() == n_);
      for (VarVector::const_iterator it=v->begin(); it!=v->end(); ++it,++i) {
        out << (*it)->getName() << "   " << x_[i] << std::endl;
      }
    } else if (problem_) {
      UInt i=0;
      for (VariableConstIterator it=problem_->varsBegin(); 
          it!=problem_->varsEnd(); ++it, ++i) {
        out << (*it)->getName() << "   " << x_[i] << std::endl;
      }
    } else {
      for (UInt i=0; i<n_; ++i) {
        out << x_[i] << std::endl;
      }
    }
  }
}


void Solution::writeDual(std::ostream &out) const
{
  out << "dual values for variables:" << std::endl;
  out << std::fixed << std::setprecision(9);
  if (dualX_) {
    if (problem_) {
      UInt i=0;
      for (VariableConstIterator it=problem_->varsBegin(); 
          it!=problem_->varsEnd(); ++it, ++i) {
        out << (*it)->getName() << "   " << dualX_[i] << std::endl;
      }
    } else {
      for (UInt i=0; i<n_; ++i) {
        out << dualX_[i] << std::endl;
      }
    } 
  } 


  if (dualCons_) {
    out << "dual values for constraints:" << std::endl;
    if (problem_) {
      UInt i=0;
      for (ConstraintConstIterator it=problem_->consBegin(); 
          it!=problem_->consEnd(); ++it, ++i) {
        out << (*it)->getName() << "   " << dualCons_[i] << std::endl;
      }
    } else {
      for (UInt i=0; i<m_; ++i) {
        out << dualCons_[i] << std::endl;
      }
    }
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
