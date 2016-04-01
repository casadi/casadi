//
//	MINOTAUR -- It's only 1/2 bull
//
//	(C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file Cut.cpp
 * \brief Update information required by CutManager for each cut
 * \author Mahdi Hamzeei
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "MinotaurConfig.h"
#include "Constraint.h"
#include "Cut.h"
#include "Environment.h"
#include "Function.h"
#include "Logger.h"
#include "Operations.h"
#include "Option.h"
#include "ProblemSize.h"
#include "Relaxation.h"
#include "Solution.h"
#include "Variable.h"

#define SPEW1 0
using namespace Minotaur;


Cut::Cut()
  : cons_(ConstraintPtr()),
    f_(FunctionPtr()),
    lb_(-INFINITY),
    logger_(LoggerPtr()),
    n_(0),
    ub_(INFINITY)
{
  initInfo_(false, false);
}


Cut::Cut(UInt n, FunctionPtr f, double lb, double ub,
         bool never_delete, bool never_disable)
  : cons_(ConstraintPtr()),
    f_(f),
    lb_(lb),
    logger_(LoggerPtr()),
    n_(n),
    ub_(ub)
{
  initInfo_(never_delete, never_disable);
}

Cut::Cut(ProblemPtr p, FunctionPtr f, double lb, double ub,
	 bool never_delete, bool never_disable)
  : cons_(ConstraintPtr()),
    f_(f),
    lb_(lb),
    logger_(LoggerPtr()),
    n_(p->getNumVars()),
    ub_(ub)
{
  cons_ = p->newConstraint(f,lb,ub);
  initInfo_(never_delete,never_disable);
}

Cut::~Cut()
{
  f_.reset();
}

void Cut::applyToProblem(ProblemPtr p)
{
  cons_ = p->newConstraint(f_,lb_,ub_);
}

double Cut::eval(const double *x, int *err)
{
  return f_->eval(x, err);
}

void Cut::evalScore(const double *x, double *vio, double *score)
{
  int error = 0;
  *vio = 0.0;
  *score = 0.0;
  double act = eval(x, &error);		
  if ( ub_ < INFINITY) 
    *vio = act - ub_;  
  else if (lb_ > -INFINITY) 
    *vio = lb_ - act;
  *score = *vio / fixedScore_;
}

void Cut::evalFixedScore_()
{
  double *gr = new double[n_];
  
  std::fill(gr,gr+n_,0.0);
  f_->getLinearFunction()->evalGradient(gr);
  fixedScore_ = sqrt(InnerProduct(gr,gr,n_));
  
}

void Cut::initInfo_(bool never_delete, bool never_disable)
{
  info_.timesEnabled = 0;
  info_.timesDisabled = 0;
  info_.lastEnabled = 0;
  info_.lastDisabled = 0;
  info_.cntSinceActive = 0;
  info_.cntSinceViol = 0;
  info_.numActive = 0;
  info_.parent_active_cnts = 0;

  info_.hash = 0;
  info_.fixedScore = 0;
  info_.varScore = 0;

  info_.neverDelete = never_delete;
  info_.neverDisable = never_disable;
  info_.inRel = false;
}


void Cut::write(std::ostream &out) const
{
  out << lb_ << " <= ";
  if (f_) {
    f_->write(out);
  } else {
    out << 0.0;
  }
  out << " <= " << ub_ << std::endl;
}


void Cut::writeStats(std::ostream &out) const
{
  out << "timesEnabled   = " << info_.timesEnabled   << std::endl
      << "timesDisabled  = " << info_.timesDisabled  << std::endl
      << "lastEnabled    = " << info_.lastEnabled    << std::endl
      << "lastDisabled   = " << info_.lastDisabled   << std::endl
      << "cntSinceActive = " << info_.cntSinceActive << std::endl
      << "cntSinceViol   = " << info_.cntSinceViol   << std::endl
      << "numActive      = " << info_.numActive      << std::endl

      << "hash           = " << info_.hash           << std::endl
      << "fixedScore     = " << info_.fixedScore     << std::endl
      << "varScore       = " << info_.varScore       << std::endl

      << "neverDelete    = " << info_.neverDelete    << std::endl
      << "neverDisable   = " << info_.neverDisable   << std::endl
  ;
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
