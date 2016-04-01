// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

/**
 * \file CutInfo.cpp
 * \brief Implement the methods of Cut class. 
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iostream>

#include "MinotaurConfig.h"
#include "Function.h"
#include "CutInfo.h"
#include "LinearFunction.h"
#include "Constraint.h"
#include "Problem.h"

using namespace Minotaur;

CutInfo::CutInfo()
  : cons_(ConstraintPtr()),
    cntSinceActive_(0),
    cntSinceViol_(0),
	 numActive_(0),
	 fixedScore_(0.0),
	 hashVal_(0.0)
{
}

CutInfo::CutInfo(ConstraintPtr cons, const double* hashval)
  : cons_(cons),
    cntSinceActive_(0),
    cntSinceViol_(0),
	 numActive_(0),
	 fixedScore_(0.0),
	 hashVal_(*hashval)
{
  evalFixedScore();
}


void CutInfo::write(std::ostream &out) const
{
  if (cons_) {
    cons_->write(out);
	 out << "cntSinceActive_ = " << cntSinceActive_;
	 out << ", cntSinceViol_ = " << cntSinceViol_;
	 out << ", numActive_ = " << numActive_;
	 out << ", fixedScore_ = " << fixedScore_;
	 out << ", hashVal_ = " << hashVal_;
	 out << std::endl;
  }
  else {
     out << "Empty cut" << std::endl;
  }
}

void CutInfo::evalFixedScore() {
  /// TODO: create a score that reflects sparsity and other numerical properties.
  fixedScore_ = 1.0;
}

void CutInfo::evalScore(const double *x, double &vio, double &score) 
{
  int error;
  vio = 0.0;
  score = 0.0;
  if (cons_) {
    double act = cons_->getActivity(x, &error);		
    if (cons_->getUb() < INFINITY) vio = std::max(0.0,act - cons_->getUb());  
    if (cons_->getLb() > -INFINITY) vio = std::max(vio, cons_->getLb() - act);
    score = fixedScore_ + vio;
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
