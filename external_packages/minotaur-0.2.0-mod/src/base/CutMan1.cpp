// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

/**
 * \file CutMan1.cpp
 * \brief Implement the methods of CutMan1 class. 
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iostream>

#include "MinotaurConfig.h"
#include "Environment.h"
#include "Function.h"
#include "Cut.h"
#include "CutMan1.h"
#include "LinearFunction.h"
#include "Logger.h"
#include "Problem.h"
#include "Solution.h"
#include "SolutionPool.h"
#include "Types.h"
#include "Variable.h"

using namespace Minotaur;

typedef std::list<CutPtr>::const_iterator CCIter;
typedef std::vector<ConstraintPtr>::const_iterator ConstIter;

const std::string CutMan1::me_ = "CutMan1: "; 

CutMan1::CutMan1()
  : absTol_(1e-6),
    env_(EnvPtr()),   // NULL
    maxDisCutAge_(3),
    maxInactCutAge_(1),
    p_(ProblemPtr())  // NULL
{
  logger_ = (LoggerPtr) new Logger(LogDebug2);
}


CutMan1::CutMan1(EnvPtr env, ProblemPtr p)
  : absTol_(1e-6),
    env_(env),
    maxDisCutAge_(1),
    maxInactCutAge_(1),
    p_(p)
{
  logger_ = (LoggerPtr) new Logger(LogDebug2);
}


CutMan1::~CutMan1()
{
  pool_.clear();
  enCuts_.clear();
}


UInt CutMan1::getNumCuts() const
{
  return getNumEnabledCuts() + getNumDisabledCuts() + getNumNewCuts();
}


UInt CutMan1::getNumEnabledCuts() const
{
  return enCuts_.size();
}


UInt CutMan1::getNumDisabledCuts() const
{
  return disCuts_.size();
}


UInt CutMan1::getNumNewCuts() const
{
  return newCuts_.size();
}


void CutMan1::postSolveUpdate(ConstSolutionPtr sol, EngineStatus)
{
  UInt n = p_->getNumVars();
  const double *x = new double[n];
  x = sol->getPrimal();
  CutPtr con;
  double viol = 0;
  bool del_const = false;
  CutList cpyrel;
  CutInfo *cinfo;

  int num_dels = 0;
  int err;

  for (CCIter it = enCuts_.begin(); it != enCuts_.end(); ++it) {
    con = *it;
    con->eval(x, &err);
    cinfo = con->getInfo();
    if (viol > -absTol_) {
      ++(cinfo->numActive);
      cinfo->cntSinceActive = 0;
    } else {
      ++(cinfo->numActive);
    }
    if (cinfo->cntSinceActive > maxInactCutAge_ &&
        false == cinfo->neverDelete) {
      addToPool_(con);
      cinfo->cntSinceViol = 0;
      p_->markDelete(con->getConstraint());
      del_const = true;
      num_dels++;
    } else {
      cpyrel.push_back(con);
    }
  }  
  enCuts_.clear();
  enCuts_ = cpyrel;

  if (del_const == true){
    p_->delMarkedCons();
  }

  cpyrel.clear();
}


void CutMan1::separate(ConstSolutionPtr sol, bool *, UInt *)
{
  UInt n = p_->getNumVars();
  const double *x = new double[n];
  x = sol->getPrimal();
  CutPtr con;
  double viol = 0;
  CutList cpypool;
  int err;
  CutInfo *cinfo;
 
  for (CCIter it=pool_.begin(); it != pool_.end(); ++it)
  {
    con = *it;
    err = 0;
    con->eval(x, &err);
    cinfo = con->getInfo();
    if (viol < -absTol_) 
    {
      ++(cinfo->cntSinceViol);
      cpypool.push_back(con);
    } else if (viol > absTol_) {
      addToRel_(con, true);
    } else if (con->getInfo()->cntSinceActive <= maxDisCutAge_ ||
               con->getInfo()->neverDelete==true) {
      cpypool.push_back(con);
    }
  }
  pool_.clear();
  pool_ = cpypool;
  cpypool.clear();
}


void CutMan1::addCut(CutPtr c)
{
  addToRel_(c, true);
}


void CutMan1::addCuts(CutVectorIter cbeg, CutVectorIter cend)
{
  for (CutVectorIter it=cbeg; it!=cend; ++it) {
    addCut(*it);
  }
}


void CutMan1::addToRel_(CutPtr cut, bool )
{
  enCuts_.push_back(cut);
  ++(cut->getInfo()->numActive);
  cut->getInfo()->cntSinceActive = 0;
}


void CutMan1::addToPool_(CutPtr cut)
{
  pool_.push_back(cut);
}


void CutMan1::write(std::ostream &out) const
{
  out << me_ << "Nothing to write" << std::endl;
}


void CutMan1::writeStats(std::ostream &out) const
{
  out << me_ << "No stats availale" << std::endl;
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
