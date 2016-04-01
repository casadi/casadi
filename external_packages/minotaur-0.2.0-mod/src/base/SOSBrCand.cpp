//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file SOSBrCand.cpp
 * \brief Define methods for SOS branching candidates (both Type 1 and 2).
 * \author Ashutosh Mahajan, IIT Bombay
 */

#include <cassert>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "MinotaurConfig.h"
#include "BrCand.h"
#include "Logger.h"
#include "SOS.h"
#include "SOSBrCand.h"
#include "Variable.h"

using namespace Minotaur;

SOSBrCand::SOSBrCand() 
: lsum_(0.0),
  rsum_(0.0),
  sos_(0)
{
  pCostIndex_ = -1;
  score_      = 0.0;
}


SOSBrCand::SOSBrCand(const SOS* sos, VarVector &left, VarVector &right,
                     double lsum, double rsum)
{
  lvars_ = left;
  rvars_ = right;
  lsum_ = lsum;
  rsum_ = rsum;
  sos_ = sos;
  pCostIndex_ = -1*sos->getId()-1;
  score_      = (double) sos->getPriority();
}


SOSBrCand::~SOSBrCand()
{
  lvars_.clear();
  rvars_.clear();
}


double SOSBrCand::getDDist()
{
  return lsum_;
}


double SOSBrCand::getLSum() const
{
  return rsum_;
}


std::string SOSBrCand::getName() const
{
  return sos_->getName();
}


double SOSBrCand::getRSum() const
{
  return rsum_;
}


double SOSBrCand::getUDist()
{
  return rsum_;
}


VariableConstIterator SOSBrCand::lVarsBegin() const
{
  return lvars_.begin();
}


VariableConstIterator SOSBrCand::lVarsEnd() const
{
  return lvars_.end();
}


VariableConstIterator SOSBrCand::rVarsBegin() const
{
  return rvars_.begin();
}


VariableConstIterator SOSBrCand::rVarsEnd() const
{
  return rvars_.end();
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
