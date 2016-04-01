//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file SOS.cpp
 * \brief SOS class stores data regarding SOS1 and SOS2 constraints.
 * \author Ashutosh Mahajan, IIT Bombay
 */

#include <cassert>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string.h> // for memcpy

#include "MinotaurConfig.h"
#include "Operations.h"
#include "SOS.h"

using namespace Minotaur;

SOS::SOS() 
: n_(0),
  priority_(0),
  type_(SOS1),
  weights_(0)
{
}


SOS::SOS(int n, SOSType type, const double *weights, const VarVector &vars,
         int priority, int id, std::string name) 
: id_(id),
  n_(n),
  priority_(priority),
  type_(type),
  name_(name)
{
  weights_ = new double[n];
  memcpy(weights_, weights, sizeof(double)*n);
  vars_ = vars;
  sort(vars_, weights_, -1);
}


SOS::~SOS()
{
  vars_.clear();
  delete [] weights_;
}


int SOS::getId() const
{
  return id_;
}


int SOS::getNz()
{
  return n_;
}


int SOS::getPriority() const
{
  return priority_;
}


SOSType SOS::getType()
{
  return type_;
}


std::string SOS::getName() const
{
  return name_;
}



const double* SOS::getWeights()
{
  return weights_;
}


VariableConstIterator SOS::varsBegin() const
{
  return vars_.begin();
}


VariableConstIterator SOS::varsEnd() const
{
  return vars_.end();
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
