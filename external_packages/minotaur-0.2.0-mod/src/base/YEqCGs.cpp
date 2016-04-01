//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file YEqCGs.cpp
 * \brief Define class for storing auxiliary variables for nonlinear functions
 * that are expressed as computational graphs.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>

#include "MinotaurConfig.h"
#include "CGraph.h"
#include "CNode.h"
#include "OpCode.h"
#include "Variable.h"
#include "YEqCGs.h"

// #define SPEW 1

using namespace Minotaur;


YEqCGs::YEqCGs()
{
  UInt numops = 35;
  for (UInt i=0; i<numops; ++i) {
    rand_.push_back((double) rand()/(RAND_MAX));
  }
}


double YEqCGs::evalHash_(const CNode* node, UInt rank)
{
  OpCode op = node->getOp();
  double hash = rand_[op]*rank;
  if (OpVar==op) {
    hash *= node->getV()->getId();
  } else if (OpInt==op || OpNum==op) {
    hash *= node->getVal();
  } else if (1==node->numChild()) {
    hash += evalHash_(node->getL(), rank+1);
  } else if (2==node->numChild()) {
    hash += evalHash_(node->getL(), rank+1);
    hash += evalHash_(node->getR(), rank+2);
  } else if (2<node->numChild()) {
    UInt i= 0;
    CNode** c1 = node->getListL();
    CNode** c2 = node->getListR();
    while (c1<c2) {
      hash += evalHash_(*c1, rank+i);
      ++i;
      ++c1;
    }
  }
  return hash;
}


VariablePtr YEqCGs::findY(CGraphPtr cg)
{
  double hash = evalHash_(cg->getOut(), 1);
  for (UInt i=0; i<y_.size(); ++i) {
    if (fabs(hash-hash_[i])<1e-10 &&
        cg->isIdenticalTo(cg_[i])) {
      return y_[i];
    }
  }
  return VariablePtr();
}


void YEqCGs::insert(VariablePtr auxvar, CGraphPtr cg)
{
  assert(auxvar && cg);
  hash_.push_back(evalHash_(cg->getOut(), 1));
  y_.push_back(auxvar);
  cg_.push_back(cg);
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
