//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file YEqUCGs.h
 * \brief Declare the class for storing auxiliary variables for univariate
 * nonlinear functions that are expressed as computational graphs.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURYEQUCGS_H
#define MINOTAURYEQUCGS_H

#include "Types.h"
#include "OpCode.h"

namespace Minotaur {
class CGraph;
class CNode;
typedef boost::shared_ptr<CGraph> CGraphPtr;

class YEqUCGs {
public:
  YEqUCGs();
  ~YEqUCGs();
  VariablePtr findY(CGraphPtr cg);
  void insert(VariablePtr auxvar, CGraphPtr cg);

private:
  std::vector<CGraphPtr> cg_;
  DoubleVector hash_;
  std::vector<OpCode> op_;
  DoubleVector rand_;
  VarVector x_;
  VarVector y_;
  double evalHash_(const CNode* node, UInt rank);
};
}
#endif

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
