//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file YEqCGs.h
 * \brief Declare the class for storing auxiliary variables for nonlinear
 * functions that are expressed as computational graphs.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURYEQCG_H
#define MINOTAURYEQCG_H

#include "Types.h"
#include "OpCode.h"

namespace Minotaur {
class CGraph;
typedef boost::shared_ptr<CGraph> CGraphPtr;

class YEqCGs {
public:
  YEqCGs();
  VariablePtr findY(CGraphPtr cg);
  void insert(VariablePtr auxvar, CGraphPtr cg);

private:
  DoubleVector hash_;
  DoubleVector rand_;
  VarVector y_;
  std::vector<CGraphPtr> cg_;
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
