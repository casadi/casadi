//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file PreSubstVars.h
 * \brief Declare the PreSubstVars class for saving and restoring 
 * variables substituted during presolve.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURPRESUBSTVARS_H
#define MINOTAURPRESUBSTVARS_H

#include "PreMod.h"
#include "Variable.h"

namespace Minotaur {

struct PreSubstVarData {
  VariablePtr vout; /// Number of nlps solved.
  UInt vinInd;      /// Number of nlps feasible.
  double rat;       /// Number of nlps infeasible.
}; 

class PreSubstVars : public PreMod {
public:
  /// Constructor.
  PreSubstVars();

  /// Destroy.
  ~PreSubstVars();

  /// Substitute variable 'vin' by variable 'vout'.
  void insert(VariablePtr vout, VariablePtr vin, double rat = 1.0);

  /// Restore x.
  void postsolveGetX(const DoubleVector &x, DoubleVector *newx);

  /// Return the number of substitutions.
  UInt getSize();

private:
  std::deque<PreSubstVarData*> vars_;

};

typedef boost::shared_ptr<PreSubstVars> PreSubstVarsPtr;
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
