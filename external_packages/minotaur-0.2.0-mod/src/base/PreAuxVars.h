//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file PreAuxVars.h
 * \brief Declare the PreAuxVars class for saving additional variables
 * variables added during presolve.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURPREAUXVARS_H
#define MINOTAURPREAUXVARS_H

#include "PreMod.h"
#include "Variable.h"

namespace Minotaur {

class PreAuxVars : public PreMod {
public:
  /// Constructor.
  PreAuxVars();

  /// Destroy.
  ~PreAuxVars();

  /// New variable that was created.
  void insert(VariablePtr v);

  /// Remove aux-vars from the solution x.
  void postsolveGetX(const DoubleVector &x, DoubleVector *newx);

  /// Return the number of additions.
  UInt getSize();

private:
  std::deque<VariablePtr> vars_;

};

typedef boost::shared_ptr<PreAuxVars> PreAuxVarsPtr;
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
