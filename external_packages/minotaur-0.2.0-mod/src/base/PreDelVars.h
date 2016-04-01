//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file PreDelVars.h
 * \brief Declare the PreDelVars class for saving and restoring variables in
 * presolve.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURPREDELVARS_H
#define MINOTAURPREDELVARS_H

#include "PreMod.h"
#include "Variable.h"

namespace Minotaur {

  class PreDelVars : public PreMod {
    public:
      /// Constructor
      PreDelVars();

      /// Destroy
      ~PreDelVars();

      /// Add a new variable to the list.
      void insert(VariablePtr v);

      /// Restore x.
      void postsolveGetX(const DoubleVector &x, DoubleVector *newx);

    private:
      /// A queue of variables deleted.
      VarQueue vars_;

  };

  typedef boost::shared_ptr<PreDelVars> PreDelVarsPtr;
  typedef boost::shared_ptr<const PreDelVars> ConstPreDelVarsPtr;
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
