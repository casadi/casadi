//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file PreMod.h
 * \brief Declare the PreMod class for saving changes in presolve.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURPREMOD_H
#define MINOTAURPREMOD_H

#include "Types.h"

namespace Minotaur {

  /** 
   * PreMod class is an abstract base class. It has a method to transform a
   * solution of the presolved problem in to a solution of a problem in which
   * some of the modifications have been undone.
   */
  class PreMod {
    public:
      /// Constructor.
      PreMod() {};

      /// Destroy.
      virtual ~PreMod() {};

      /// Restore x.
      virtual void postsolveGetX(const DoubleVector &x, DoubleVector *newx) = 0;

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
