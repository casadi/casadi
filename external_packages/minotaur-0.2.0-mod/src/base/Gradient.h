// 
//      MINOTAUR -- It's only 1/2 bull
// 
//      (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file Gradient.h
 * \brief Get information about gradient of a function.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURGRADIENT_H
#define MINOTAURGRADIENT_H

#include "linalg/Vector.h"
#include "Variable.h"

namespace Minotaur {

  class Gradient {
    public:
      Gradient() {};

      // only the following coordinates have non-zero components in the
      // gradient.
      //virtual VariableIterator varsBegin() { return vars_.begin(); }
      //virtual VariableIterator varsEnd() { return vars_.end(); }

      // initialize
      virtual void initialize() {};

      // get no. of nonzeros in the gradient vector
      virtual UInt getNumNz() const {};

      // evaluate gradient at x
      virtual LinearAlgebra::ConstDoubleVectorPtr 
        eval(const std::vector<double> &x) const { assert(!"implement me!"); }

      // evaluate gradient at x
      virtual LinearAlgebra::ConstDoubleVectorPtr 
        eval (LinearAlgebra::ConstDoubleVectorPtr x) { assert(!"implement me!"); }
      //virtual bool isZero() const { return true; }

    protected:
      // LinearGradientPtr lgrad_;
      // QuadraticGradientPtr qgrad_;
  };

  typedef boost::shared_ptr<Gradient> GradientPtr;
  typedef boost::shared_ptr<const Gradient> ConstGradientPtr;  
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
