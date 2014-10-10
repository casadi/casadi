/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#ifndef CASADI_FUNCTOR_HPP
#define CASADI_FUNCTOR_HPP

#include "casadi_types.hpp"

#include "shared_object.hpp"

namespace casadi {

  /// Function pointer to a derivative generator function
  typedef Function (*DerivativeGeneratorCPtr)(Function& fcn, int nfwd, int nadj, void* user_data);

  /// Wrapper around functions
  typedef void (*CustomEvaluateCPtr)(CustomFunction &f, void* user_data);

  /// Wrapper around callback
  typedef int (*CallbackCPtr)(Function &f, void* user_data);

  /** \brief  Functor
      \author Joris Gillis
      \date 2013
  */
  class CASADI_CORE_EXPORT Functor : public SharedObject {
  };

  /** \brief Derivative Generator Functor
  *
  * In C++, supply a DerivativeGeneratorCPtr function pointer
  *
  * In python, supply a callable, annotated with derivativegenerator decorator
  * \code
  *
  *   @derivativegenerator
  *   def c(f, nadj, nadir):
  *     print f
  *
  *   ff.setOption("derivative_generator", c)
  * \endcode
  *
  */
  class CASADI_CORE_EXPORT DerivativeGenerator : public Functor {
    public:
      /// Default constructor
      DerivativeGenerator() { }

      /// Construct from C pointer
      DerivativeGenerator(DerivativeGeneratorCPtr ptr);

      /// Call
      Function operator() (Function& fcn, int nfwd, int nadj, void* user_data);
  };

  /** \brief CustomEvaluate
  *
  * In C++, supply a CustomEvaluateCPtr function pointer
  *
  * In python, supply a callable, annotated with pyevaluate decorator
  * \code
  *
  *   @pyevaluate
  *   def c(f, nfdir, nadir):
  *     print f
  *
  *   f = CustomFunction(c, ...)
  * \endcode
  *
  */
  class CASADI_CORE_EXPORT CustomEvaluate : public Functor {
    public:
      /// Default constructor
      CustomEvaluate() {}

      /// Construct from C pointer
      CustomEvaluate(CustomEvaluateCPtr ptr);

      /// Call
      void operator() (CustomFunction& fcn, void* user_data);
  };

  /** \brief Callback
  *
  * In C++, supply a CallbackCPtr function pointer
  * When the callback function returns a non-zero integer, the host is signalled of a problem.
  * E.g. an NlpSolver may halt iterations if the Callback is something else than 0
  *
  * In python, supply a callable, annotated with pycallback decorator
  * \code
  *
  *   @pycallback
  *   def c(f):
  *     print f
  *     return 0
  *
  *   solver.setOption("iteration_callback", c)
  * \endcode
  *
  */
  class CASADI_CORE_EXPORT Callback : public Functor {
    public:
      /// Default constructor
      Callback() {}

      /// Construct from C pointer
      Callback(CallbackCPtr ptr);

      /// Call
      int operator() (Function& fcn, void* user_data);
  };

} // namespace casadi


#endif // CASADI_FUNCTOR_HPP
