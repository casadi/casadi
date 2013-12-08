/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef FUNCTOR_HPP
#define FUNCTOR_HPP

#include "casadi_types.hpp"

#include "shared_object.hpp"

namespace CasADi{

  /// Wrapper around functions
  typedef void (*CustomEvaluateCPtr)(CustomFunction &f, void* user_data);

  /// Wrapper around callback
  typedef int (*CallbackCPtr)(FX &f, void* user_data);
  
  class Functor;
  
  /** \brief Internal class for Functor
      \author Joris Gillis 
      \date 2013
  */
  class Functor : public SharedObject {
    //Callback();
    
  };    
    
  /** \brief CustomEvaluate
  *
  * In C++, supply a CustomEvaluateCPtr function pointer
  *
  * In python, supply a callable, annotated with pyevaluate decorator
  * \code
  *  
  *   @pyevaluate
  *   def c(f,nfdir,nadir):
  *     print f
  *
  *   f = CustomFunction(c,...)
  * \endcode
  *
  */
  class CustomEvaluate : public Functor {
    public:
      /// Default constructor
      CustomEvaluate() {};
      /// Construct from C pointer
      CustomEvaluate(CustomEvaluateCPtr ptr);
      /// Call
      virtual void operator() (CustomFunction& fcn, void* user_data);
  };

  /** \brief Callback
  * 
  * In C++, supply a CallbackCPtr function pointer
  * When the callback function returns a non-zero integer, the host is signalled of a problem.
  * E.g. an NLPSolver may halt iterations if the Callback is something else than 0
  *
  * In python, supply a callable, annotated with pycallback decorator
  * \code
  *  
  *   @pycallback
  *   def c(f):
  *     print f
  *     return 0
  *
  *   solver.setOption("iteration_callback",c)
  * \endcode
  *
  */
  class Callback : public Functor {
    public:
      /// Default constructor
      Callback() {};
      /// Construct from C pointer
      Callback(CallbackCPtr ptr);
      /// Call
      virtual int operator() (FX& fcn, void* user_data);
  };

} // namespace CasADi


#endif // FUNCTOR_HPP
