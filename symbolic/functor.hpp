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

  /// Function pointer to a Jacobian generator function
  typedef FX (*JacobianGeneratorCPtr)(FX& fcn, int iind, int oind, void* user_data);
  
  /// Function pointer to a sparsity generator function
  typedef CRSSparsity (*SparsityGeneratorCPtr)(FX& fcn, int iind, int oind, void* user_data);
  
  /// Wrapper around functions
  typedef void (*CustomEvaluateCPtr)(CustomFunction &f, int nfdir, int nadir, void* user_data);

  /// Wrapper around callback
  typedef int (*CallbackCPtr)(const FX &f, void* user_data);
  
  class Functor;
  
  /** \brief Internal class for Functor
      \author Joris Gillis 
      \date 2013
  */
  class Functor : public SharedObject {
    //Callback();
    
  };
  
  
  /** \brief Sparsity Generator Functor */
  class SparsityGenerator : public Functor {
    public:
      /// Call
      virtual CRSSparsity operator() (FX& fcn, int iind, int oind, void* user_data);
  };

  /** \brief Sparsity Generator Functor C Implementation */
  class SparsityGeneratorC : public SparsityGenerator {
    /// Default constructor
    //SparsityGeneratorC();
    public:
      /// Construct from C pointer
      SparsityGeneratorC(SparsityGeneratorCPtr ptr);
  };
  

  /** \brief Jacobian Generator Functor */
  class JacobianGenerator : public Functor {
    public:
      /// Call
      virtual FX operator() (FX& fcn, int iind, int oind, void* user_data);
  };

  /** \brief Sparsity Generator Functor C Implementation */
  class JacobianGeneratorC : public JacobianGenerator {
    public:
      /// Construct from C pointer
      JacobianGeneratorC(JacobianGeneratorCPtr ptr);
  };
  
  
    
  /** \brief CustomEvaluate */
  class CustomEvaluate : public Functor {
    public:
      /// Call
      virtual void operator() (CustomFunction& fcn, int nfdir, int nadir, void* user_data);
  };

  /** \brief CustomEvaluate C Implementation */
  class CustomEvaluateC : public SparsityGenerator {
    public:
      /// Construct from C pointer
      CustomEvaluateC(CustomEvaluateCPtr ptr);
  };

  /** \brief Callback */
  class Callback : public Functor {
    public:
      /// Call
      virtual int operator() (const FX& fcn, void* user_data);
  };

  /** \brief Callback C Implementation */
  class CallbackC : public Callback {
    public:
      /// Construct from C pointer
      CallbackC(CallbackCPtr ptr);
  };

} // namespace CasADi


#endif // FUNCTOR_HPP
