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

#ifndef FUNCTOR_INTERNAL_HPP
#define FUNCTOR_INTERNAL_HPP

#include "functor.hpp"
/// \cond INTERNAL

namespace CasADi{

  
  /** \brief Internal class for Callback
      \author Joris Gillis 
      \date 2013
  */
  class FunctorInternal : public SharedObjectNode {
    friend class Functor;
    
  };
  
  template<typename P>
  class FunctorCInternal {
    public:
      FunctorCInternal(P ptr) : ptr_(ptr) {};
    protected:
      P ptr_;
  };
  
  class DerivativeGeneratorInternal : public FunctorInternal {
    friend class DerivativeGenerator;
    virtual Function call(Function& fcn, int nfwd, int nadj, void* user_data)=0;
  };

  class DerivativeGeneratorCInternal : public DerivativeGeneratorInternal, FunctorCInternal<DerivativeGeneratorCPtr> {
    friend class DerivativeGenerator;
    
    DerivativeGeneratorCInternal(DerivativeGeneratorCPtr ptr);
    virtual Function call(Function& fcn, int nfwd, int nadj, void* user_data);
    virtual DerivativeGeneratorCInternal* clone() const;
  };

  class CustomEvaluateInternal : public FunctorInternal {
    friend class CustomEvaluate;
    virtual void call(CustomFunction& fcn, void* user_data)=0;
  };

  class CustomEvaluateCInternal : public CustomEvaluateInternal, FunctorCInternal<CustomEvaluateCPtr> {
    friend class CustomEvaluate;
    
    CustomEvaluateCInternal(CustomEvaluateCPtr ptr);
    virtual void call(CustomFunction& fcn, void* user_data);
    virtual CustomEvaluateCInternal* clone() const;
  };
  
  class CallbackInternal : public FunctorInternal {
    friend class Callback;
    virtual int call(Function& fcn, void* user_data)=0;
  };

  class CallbackCInternal : public CallbackInternal, FunctorCInternal<CallbackCPtr> {
    friend class Callback;
    
    CallbackCInternal(CallbackCPtr ptr);
    virtual int call(Function& fcn, void* user_data);
    virtual CallbackCInternal* clone() const;
  };

} // namespace CasADi
/// \endcond

#endif // FUNCTOR_INTERNAL_HPP
