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

namespace CasADi{

  
  /** \brief Internal class for Callback
      \author Joris Gillis 
      \date 2013
  */
  class FunctorInternal : public SharedObjectNode {
    friend class Functor;
    
  };
  
  class SparsityGeneratorInternal : public FunctorInternal {
    friend class SparsityGenerator;
    virtual CRSSparsity call(FX& fcn, int iind, int oind, void* user_data)=0;
  };
  
  template<typename P>
  class FunctorCInternal {
    public:
      FunctorCInternal(P ptr) : ptr_(ptr) {};
    protected:
      P ptr_;
  };

  class SparsityGeneratorCInternal : public SparsityGeneratorInternal, FunctorCInternal<SparsityGeneratorCPtr> {
    friend class SparsityGeneratorC;
    
    SparsityGeneratorCInternal(SparsityGeneratorCPtr ptr);
    virtual CRSSparsity call(FX& fcn, int iind, int oind, void* user_data);
    virtual SparsityGeneratorCInternal* clone() const;

  };

  class JacobianGeneratorInternal : public FunctorInternal {
    friend class JacobianGenerator;
    virtual FX call(FX& fcn, int iind, int oind, void* user_data)=0;
  };

  class JacobianGeneratorCInternal : public JacobianGeneratorInternal, FunctorCInternal<JacobianGeneratorCPtr> {
    friend class JacobianGeneratorC;
    
    JacobianGeneratorCInternal(JacobianGeneratorCPtr ptr);
    virtual FX call(FX& fcn, int iind, int oind, void* user_data);
    virtual JacobianGeneratorCInternal* clone() const;
  };
  
  class CustomEvaluateInternal : public FunctorInternal {
    friend class CustomEvaluate;
    virtual void call(CustomFunction& fcn, int nfdir, int nadir, void* user_data)=0;
  };

  class CustomEvaluateCInternal : public CustomEvaluateInternal, FunctorCInternal<CustomEvaluateCPtr> {
    friend class CustomEvaluateC;
    
    CustomEvaluateCInternal(CustomEvaluateCPtr ptr);
    virtual void call(CustomFunction& fcn, int nfdir, int nadir, void* user_data);
    virtual CustomEvaluateCInternal* clone() const;
  };
  
  class CallbackInternal : public FunctorInternal {
    friend class Callback;
    virtual int call(FX& fcn, void* user_data)=0;
  };

  class CallbackCInternal : public CallbackInternal, FunctorCInternal<CallbackCPtr> {
    friend class CallbackC;
    
    CallbackCInternal(CallbackCPtr ptr);
    virtual int call(FX& fcn, void* user_data);
    virtual CallbackCInternal* clone() const;
  };

} // namespace CasADi


#endif // FUNCTOR_INTERNAL_HPP
