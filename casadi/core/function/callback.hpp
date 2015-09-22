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

#ifndef CASADI_CALLBACK_HPP
#define CASADI_CALLBACK_HPP

#include "function_internal.hpp"
#include "../functor_internal.hpp"

namespace casadi {


class CASADI_EXPORT Callback2 {

  public:
    Callback2();

    //virtual Sparsity jacSparsity(int iind, ind oind) { }
    virtual std::vector<DMatrix> operator()(const std::vector<DMatrix>& arg);
    //virtual Function derForward(const std::string& name, int nfwd, const Dict& options) {}

    /** \brief Number of input arguments
    *
    * Specify the number of input arguments that a specific instance can handle.
    * The number must not be changed over the lifetime of the object
    *
    * Default implementation: 1
    *
    */
    virtual int nIn() { return 1;}
    /** \brief Number of output arguments
    *
    * Specify the number of output arguments that a specific instance can handle.
    * The number must not be changed over the lifetime of the object
    *
    * Default implementation: 1
    */
    virtual int nOut() { return 1;}
    /** \brief Specify input sparsity
    *
    * Specify the sparsity corresponding to a given input.
    * The sparsity must not be changed over the lifetime of the object
    *
    * Default implementation: dense using inputShape
    *
    */
    virtual Sparsity inputSparsity(int i) { return Sparsity::dense(inputShape(i)); }
    /** \brief Specify output sparsity
    *
    * Specify the sparsity corresponding to a given output.
    * The sparsity must not be changed over the lifetime of the object
    *
    * Default implementation: dense using outputShape
    *
    */
    virtual Sparsity outputSparsity(int i) { return Sparsity::dense(outputShape(i)); }
    /** \brief Specify input shape
    *
    * Specify the shape corresponding to a given input.
    * The shape must not be changed over the lifetime of the object
    *
    * Default implementation: scalar (1,1)
    *
    */
    virtual std::pair<int, int> inputShape(int i) { return std::pair<int, int>(1, 1); }
    /** \brief Specify output shape
    *
    * Specify the shape corresponding to a given output.
    * The shape must not be changed over the lifetime of the object
    *
    * Default implementation: scalar (1,1)
    *
    */
    virtual std::pair<int, int> outputShape(int i) { return std::pair<int, int>(1, 1); }
    /** \brief Specify the name of the object
    */
    virtual std::string name() { return "Custom callback"; }

    /** \brief Specify the options of the object
    */
    virtual Dict options() { return Dict(); }

    Function create();

    /** \brief  Destructor */
    virtual ~Callback2();


};

class CASADI_EXPORT DerivativeGenerator2 {

  public:
    DerivativeGenerator2();

    virtual Function operator()(Function& fcn, int ndir);

    /// Computes the derivative as if this derivative generator does not exist
    Function original(Function& fcn, int ndir, bool fwd);

    DerivativeGenerator create();

    /** \brief  Destructor */
    virtual ~DerivativeGenerator2();
};

#ifndef SWIG

class DerivativeGeneratorInternal2 : public DerivativeGeneratorInternal {
public:
  /** \brief  Create a function */
  explicit DerivativeGeneratorInternal2(DerivativeGenerator2 &callback);

  /** \brief  Destructor */
  virtual ~DerivativeGeneratorInternal2();

  /** \brief  Cloning */
  virtual DerivativeGeneratorInternal2* clone() const {
    return new DerivativeGeneratorInternal2(*this);
  }

  virtual Function call(Function& fcn, int ndir, void* user_data) { return callback_(fcn, ndir); }

  DerivativeGenerator2& callback_;
};


class CASADI_EXPORT CallbackFunctionInternal : public FunctionInternal {
  friend class CallbackFunction;
  public:


    /** \brief  Create a function */
    explicit CallbackFunctionInternal(Callback2 &callback);

    /** \brief  Destructor */
    virtual ~CallbackFunctionInternal();

    /** \brief  Cloning */
    virtual CallbackFunctionInternal* clone() const { return new CallbackFunctionInternal(*this);}

    virtual void evalD(const double** arg,
                               double** res, int* iw, double* w);

    /** \brief  Initialize */
    virtual void init();

    Callback2& callback_;

    double eps_;

};
#endif


} // namespace casadi

#endif // CASADI_CALLBACK_HPP
