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


#ifndef CASADI_EXTERNAL_FUNCTION_INTERNAL_HPP
#define CASADI_EXTERNAL_FUNCTION_INTERNAL_HPP

#include "external_function.hpp"
#include "function_internal.hpp"

#ifdef WITH_DL
#ifdef _WIN32 // also for 64-bit
#define NOMINMAX
#include <windows.h>
#else // _WIN32
#include <dlfcn.h>
#endif // _WIN32
#endif // WITH_DL

/// \cond INTERNAL

namespace casadi {

  /** \brief Structure with information about the library */
  template<typename LibType>
  class LibInfo {};

  /** \brief Library given as a dynamically linked library */
  template<>
  class LibInfo<std::string> {
  private:
#if defined(WITH_DL) && defined(_WIN32) // also for 64-bit
    typedef HINSTANCE handle_t;
#else
    typedef void* handle_t;
#endif

    std::string bin_name_;
    std::string name_;
    handle_t handle_;

  public:
    int n_in, n_out, sz_arg, sz_res;

    // Default constructor
    LibInfo() : handle_(0) {}

    // Constructor
    explicit LibInfo(const std::string& bin_name, const std::string& name);

    // Clear memory
    void clear();

    // Automatic type conversion
    operator const std::string&() const { return bin_name_;}

    // Function name
    const std::string& name() const { return name_;}

    // Get function pointer
    template<typename FcnPtr>
    void get(FcnPtr& fcnPtr, const std::string& sym);
  };

  /** \brief Library that has been just-in-time compiled */
  template<>
  class LibInfo<Compiler> {
  private:
    Compiler compiler_;
    std::string name_;

  public:
    int n_in, n_out, sz_arg, sz_res;

    // Default constructor
    LibInfo() {}

    // Constructor
    explicit LibInfo(const Compiler& compiler, const std::string& name);

    // Clear memory
    void clear() {}

    // Automatic type conversion
    operator const Compiler&() const { return compiler_;}

    // Function name
    const std::string& name() const { return name_;}

    // Get function pointer
    template<typename FcnPtr>
    void get(FcnPtr& fcnPtr, const std::string& sym);
  };

  class CASADI_EXPORT ExternalFunctionInternal : public FunctionInternal {
  public:
    /** \brief Creator function, dynamically linked library */
    static ExternalFunctionInternal*
      create(const std::string& bin_name, const std::string& name);

    /** \brief Creator function, just-in-time compiled library */
    static ExternalFunctionInternal*
      create(const Compiler& compiler, const std::string& name);

    /** \brief  clone function */
    virtual ExternalFunctionInternal* clone() const;

    /** \brief  Destructor */
    virtual ~ExternalFunctionInternal() = 0;

  private:
    /** \brief Creator function, use this for creating instances of the class */
    template<typename LibType>
      static ExternalFunctionInternal*
      createGeneric(const LibType& bin_name, const std::string& name);
  };

  template<typename LibType>
  class CASADI_EXPORT CommonExternal : public ExternalFunctionInternal {
  protected:

    /** \brief Information about the library */
    LibInfo<LibType> li_;

    /** \brief  constructor is protected */
    CommonExternal(const LibInfo<LibType>& li);
  public:
    /** \brief  Destructor */
    virtual ~CommonExternal() = 0;

    ///@{
    /** \brief Forward mode derivatives */
    virtual Function getDerForward(const std::string& name, int nfwd, Dict& opts);
    virtual int numDerForward() const;
    ///@}

    ///@{
    /** \brief Reverse mode derivatives */
    virtual Function getDerReverse(const std::string& name, int nadj, Dict& opts);
    virtual int numDerReverse() const;
    ///@}

    ///@{
    /** \brief Full Jacobian */
    virtual bool hasFullJacobian() const;
    virtual Function getFullJacobian(const std::string& name, const Dict& opts);
    ///@}
  };

  template<typename LibType>
  class CASADI_EXPORT SimplifiedExternal : public CommonExternal<LibType> {
    friend class ExternalFunctionInternal;
  private:
    /** \brief Constructor is called from factory */
    SimplifiedExternal(const LibInfo<LibType>& li);
  public:
    using CommonExternal<LibType>::li_;

    /** \brief  Destructor */
    virtual ~SimplifiedExternal() {}

    /** \brief  Evaluate numerically, work vectors given */
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /** \brief Add a dependent function */
    virtual void addDependency(CodeGenerator& g) const;

    /** \brief Generate a call to a function (simplified signature) */
    virtual std::string generateCall(const CodeGenerator& g,
                                     const std::string& arg, const std::string& res) const;

    /** \brief Use simplified signature */
    virtual bool simplifiedCall() const { return true;}
  protected:
    /** \brief  Function pointers */
    simplifiedPtr eval_;
  };

  template<typename LibType>
  class CASADI_EXPORT GenericExternal : public CommonExternal<LibType> {
    friend class ExternalFunctionInternal;
  private:
    /** \brief Constructor is called from factory */
    GenericExternal(const LibInfo<LibType>& li);
  public:
    using CommonExternal<LibType>::li_;

    /** \brief  Destructor */
    virtual ~GenericExternal() {}

    /** \brief  Evaluate numerically, work vectors given */
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /** \brief Add a dependent function */
    virtual void addDependency(CodeGenerator& g) const;

    /** \brief Generate a call to a function (generic signature) */
    virtual std::string generateCall(const CodeGenerator& g,
                                     const std::string& arg, const std::string& res,
                                     const std::string& iw, const std::string& w) const;

    /** \brief All inputs and outputs are scalar (default if sparsity not defined) */
    static int scalarSparsity(int i, int *n_row, int *n_col,
                              const int **colind, const int **row);

    /** \brief  Function pointers */
    evalPtr eval_;
  };


} // namespace casadi
/// \endcond

#endif // CASADI_EXTERNAL_FUNCTION_INTERNAL_HPP
