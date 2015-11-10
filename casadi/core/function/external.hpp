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


#ifndef CASADI_EXTERNAL_HPP
#define CASADI_EXTERNAL_HPP

#include "function_internal.hpp"

#ifdef WITH_DL
#ifdef _WIN32 // also for 64-bit
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0502
#endif
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
    handle_t handle_;

  public:
    // Default constructor
    LibInfo() : handle_(0) {}

    // Constructor
    explicit LibInfo(const std::string& bin_name);

    // Clear memory
    void clear();

    // Automatic type conversion
    operator const std::string&() const { return bin_name_;}

    // Get function pointer
    template<typename FcnPtr>
    void get(FcnPtr& fcnPtr, const std::string& sym);
  };

  /** \brief Library that has been just-in-time compiled */
  template<>
  class LibInfo<Compiler> {
  private:
    Compiler compiler_;

  public:
    // Default constructor
    LibInfo() {}

    // Constructor
    explicit LibInfo(const Compiler& compiler);

    // Clear memory
    void clear() {}

    // Automatic type conversion
    operator const Compiler&() const { return compiler_;}

    // Get function pointer
    template<typename FcnPtr>
    void get(FcnPtr& fcnPtr, const std::string& sym);
  };

  class CASADI_EXPORT External : public FunctionInternal {
  public:
    /** \brief Creator function, dynamically linked library */
    static External*
      create(const std::string& bin_name, const std::string& name);

    /** \brief Creator function, just-in-time compiled library */
    static External*
      create(const Compiler& compiler, const std::string& name);

    /** \brief Constructor */
    External(const std::string& name);

    /** \brief Destructor */
    virtual ~External() = 0;

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int ind) const;
    virtual Sparsity get_sparsity_out(int ind) const;
    virtual Sparsity get_sparsity(int ind) const = 0;
    /// @}

  private:
    /** \brief Creator function, use this for creating instances of the class */
    template<typename LibType>
      static External*
      createGeneric(const LibType& bin_name, const std::string& name);
  };

  template<typename LibType>
  class CASADI_EXPORT CommonExternal : public External {
  protected:

    /** \brief Information about the library */
    LibInfo<LibType> li_;

    /** \brief Function */
    void *mem_;

    /** \brief Number of inputs and outputs */
    int n_in_, n_out_;

    /// @{
    /** \brief Retreive sparsities */
    sparsity_t sparsity_;
    virtual Sparsity get_sparsity(int ind) const;
    /// @}

    /** \brief Free memory */
    freemem_t freemem_;

    /** \brief  constructor is protected */
    CommonExternal(const std::string& name, const LibInfo<LibType>& li);
  public:
    /** \brief  Destructor */
    virtual ~CommonExternal() = 0;

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() const { return n_in_;}
    virtual size_t get_n_out() const { return n_out_;}
    ///@}

    ///@{
    /** \brief Forward mode derivatives */
    virtual Function get_forward(const std::string& name, int nfwd, Dict& opts);
    virtual int get_n_forward() const;
    ///@}

    ///@{
    /** \brief Reverse mode derivatives */
    virtual Function get_reverse(const std::string& name, int nadj, Dict& opts);
    virtual int get_n_reverse() const;
    ///@}

    ///@{
    /** \brief Full Jacobian */
    virtual bool hasFullJacobian() const;
    virtual Function getFullJacobian(const std::string& name, const Dict& opts);
    ///@}
  };

  template<typename LibType>
  class CASADI_EXPORT SimplifiedExternal : public CommonExternal<LibType> {
    friend class External;
  private:
    /** \brief Constructor is called from factory */
    SimplifiedExternal(const std::string& name, const LibInfo<LibType>& li);
  public:
    using CommonExternal<LibType>::li_;

    /** \brief  Destructor */
    virtual ~SimplifiedExternal() {}

    /** \brief  Evaluate numerically, work vectors given */
    virtual void eval(const double** arg, double** res, int* iw, double* w, void* mem);

    /** \brief Add a dependent function */
    virtual void addDependency(CodeGenerator& g) const;

    /** \brief Generate a call to a function (simplified signature) */
    virtual std::string simple_call(const CodeGenerator& g,
                                    const std::string& arg, const std::string& res) const;

    /** \brief Use simplified signature */
    virtual bool simplifiedCall() const { return true;}
  protected:
    /** \brief  Function pointers */
    simple_t eval_;
  };

  template<typename LibType>
  class CASADI_EXPORT GenericExternal : public CommonExternal<LibType> {
    friend class External;
  private:
    /** \brief Constructor is called from factory */
    GenericExternal(const std::string& name, const LibInfo<LibType>& li);
  public:
    using CommonExternal<LibType>::li_;

    /** \brief  Destructor */
    virtual ~GenericExternal() {}

    /** \brief  Evaluate numerically, work vectors given */
    virtual void eval(const double** arg, double** res, int* iw, double* w, void* mem);

    /** \brief Add a dependent function */
    virtual void addDependency(CodeGenerator& g) const;

    /** \brief Generate a call to a function (generic signature) */
    virtual std::string generic_call(const CodeGenerator& g, const std::string& mem,
                                     const std::string& arg, const std::string& res,
                                     const std::string& iw, const std::string& w) const;

    /** \brief All inputs and outputs are scalar (default if sparsity not defined) */
    static int scalarSparsity(void* mem, int i, int *n_row, int *n_col,
                              const int **colind, const int **row);

    /** \brief  Function pointers */
    eval_t eval_;
  };


} // namespace casadi
/// \endcond

#endif // CASADI_EXTERNAL_HPP
