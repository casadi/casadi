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


#ifndef CASADI_EXTERNAL_IMPL_HPP
#define CASADI_EXTERNAL_IMPL_HPP

#include "external.hpp"
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
  class CASADI_EXPORT External : public FunctionInternal {
  protected:
    /** \brief Information about the library */
    Library li_;

    /** \brief Initialize (called before every other function) */
    init_t init_;

    /** \brief Number of inputs and outputs */
    getint_t n_in_, n_out_;

    /** \brief Work vector sizes */
    work_t work_;

    /** \brief Allocate memory */
    checkout_t checkout_;

    /** \brief Free memory */
    release_t release_;

    ///@{
    /** \brief Data vectors */
    std::vector<int> int_data_;
    std::vector<double> real_data_;
    std::string string_data_;
    ///@}
  public:

    /** \brief Constructor */
    External(const std::string& name, const Library& li);

    /** \brief Destructor */
    virtual ~External() = 0;

    /** \brief Get type name */
    virtual std::string type_name() const { return "external";}

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /// Initialize
    virtual void init(const Dict& opts);

    /** \brief Add a dependent function */
    virtual void addDependency(CodeGenerator& g) const;

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() const;
    virtual size_t get_n_out() const;
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

    /** \brief Create memory block */
    virtual void* alloc_memory() const {
      if (checkout_) {
        int* m = new int;
        *m = checkout_();
        return m;
      } else if (!derivative_of_.is_null()) {
        return derivative_of_->alloc_memory();
      } else {
        return 0;
      }
    }

    /** \brief Free memory block */
    virtual void free_memory(void *mem) const {
      if (release_) {
        int* m = static_cast<int*>(mem);
        release_(*m);
      } else if (!derivative_of_.is_null()) {
        derivative_of_->free_memory(mem);
      }
    }
  };

  class CASADI_EXPORT SimplifiedExternal : public External {
    friend class External;
  public:
    /** \brief Constructor */
    SimplifiedExternal(const std::string& name, const Library& li);

    /** \brief  Destructor */
    virtual ~SimplifiedExternal() { this->clear_memory();}

    /// Initialize
    virtual void init(const Dict& opts);

    /** \brief  Evaluate numerically */
    virtual void simple(const double* arg, double* res);

    /** \brief Generate a call to a function (simplified signature) */
    virtual std::string simple_call(const CodeGenerator& g,
                                    const std::string& arg, const std::string& res) const;

    /** \brief Use simplified signature */
    virtual bool simplifiedCall() const { return true;}

    /// @{
    /** \brief Retreive sparsities */
    virtual Sparsity get_sparsity_in(int ind) const { return Sparsity::scalar();}
    virtual Sparsity get_sparsity_out(int ind) const { return Sparsity::scalar();}
    /// @}
  protected:
    /** \brief  Function pointers */
    simple_t eval_;
  };

  class CASADI_EXPORT GenericExternal : public External {
    friend class External;
  public:
    /** \brief Constructor */
    GenericExternal(const std::string& name, const Library& li);

    /** \brief  Destructor */
    virtual ~GenericExternal() { this->clear_memory();}

    /// Initialize
    virtual void init(const Dict& opts);

    /// @{
    /** \brief Retreive sparsities */
    virtual Sparsity get_sparsity_in(int ind) const;
    virtual Sparsity get_sparsity_out(int ind) const;
    /// @}

    /** \brief  Evaluate numerically, work vectors given */
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /** \brief Generate a call to a function (generic signature) */
    virtual std::string generic_call(const CodeGenerator& g, const std::string& arg,
                                     const std::string& res, const std::string& iw,
                                     const std::string& w, const std::string& mem) const;
    // Sparsities
    sparsity_t sparsity_in_, sparsity_out_;

    /** \brief  Function pointers */
    eval_t eval_;
  };


} // namespace casadi
/// \endcond

#endif // CASADI_EXTERNAL_IMPL_HPP
