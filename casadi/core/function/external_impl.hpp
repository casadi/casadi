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
    Importer li_;

    /** \brief Increase/decrease reference counter */
    signal_t incref_, decref_;

    /** \brief Number of inputs and outputs */
    getint_t n_in_, n_out_;

    /** \brief Names of inputs and outputs */
    name_t name_in_, name_out_;

    /** \brief Work vector sizes */
    work_t work_;

    ///@{
    /** \brief Data vectors */
    std::vector<int> int_data_;
    std::vector<double> real_data_;
    std::string string_data_;
    ///@}
  public:

    /** \brief Constructor */
    External(const std::string& name, const Importer& li);

    /** \brief Destructor */
    virtual ~External() = 0;

    // Factory
    virtual Function factory(const std::string& name,
                             const std::vector<std::string>& s_in,
                             const std::vector<std::string>& s_out,
                             const Function::AuxOut& aux,
                             const Dict& opts) const;

    /** \brief Get type name */
    virtual std::string type_name() const { return "external";}

    /// Initialize
    virtual void init(const Dict& opts);

    /** \brief Add a dependent function */
    virtual void addDependency(CodeGenerator& g) const;

    /** \brief Generate code the function */
    virtual void generateFunction(CodeGenerator& g, const std::string& fname,
                                  bool decl_static) const;

    /** \brief Get name in codegen */
    virtual std::string codegen_name(const CodeGenerator& g) const;

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in();
    virtual size_t get_n_out();
    ///@}

    ///@{
    /** \brief Names of function input and outputs */
    virtual std::string get_name_in(int i);
    virtual std::string get_name_out(int i);
    /// @}

    ///@{
    /** \brief Forward mode derivatives */
    virtual Function get_forward(const std::string& name, int nfwd,
                                 const std::vector<std::string>& i_names,
                                 const std::vector<std::string>& o_names,
                                 const Dict& opts);
    virtual int get_n_forward() const;
    ///@}

    ///@{
    /** \brief Reverse mode derivatives */
    virtual Function get_reverse(const std::string& name, int nadj,
                                 const std::vector<std::string>& i_names,
                                 const std::vector<std::string>& o_names,
                                 const Dict& opts);
    virtual int get_n_reverse() const;
    ///@}

    ///@{
    /** \brief Full Jacobian */
    virtual bool hasFullJacobian() const;
    virtual Function getFullJacobian(const std::string& name, const Dict& opts);
    ///@}
  };

  class CASADI_EXPORT SimplifiedExternal : public External {
  public:
    /** \brief Constructor */
    SimplifiedExternal(const std::string& name, const Importer& li);

    /** \brief  Destructor */
    virtual ~SimplifiedExternal() { this->clear_memory();}

    /// Initialize
    virtual void init(const Dict& opts);

    /** \brief Use simplified signature */
    virtual bool simplifiedCall() const { return true;}

    /// @{
    /** \brief Retreive sparsities */
    virtual Sparsity get_sparsity_in(int i) { return Sparsity::scalar();}
    virtual Sparsity get_sparsity_out(int i) { return Sparsity::scalar();}
    /// @}
  };

  class CASADI_EXPORT GenericExternal : public External {
    // Sparsities
    sparsity_t sparsity_in_, sparsity_out_;

    // Maximum number of memory objects
    int n_mem_;

  public:
    /** \brief Constructor */
    GenericExternal(const std::string& name, const Importer& li);

    /** \brief  Destructor */
    virtual ~GenericExternal() { this->clear_memory();}

    /// Initialize
    virtual void init(const Dict& opts);

    /// @{
    /** \brief Retreive sparsities */
    virtual Sparsity get_sparsity_in(int i);
    virtual Sparsity get_sparsity_out(int i);
    /// @}

    /** \brief Maximum number of memory objects */
    virtual int n_mem() const { return n_mem_;}
  };


} // namespace casadi
/// \endcond

#endif // CASADI_EXTERNAL_IMPL_HPP
