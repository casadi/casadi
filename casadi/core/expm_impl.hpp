/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
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

#ifndef CASADI_EXPM_IMPL_HPP
#define CASADI_EXPM_IMPL_HPP

#include "expm.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL
namespace casadi {
  /// Internal class
  class CASADI_EXPORT Expm : public FunctionInternal, public PluginInterface<Expm> {
  public:

    // Constructor
    Expm(const std::string& name, const Sparsity& A);
    // Destructor
    ~Expm() override = 0;

    ///@{
    /** \brief Number of function inputs and outputs

        \identifier{123} */
    size_t get_n_in() override { return 2;}
    size_t get_n_out() override { return 1;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs

        \identifier{124} */
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;
    /// @}

    ///@{
    /** \brief Options

        \identifier{125} */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize
    void init(const Dict& opts) override;

    ///@{
    /** \brief Generate a function that calculates \a nfwd forward derivatives

        \identifier{126} */
    Function get_forward(casadi_int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    bool has_forward(casadi_int nfwd) const override { return true;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives

        \identifier{127} */
    Function get_reverse(casadi_int nadj, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    bool has_reverse(casadi_int nadj) const override { return true;}
    ///@}

    /** \brief Which inputs are differentiable

        \identifier{128} */
    bool get_diff_in(casadi_int i) override { return i != 0 || !const_A_;}

    /// Generate the sparsity of a Jacobian block
    Sparsity get_jac_sparsity(casadi_int oind, casadi_int iind, bool symmetric) const override;
    bool has_jac_sparsity(casadi_int oind, casadi_int iind) const override { return true;}

    // Creator function for internal class
    typedef Expm* (*Creator)(const std::string& name, const Sparsity& A);

    // No static functions exposed
    struct Exposed{   };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "expm";}

  protected:
    Sparsity A_;
    bool const_A_;

  };


} // namespace casadi
/// \endcond
#endif // CASADI_EXPM_IMPL_HPP
