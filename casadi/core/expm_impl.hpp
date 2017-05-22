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
    /** \brief Number of function inputs and outputs */
    size_t get_n_in() override { return 2;}
    size_t get_n_out() override { return 1;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    Sparsity get_sparsity_in(int i) override;
    Sparsity get_sparsity_out(int i) override;
    /// @}

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize
    void init(const Dict& opts) override;

    ///@{
    /** \brief Generate a function that calculates \a nfwd forward derivatives */
    Function get_forward(const std::string& name, int nfwd,
                                 const std::vector<std::string>& i_names,
                                 const std::vector<std::string>& o_names,
                                 const Dict& opts) const override;
    int get_n_forward() const override { return 64;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives */
    Function get_reverse(const std::string& name, int nadj,
                                 const std::vector<std::string>& i_names,
                                 const std::vector<std::string>& o_names,
                                 const Dict& opts) const override;
    int get_n_reverse() const override { return 64;}
    ///@}

    /// Generate the sparsity of a Jacobian block
    Sparsity getJacSparsity(int iind, int oind, bool symmetric) const override;

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

    /** \brief Get type name */
    std::string type_name() const override {
      return std::string("expm_") + plugin_name();
    }

    /** \brief Get default input value */
    double default_in(int ind) const override;

  protected:
    Sparsity A_;
    bool const_A_;

  };


} // namespace casadi
/// \endcond
#endif // CASADI_EXPM_IMPL_HPP
