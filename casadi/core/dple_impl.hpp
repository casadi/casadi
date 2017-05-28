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


#ifndef CASADI_DPLE_IMPL_HPP
#define CASADI_DPLE_IMPL_HPP

#include "dple.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL
namespace casadi {
  /// Internal class
  class CASADI_EXPORT Dple : public FunctionInternal, public PluginInterface<Dple> {
  public:

    // Constructor
    Dple(const std::string& name, const SpDict &st);

    // Destructor
    ~Dple() override = 0;

    ///@{
    /** \brief Number of function inputs and outputs */
    size_t get_n_in() override { return DPLE_NUM_IN;}
    size_t get_n_out() override { return DPLE_NUM_OUT;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    Sparsity get_sparsity_in(int i) override;
    Sparsity get_sparsity_out(int i) override;
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    std::string get_name_in(int i) override { return dple_in(i);}
    std::string get_name_out(int i) override { return dple_out(i);}
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
    Function get_forward(int nfwd, const std::string& name,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const override;
    bool has_forward(int nfwd) const override { return true;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives */
    Function get_reverse(const std::string& name, int nadj,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const override;
    bool has_reverse(int nadj) const override { return true;}
    ///@}

    // Creator function for internal class
    typedef Dple* (*Creator)(const std::string& name,
                              const std::map<std::string, Sparsity>& st);

    // No static functions exposed
    struct Exposed{   };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "dple";}

    /** \brief Get type name */
    std::string type_name() const override {
      return std::string("dple_") + plugin_name();
    }

    /** \brief Get default input value */
    double default_in(int ind) const override;

  protected:

    /// List of sparsities of A_i
    Sparsity A_;

    /// List of sparsities of V_i
    Sparsity V_;

    /// Period
    int K_;

    /// Constant dimensions
    bool const_dim_;

    /// Assume positive definiteness of P_i
    bool pos_def_;

    /// Throw an error when system is unstable
    bool error_unstable_;

    /// Margin for instability detection
    double eps_unstable_;

    /// Number of right hand sides
    int nrhs_;

  };


} // namespace casadi
/// \endcond
#endif // CASADI_DPLE_IMPL_HPP
