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
    virtual ~Dple() = 0;

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() { return DPLE_NUM_IN;}
    virtual size_t get_n_out() { return DPLE_NUM_OUT;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int i);
    virtual Sparsity get_sparsity_out(int i);
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    virtual std::string get_name_in(int i) { return dple_in(i);}
    virtual std::string get_name_out(int i) { return dple_out(i);}
    /// @}

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    // Initialize
    virtual void init(const Dict& opts);

    ///@{
    /** \brief Generate a function that calculates \a nfwd forward derivatives */
    virtual Function get_forward(const std::string& name, int nfwd,
                                 const std::vector<std::string>& i_names,
                                 const std::vector<std::string>& o_names,
                                 const Dict& opts);
    virtual int get_n_forward() const { return 64;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives */
    virtual Function get_reverse(const std::string& name, int nadj,
                                 const std::vector<std::string>& i_names,
                                 const std::vector<std::string>& o_names,
                                 const Dict& opts);
    virtual int get_n_reverse() const { return 64;}
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

    /** \brief Get default input value */
    virtual double default_in(int ind) const;

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
