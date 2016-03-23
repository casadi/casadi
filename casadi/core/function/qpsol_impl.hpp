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


#ifndef CASADI_QPSOL_IMPL_HPP
#define CASADI_QPSOL_IMPL_HPP

#include "qpsol.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL
namespace casadi {
  /// Internal class
  class CASADI_EXPORT Qpsol : public FunctionInternal, public PluginInterface<Qpsol> {
  public:

    // Constructor
    Qpsol(const std::string& name, const std::map<std::string, Sparsity> &st);

    // Destructor
    virtual ~Qpsol() = 0;

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() { return QPSOL_NUM_IN;}
    virtual size_t get_n_out() { return QPSOL_NUM_OUT;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int i);
    virtual Sparsity get_sparsity_out(int i);
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    virtual std::string get_name_in(int i) { return qpsol_in(i);}
    virtual std::string get_name_out(int i) { return qpsol_out(i);}
    /// @}

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    // Initialize
    virtual void init(const Dict& opts);

    /// \brief Check if the numerical values of the supplied bounds make sense
    virtual void checkInputs(const double* lbx, const double* ubx,
                             const double* lba, const double* uba) const;

    /** Generate native code in the interfaced language for debugging */
    virtual void generateNativeCode(std::ostream& file) const;

    // Creator function for internal class
    typedef Qpsol* (*Creator)(const std::string& name,
                              const std::map<std::string, Sparsity>& st);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "qpsol";}

    /** \brief Get default input value */
    virtual double default_in(int ind) const;

    /// Can discrete variables be treated
    virtual bool integer_support() const { return false;}

  protected:
    /// Options
    std::vector<bool> discrete_;

    /// Problem structure
    Sparsity H_, A_;

    /// Number of decision variables
    int n_;

    /// The number of constraints (counting both equality and inequality) == A.size1()
    int nc_;
  };


} // namespace casadi
/// \endcond
#endif // CASADI_QPSOL_IMPL_HPP
