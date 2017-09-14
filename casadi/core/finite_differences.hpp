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


#ifndef CASADI_FINITE_DIFFERENCES_HPP
#define CASADI_FINITE_DIFFERENCES_HPP

#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  /** Calculate derivative using central differences
    * The algorithm is based on the package specification for TD12 in the HSL archive
    * \author Joel Andersson
    * \date 2017
  */
  class CASADI_EXPORT CentralDiff : public FunctionInternal {
  public:
    // Create function (use instead of constructor)
    static Function create(const std::string& name, int n, const Dict& opts);

    /** \brief Destructor */
    ~CentralDiff() override;

    /** \brief Get type name */
    std::string class_name() const override {return "CentralDiff";}

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    Sparsity get_sparsity_in(int i) override;
    Sparsity get_sparsity_out(int i) override;
    /// @}

    /** \brief Get default input value */
    double default_in(int ind) const override;

    ///@{
    /** \brief Number of function inputs and outputs */
    size_t get_n_in() override;
    size_t get_n_out() override;
    ///@}

    ///@{
    /** \brief Names of function input and outputs */
    std::string get_name_in(int i) override;
    std::string get_name_out(int i) override;
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    // Evaluate numerically
    int eval(const double** arg, double** res, int* iw, double* w, void* mem) const override;

    /** \brief Is the scheme using the (nondifferentiated) output? */
    bool uses_output() const override {return true;}

    /** \brief Is codegen supported? */
    bool has_codegen() const override { return true;}

    /** \brief Generate code for the declarations of the C function */
    void codegen_declarations(CodeGenerator& g) const override;

    /** \brief Generate code for the body of the C function */
    void codegen_body(CodeGenerator& g) const override;

    ///@{
    /** \brief Second order derivatives */
    bool has_forward(int nfwd) const override { return true;}
    Function get_forward(int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;

  protected:
    // Constructor (protected, use create function)
    CentralDiff(const std::string& name, int n);

    // Dimensions
    int n_z_, n_y_;

    // Number of directional derivatives
    int n_;

    // Maximum allowed precision
    double h_max_;

    // Input precision
    double eps_in_;

    // Output precision
    double eps_out_;

    // Ratio of roundoff error to truncation error
    double u_aim_;

    // Perturbation
    double h_, h2_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_FINITE_DIFFERENCES_HPP
