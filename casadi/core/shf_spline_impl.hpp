/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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

#ifndef CASADI_SHF_SPLINE_IMPL_HPP
#define CASADI_SHF_SPLINE_IMPL_HPP

#include "shf_spline.hpp"
#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {
  /** \brief Smooth Hat Function (SHF) spline

      Approximates the piecewise-linear interpolant of a look-up table by a
      k-times continuously differentiable spline that coincides with the table
      exactly outside epsilon-neighbourhoods of the interior grid points, and
      converges to it as epsilon goes to zero.

      The coefficients are the table values themselves; there is no fitting
      step. Support is two cells wide independently of k, so the evaluation
      stencil is 2 (padded to 3) per axis rather than 2k+2.

      Inputs are x, the table values C when parametric, and epsilon, one
      absolute half-width per axis that must stay below half the smallest grid
      spacing of the axis it applies to. Only x is differentiable.

      A function of derivative order p outputs every mixed partial of total
      order p (m-by-nb), all from one traversal of the coefficient stencil; its
      Jacobian is a selection out of the order p+1 function.

      \author Joris Gillis
      \date 2026
  */
  class CASADI_EXPORT SHFSplineFunction : public FunctionInternal {
  public:
    /// Constructor, fixed table values
    SHFSplineFunction(const std::string& name, const std::vector<double>& grid,
        const std::vector<casadi_int>& offset, const std::vector<double>& values,
        casadi_int k, casadi_int m, casadi_int order = 0);

    /// Constructor, parametric table values
    SHFSplineFunction(const std::string& name, const std::vector<double>& grid,
        const std::vector<casadi_int>& offset, casadi_int k, casadi_int m,
        casadi_int order = 0);

    /// Destructor
    ~SHFSplineFunction() override;

    /// Get type name
    std::string class_name() const override { return "SHFSplineFunction";}

    ///@{
    /// Options
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Initialize
    void init(const Dict& opts) override;

    ///@{
    /// Number of function inputs and outputs
    size_t get_n_in() override { return 2 + parametric_; }
    size_t get_n_out() override { return 1; }
    ///@}

    /// Only x is differentiable
    bool get_diff_in(casadi_int i) override { return i==0; }

    ///@{
    /// Sparsities of function inputs and outputs
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;
    ///@}

    ///@{
    /// Names of function input and outputs
    std::string get_name_in(casadi_int i) override;
    std::string get_name_out(casadi_int i) override;
    ///@}

    /// True when the table values arrive as an input rather than being stored
    bool parametric() const { return parametric_; }

    /// Index of the C input (only valid when parametric)
    casadi_int arg_c() const { return 1; }

    /// Index of the epsilon input
    casadi_int arg_eps() const { return 1 + parametric_; }

    /// Number of axes
    casadi_int n_dims() const { return offset_.size()-1; }

    /// Number of mixed partials at this order
    casadi_int nb() const { return multi_.size()/n_dims(); }

    /// Evaluate numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    /// Is codegen supported?
    bool has_codegen() const override { return true;}

    /// Generate code for the function body
    void codegen_body(CodeGenerator& g) const override;

    ///@{
    /// Jacobian of all outputs with respect to all inputs
    bool has_jacobian() const override { return true; }
    Function get_jacobian(const std::string& name,
                          const std::vector<std::string>& inames,
                          const std::vector<std::string>& onames,
                          const Dict& opts) const override;
    ///@}

    ///@{
    /// Forward and reverse sweeps: a few axpys on the next-order output
    bool has_forward(casadi_int nfwd) const override { return true; }
    bool has_reverse(casadi_int nadj) const override { return true; }
    Function get_forward(casadi_int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    Function get_reverse(casadi_int nadj, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    /// Same table at total derivative order order_+1, cached
    Function next_order() const;

    /// d f / d x_s for every axis s, each m-by-nb, all from one call of the next-order function
    std::vector<MX> jac_cols(const std::vector<MX>& arg) const;

    /// Selection matrix picking column multi_[b] + e_s out of the next-order output
    DM selector(casadi_int s) const;

    /// Multi-indices of total order p in n_dims axes, canonically ordered so that
    /// order 1 comes out as e_0, e_1, ... and the jacobian columns land in axis order
    static std::vector<casadi_int> multi_index(casadi_int n_dims, casadi_int p);

    static casadi_int get_coeff_size(casadi_int m, const std::vector<casadi_int>& offset);

    static void prepare(casadi_int m, const std::vector<casadi_int>& offset,
      casadi_int& coeffs_size, std::vector<casadi_int>& coeffs_dims,
      std::vector<casadi_int>& strides);

    /** Per-axis constants that are fixed once the grid is known: reciprocal
        spacings, the padded stencil width, and the smallest grid spacing
        (epsilon must stay below half of it for the epsilon-neighbourhoods of
        neighbouring grid points to stay disjoint). */
    static void derive(const std::vector<double>& grid,
      const std::vector<casadi_int>& offset, std::vector<double>& inv_h,
      std::vector<casadi_int>& width, std::vector<double>& min_h);

    /// order+1 rows of Bernstein control points: row p holds those of s_k^(p), the
    /// p-th derivative of the smoothing polynomial of order k, leading zeros dropped
    static std::vector<double> bernstein_ctrl(casadi_int k, casadi_int order);

    static size_t n_iw(casadi_int n_dims, casadi_int nb);
    static size_t n_w(casadi_int n_dims, casadi_int k, casadi_int order,
      casadi_int nb);

    /// Serialize an object without type information
    void serialize_body(SerializingStream &s) const override;

    /// Deserialize without type information
    static ProtoFunction* deserialize(DeserializingStream& s);

    /// String used to identify the immediate FunctionInternal subclass
    std::string serialize_base_function() const override { return "SHFSplineFunction"; }

    std::vector<double> grid_;
    std::vector<casadi_int> offset_;
    casadi_int k_;
    casadi_int m_;

    /// Table values; empty iff parametric()
    std::vector<double> values_;
    bool parametric_;

    /// Total derivative order of this function's output
    casadi_int order_;

    std::vector<std::string> lookup_modes_;
    std::vector<casadi_int> lookup_mode_;

    /// nb-by-n_dims table of multi-indices, all of total order order_
    std::vector<casadi_int> multi_;

    // Derived fields
    std::vector<double> inv_h_;
    std::vector<casadi_int> width_;
    std::vector<casadi_int> strides_;
    std::vector<casadi_int> coeffs_dims_;
    casadi_int coeffs_size_;
    std::vector<double> min_h_;
    std::vector<double> ctrl_;

  protected:
    /// Shared tail of the constructors
    void init_derived_members();

    /// Deserializing constructor
    explicit SHFSplineFunction(DeserializingStream& s);
  };

} // namespace casadi
/// \endcond

#endif // CASADI_SHF_SPLINE_IMPL_HPP
