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

#ifndef CASADI_SHF_INTERPOLANT_HPP
#define CASADI_SHF_INTERPOLANT_HPP

#include "casadi/core/interpolant_impl.hpp"
#include "casadi/core/shf_spline.hpp"
#include <casadi/solvers/casadi_interpolant_shf_export.h>

/** \defgroup plugin_Interpolant_shf Title
    \par

    Smooth hat function (SHF) approximation of the piecewise-linear interpolant:
    C^k, coincides with the look-up table outside epsilon-balls around the interior
    grid points, no fitting step. Wraps shf_spline.
*/

/** \pluginsection{Interpolant,shf} */

/// \cond INTERNAL

namespace casadi {

  /** \brief \pluginbrief{Interpolant,shf}

    @copydoc Interpolant_doc
    @copydoc plugin_Interpolant_shf
    \author Joris Gillis
    \date 2026
  */
  class CASADI_INTERPOLANT_SHF_EXPORT SHFInterpolant : public Interpolant {
  public:
    // Constructor
    SHFInterpolant(const std::string& name,
                   const std::vector<double>& grid,
                   const std::vector<casadi_int>& offset,
                   const std::vector<double>& values,
                   casadi_int m);

    // Destructor
    ~SHFInterpolant() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "shf";}

    // Get name of the class
    std::string class_name() const override { return "SHFInterpolant";}

    /** \brief  Create a new Interpolant */
    static Interpolant* creator(const std::string& name,
                                const std::vector<double>& grid,
                                const std::vector<casadi_int>& offset,
                                const std::vector<double>& values,
                                casadi_int m) {
      return new SHFInterpolant(name, grid, offset, values, m);
    }

    // Initialize
    void init(const Dict& opts) override;

    /// Evaluate numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    /// Evaluate symbolically: inline the wrapper, so only the shf_spline calls survive in the graph
    void eval_mx(const MXVector& arg, MXVector& res,
                 bool always_inline, bool never_inline) const override;

    ///@{
    /** \brief Full Jacobian */
    bool has_jacobian() const override { return true;}
    Function get_jacobian(const std::string& name,
                          const std::vector<std::string>& inames,
                          const std::vector<std::string>& onames,
                          const Dict& opts) const override;
    ///@}

    ///@{
    /** \brief Return function that calculates forward derivatives */
    bool has_forward(casadi_int nfwd) const override { return true; }
    Function get_forward(casadi_int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    ///@{
    /** \brief Return function that calculates adjoint derivatives */
    bool has_reverse(casadi_int nadj) const override { return true; }
    Function get_reverse(casadi_int nadj, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    /** \brief Is codegen supported? */
    bool has_codegen() const override { return true;}

    /** \brief Generate code for the body of the C function */
    void codegen_body(CodeGenerator& g) const override;

    /** \brief Generate code for the declarations of the C function */
    void codegen_declarations(CodeGenerator& g) const override;

    /// Inputs: x, the values when parametric, epsilon when parametric
    size_t get_n_in() override { return Interpolant::get_n_in() + epsilon_parametric_; }
    Sparsity get_sparsity_in(casadi_int i) override;
    std::string get_name_in(casadi_int i) override;

    /// Index of the epsilon input
    casadi_int arg_epsilon() const { return 1 + has_parametric_values() + has_parametric_grid(); }

    /// A documentation string
    static const std::string meta_doc;

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Spline Function
    Function S_;

    // Get all embedded functions, recursively
    void find(std::map<FunctionInternal*, std::pair<Function, size_t> >& all_fun,
      casadi_int max_depth) const override;

    /** \brief  Propagate sparsity forward */
    int sp_forward(const bvec_t** arg, bvec_t** res,
                    casadi_int* iw, bvec_t* w, void* mem) const override {
      return S_->sp_forward(arg, res, iw, w, mem);
    }

    /** \brief  Propagate sparsity backwards */
    int sp_reverse(bvec_t** arg, bvec_t** res,
        casadi_int* iw, bvec_t* w, void* mem) const override {
      return S_->sp_reverse(arg, res, iw, w, mem);
    }

    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new SHFInterpolant(s); }

  protected:
     /** \brief Deserializing constructor */
    explicit SHFInterpolant(DeserializingStream& s);

    /// Whether epsilon is an input rather than a constant
    bool epsilon_parametric_;

    /// Only used during init, no need to serialize these
    casadi_int smoothness_order_;
    std::vector<double> epsilon_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_SHF_INTERPOLANT_HPP
