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


#ifndef CASADI_BSPLINE_INTERPOLANT_HPP
#define CASADI_BSPLINE_INTERPOLANT_HPP

#include "casadi/core/interpolant_impl.hpp"
#include <casadi/solvers/casadi_interpolant_bspline_export.h>

/** \defgroup plugin_Interpolant_bspline
*/

/** \pluginsection{Interpolant,bspline} */

/// \cond INTERNAL

namespace casadi {
  class BSplineCommon;
  /** \brief \pluginbrief{Interpolant,bspline}

    N-dimensional BSpline interpolator

    Uses not-a-knot conditions.
    For 1D and 2D cases, this code is equivalent to fitpack

    @copydoc Interpolant_doc
    @copydoc plugin_Interpolant_bspline
    \author Joris Gillis
    \date 2017
  */
  class CASADI_INTERPOLANT_BSPLINE_EXPORT BSplineInterpolant : public Interpolant {
  public:
    // Constructor
    BSplineInterpolant(const std::string& name,
                      const std::vector<double>& grid,
                      const std::vector<casadi_int>& offset,
                      const std::vector<double>& values,
                      casadi_int m);

    // Destructor
    ~BSplineInterpolant() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "bspline";}

    // Get name of the class
    std::string class_name() const override { return "BSplineInterpolant";}

    /** \brief  Create a new Interpolant */
    static Interpolant* creator(const std::string& name,
                                const std::vector<double>& grid,
                                const std::vector<casadi_int>& offset,
                                const std::vector<double>& values,
                                casadi_int m) {
      return new BSplineInterpolant(name, grid, offset, values, m);
    }

    // Is differentiable? Deferred to bspline
    bool is_diff_in(casadi_int i) override { return true; }

    // Initialize
    void init(const Dict& opts) override;

    /// Evaluate numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    ///@{
    /** \brief Full Jacobian */
    bool has_jacobian() const override { return true;}
    Function get_jacobian(const std::string& name,
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

    /// A documentation string
    static const std::string meta_doc;

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Spline Function
    Function S_;

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

    static DM greville_points(const DM& x, casadi_int deg);

    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new BSplineInterpolant(s); }

  protected:
     /** \brief Deserializing constructor */
    explicit BSplineInterpolant(DeserializingStream& s);

    static std::vector<double> not_a_knot(const std::vector<double>& x, casadi_int k);

    template <typename M, typename Mk>
    MX construct_graph(const MX& x, const M& values, const Mk& g, const Dict& linsol_options, const Dict& opts);

    enum FittingAlgorithm {ALG_NOT_A_KNOT, ALG_SMOOTH_LINEAR};

    /// Only used during init, no need to serialize these
    std::string linear_solver_;
    FittingAlgorithm algorithm_;
    double smooth_linear_frac_;
    std::vector<casadi_int> degree_;
  };


  template <typename M, typename Mk>
  MX BSplineInterpolant::construct_graph(const MX& x, const M& values, const Mk& g,
      const Dict& linsol_options, const Dict& opts) {

    std::vector< Mk > grid;
    for (casadi_int k=0;k<degree_.size();++k) {
      grid.push_back(g(range(offset_[k], offset_[k+1])));
    }

    bool do_inline = false;
    for (auto&& op : opts) {
      if (op.first=="inline") {
        do_inline = op.second;
      }
    }

    Dict opts_bspline;
    opts_bspline["lookup_mode"] = lookup_modes_;
    opts_bspline["inline"] = do_inline;

    switch (algorithm_) {
      case ALG_NOT_A_KNOT:
        {
          std::vector< std::vector<double> > knots;
          for (casadi_int k=0;k<degree_.size();++k)
            knots.push_back(not_a_knot(Interpolant::parse_grid(grid[k]), degree_[k]));
          Dict opts_dual;
          opts_dual["lookup_mode"] = lookup_modes_;

          DM J = MX::bspline_dual(meshgrid(Interpolant::parse_grid(grid)), knots, degree_, opts_dual);

          casadi_assert_dev(J.size1()==J.size2());

          M V = M::reshape(values, m_, -1).T();
          M C_opt = solve(J, V, linear_solver_, linsol_options);

          if (!has_parametric_values()) {
            double fit = static_cast<double>(norm_1(mtimes(J, C_opt) - V));
            if (verbose_) casadi_message("Lookup table fitting error: " + str(fit));
          }

          return MX::bspline(x, C_opt.T(), knots, degree_, m_, opts_bspline);
        }
      case ALG_SMOOTH_LINEAR:
        {
          casadi_int n_dim = degree_.size();
          // Linear fit
          Function linear;
          if (has_parametric_values()) {
            linear = interpolant("linear", "linear", Interpolant::parse_grid(grid), m_);
          } else {
            linear = interpolant("linear", "linear", Interpolant::parse_grid(grid), values_);
          }

          std::vector< Mk > egrid;
          std::vector< Mk > new_grid;

          for (casadi_int k=0;k<n_dim;++k) {
            casadi_assert(degree_[k]==3, "Only degree 3 supported for 'smooth_linear'.");

            // Add extra knots
            Mk& g = grid[k];

            // Determine smallest gap.
            Mk m = mmin(g(range(1, g.numel()))-g(range(g.numel()-1)));
            Mk step = smooth_linear_frac_*m;

            // Add extra knots
            std::vector<Mk> new_g;

            std::vector<Mk> g_parts = vertsplit(g,{0, 1, g.numel()-1, g.numel()});
            Mk g0 = g_parts[0];
            Mk gm = g_parts[1];
            Mk gN = g_parts[2];
            g0 = vertcat(repmat(g0, degree_[k]+1, 1), g0+step);
            gm = vec(horzcat(gm-step, gm, gm+step).T());
            gN = vertcat(gN-step, repmat(gN, degree_[k]+1, 1));

            g = vertcat(g0, gm, gN);

            // Compute greville points
            egrid.push_back(greville_points(g, degree_[k]));
          }

          Mk mg = meshgrid(egrid);
          uout() << mg << std::endl;
          casadi_int N = mg.numel()/n_dim;

          // Evaluate linear interpolation on greville grid
          Mk arg = Mk::reshape(mg, n_dim, N);
          std::vector<M> res;
          if (has_parametric_values()) {
            res = linear(std::vector<M>{M(arg), values});
          } else {
            res = linear(std::vector<M>{M(arg)});
          }

          return MX::bspline(x, res[0], Interpolant::parse_grid(grid), degree_, m_, opts_bspline);
        }
      default:
        casadi_assert_dev(false);
      }
    }

} // namespace casadi

/// \endcond
#endif // CASADI_BSPLINE_INTERPOLANT_HPP
