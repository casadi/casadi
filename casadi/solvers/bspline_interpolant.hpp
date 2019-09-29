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

    static std::vector<double> greville_points(const std::vector<double>& x, casadi_int deg);

    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new BSplineInterpolant(s); }

  protected:
     /** \brief Deserializing constructor */
    explicit BSplineInterpolant(DeserializingStream& s);

    static std::vector<double> not_a_knot(const std::vector<double>& x, casadi_int k);

    template <typename M>
    MX construct_graph(const MX& x, const M& values, const Dict& opts);

    enum FittingAlgorithm {ALG_NOT_A_KNOT, ALG_SMOOTH_LINEAR};

    /// Only used during init, no need to serialize these
    std::string linear_solver_;
    FittingAlgorithm algorithm_;
    double smooth_linear_frac_;
    std::vector<casadi_int> degree_;
  };


  template <typename M>
  MX BSplineInterpolant::construct_graph(const MX& x, const M& values, const Dict& opts) {

    std::vector< std::vector<double> > grid;
    for (casadi_int k=0;k<degree_.size();++k) {
      std::vector<double> local_grid(grid_.begin()+offset_[k], grid_.begin()+offset_[k+1]);
      grid.push_back(local_grid);
    }

    Dict opts_bspline;
    opts_bspline["lookup_mode"] = lookup_modes_;

    switch (algorithm_) {
      case ALG_NOT_A_KNOT:
        {
          std::vector< std::vector<double> > knots;
          for (casadi_int k=0;k<degree_.size();++k)
            knots.push_back(not_a_knot(grid[k], degree_[k]));
          Dict opts_dual;
          opts_dual["lookup_mode"] = lookup_modes_;

          DM J = MX::bspline_dual(meshgrid(grid), knots, degree_, opts_dual);

          casadi_assert_dev(J.size1()==J.size2());

          M V = M::reshape(values, m_, -1).T();
          M C_opt = solve(J, V, linear_solver_);

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
            linear = interpolant("linear", "linear", grid, m_);
          } else {
            linear = interpolant("linear", "linear", grid, values_);
          }

          std::vector< std::vector<double> > egrid;
          std::vector< std::vector<double> > new_grid;

          for (casadi_int k=0;k<n_dim;++k) {
            casadi_assert(degree_[k]==3, "Only degree 3 supported for 'smooth_linear'.");

            // Add extra knots
            const std::vector<double>& g = grid[k];

            // Determine smallest gap.
            double m = inf;
            for (casadi_int i=0;i<g.size()-1;++i) {
              double delta = g[i+1]-g[i];
              if (delta<m) m = delta;
            }
            double step = smooth_linear_frac_*m;

            // Add extra knots
            std::vector<double> new_g;
            new_g.push_back(g.front());
            new_g.push_back(g.front()+step);
            for (casadi_int i=1;i<g.size()-1;++i) {
              new_g.push_back(g[i]-step);
              new_g.push_back(g[i]);
              new_g.push_back(g[i]+step);
            }
            new_g.push_back(g.back()-step);
            new_g.push_back(g.back());
            new_grid.push_back(new_g);

            // Correct multiplicity
            double v1 = new_g.front();
            double vend = new_g.back();
            new_g.insert(new_g.begin(), degree_[k], v1);
            new_g.insert(new_g.end(), degree_[k], vend);

            grid[k] = new_g;

            // Compute greville points
            egrid.push_back(greville_points(new_g, degree_[k]));
          }

          std::vector<double> mg = meshgrid(egrid);
          casadi_int N = mg.size()/n_dim;

          // Evaluate linear interpolation on greville grid
          DM arg = DM::reshape(mg, n_dim, N);
          std::vector<M> res;
          if (has_parametric_values()) {
            res = linear(std::vector<M>{M(arg), values});
          } else {
            res = linear(std::vector<M>{M(arg)});
          }

          return MX::bspline(x, res[0], grid, degree_, m_, opts_bspline);
        }
      default:
        casadi_assert_dev(false);
      }
    }

} // namespace casadi

/// \endcond
#endif // CASADI_BSPLINE_INTERPOLANT_HPP
