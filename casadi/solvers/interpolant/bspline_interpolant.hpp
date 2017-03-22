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

#include "casadi/core/function/interpolant_impl.hpp"
#include <casadi/solvers/interpolant/casadi_interpolant_bspline_export.h>

/** \defgroup plugin_Interpolant_bspline
*/

/** \pluginsection{Interpolant,bspline} */

/// \cond INTERNAL

namespace casadi {
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
                      const std::vector<int>& offset,
                      const std::vector<double>& values);

    // Destructor
    virtual ~BSplineInterpolant();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "bspline";}

    /** \brief  Create a new Interpolant */
    static Interpolant* creator(const std::string& name,
                                const std::vector<double>& grid,
                                const std::vector<int>& offset,
                                const std::vector<double>& values) {
      return new BSplineInterpolant(name, grid, offset, values);
    }

    // Initialize
    virtual void init(const Dict& opts);

    /// Evaluate numerically
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    ///@{
    /** \brief Full Jacobian */
    virtual bool hasFullJacobian() const { return true;}
    virtual Function getFullJacobian(const std::string& name,
                                      const std::vector<std::string>& i_names,
                                      const std::vector<std::string>& o_names,
                                      const Dict& opts);
    ///@}

    /** \brief Is codegen supported? */
    virtual bool has_codegen() const { return true;}

    /** \brief Generate code for the body of the C function */
    virtual void generateBody(CodeGenerator& g) const;

    /// A documentation string
    static const std::string meta_doc;

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /// Degree of the spline
    std::vector<int> degree_;

    /// Linear solvere
    std::string linear_solver_;

    // Spline Function
    Function S_;
  };


} // namespace casadi

/// \endcond
#endif // CASADI_BSPLINE_INTERPOLANT_HPP
