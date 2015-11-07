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


#ifndef CASADI_OOQP_INTERFACE_HPP
#define CASADI_OOQP_INTERFACE_HPP

#include "casadi/core/function/qpsol.hpp"
#include <casadi/interfaces/ooqp/casadi_qpsol_ooqp_export.h>

/** \defgroup plugin_Qpsol_ooqp
 Interface to the OOQP Solver for quadratic programming
  The current implementation assumes that OOQP is configured with the MA27 sparse linear solver.

  NOTE: when doing multiple calls to evaluate(), check if you need to reInit();
*/

/** \pluginsection{Qpsol,ooqp} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{Qpsol,ooqp}

      @copydoc Qpsol_doc
      @copydoc plugin_Qpsol_ooqp

  */
  class CASADI_QPSOL_OOQP_EXPORT OoqpInterface : public Qpsol {
  public:
    /** \brief  Create a new Solver */
    explicit OoqpInterface(const std::string& name,
                           const std::map<std::string, Sparsity>& st);

    /** \brief  Create a new QP Solver */
    static Qpsol* creator(const std::string& name,
                                     const std::map<std::string, Sparsity>& st) {
      return new OoqpInterface(name, st);
    }

    /** \brief  Destructor */
    virtual ~OoqpInterface();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "ooqp";}

    /** \brief  Initialize */
    virtual void init();

    /// Solve the QP
    virtual void evalD(const double** arg, double** res, int* iw, double* w, void* mem);

    /// Throw error
    static const char* errFlag(int flag);

    /// Print an OOQP bounds vector
    static std::string printBounds(const std::vector<double>& b,
                                   const std::vector<char>& ib, int n, const char *sign);

    /// Problem data (vectors)
    std::vector<double> c_, bA_, xlow_, xupp_, clow_, cupp_, x_, gamma_, phi_, y_, z_, lambda_, pi_;

    /// Type of bounds
    std::vector<char> ixlow_, ixupp_, iclow_, icupp_;

    /// Problem data (matrices)
    std::vector<double> dQ_, dA_, dC_;

    /// Sparsities of matrices
    std::vector<int> irowQ_, jcolQ_, irowA_, jcolA_, irowC_, jcolC_;

    // Variable/constraint index
    std::vector<int> x_index_, c_index_;

    // Parameters
    std::vector<double> p_;

    // Transpose of linear constraints
    DMatrix AT_;
    std::vector<int> AT_tmp_;

    // Print level
    int print_level_;

    // Tolerances
    double mutol_, artol_;

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_OOQP_INTERFACE_HPP

