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

#ifndef CASADI_CPLEX_INTERFACE_HPP
#define CASADI_CPLEX_INTERFACE_HPP

#include "casadi/core/function/qpsol_impl.hpp"
#include <casadi/interfaces/cplex/casadi_qpsol_cplex_export.h>
#include "ilcplex/cplex.h"

#include <string>

/** \defgroup plugin_Qpsol_cplex

      Interface to Cplex solver for sparse Quadratic Programs
*/

/** \pluginsection{Qpsol,cplex} */

/// \cond INTERNAL

namespace casadi {

  struct CASADI_QPSOL_CPLEX_EXPORT CplexMemory : public Memory {
    /// Indicates if we have to warm-start
    bool is_warm;

    /// Nature of problem (always minimization)
    int objsen;

    /// Determines relation >,<, = in the linear constraints
    std::vector<char> sense;

    /// Coefficients of matrix A (constraint Jacobian)
    std::vector<int> matcnt;

    /// Right-hand side of constraints
    std::vector<double> rhs;

    /// Range of constraints
    std::vector<double> rngval;

    /// Coefficients of matrix H (objective Hessian)
    std::vector<int> qmatcnt;

    /// Storage for basis info of primal variables
    std::vector<int> cstat;

    /// Storage for basis info of slack variables
    std::vector<int> rstat;

    /// CPLEX environment
    CPXENVptr env;
    CPXLPptr lp;

    /// Constructor
    CplexMemory();

    /// Destructor
    virtual ~CplexMemory();
  };

  /** \brief \pluginbrief{Qpsol,cplex}

      @copydoc Qpsol_doc
      @copydoc plugin_Qpsol_cplex

      \author Attila Kozma, Joel Andersson
      \date 2012
  */
  class CASADI_QPSOL_CPLEX_EXPORT CplexInterface : public Qpsol {
  public:
    /** \brief  Create a new QP Solver */
    static Qpsol* creator(const std::string& name,
                                     const std::map<std::string, Sparsity>& st) {
      return new CplexInterface(name, st);
    }

    /// Constructor using sparsity patterns
    explicit CplexInterface(const std::string& name,
                            const std::map<std::string, Sparsity>& st);

    /// Destructor
    virtual ~CplexInterface();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "cplex";}

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    // Initialize the solver
    virtual void init(const Dict& opts);

    /** \brief Create memory block */
    virtual Memory* memory() const { return new CplexMemory();}

    /** \brief Initalize memory block */
    virtual void init_memory(Memory& mem) const;

    // Solve the QP
    virtual void eval(Memory& mem, const double** arg, double** res, int* iw, double* w) const;

    ///@{
    /// Options
    int qp_method_;
    bool dump_to_file_;
    std::string dump_filename_;
    double tol_;
    int dep_check_;
    int simplex_maxiter_, barrier_maxiter_;
    bool warm_start_;
    bool convex_;
    ///@}

    /// A documentation string
    static const std::string meta_doc;

  };
} // end namespace casadi
/// \endcond
#endif // CASADI_CPLEX_INTERFACE_HPP
