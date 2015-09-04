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


#ifndef CASADI_SOCP_SOLVER_INTERNAL_HPP
#define CASADI_SOCP_SOLVER_INTERNAL_HPP

#include "socp_solver.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL

namespace casadi {

  /// Internal class
  class CASADI_EXPORT
  SocpSolverInternal : public FunctionInternal,
                       public PluginInterface<SocpSolverInternal> {
  public:

    // Constructor
    SocpSolverInternal(const std::map<std::string, Sparsity>& st);

    // Destructor
    virtual ~SocpSolverInternal() = 0;

    // Initialize
    virtual void init();

    // Solve the system of equations
    virtual void evaluate();

    // Solve the system of equations
    virtual void solve();

    /// \brief Check if the numerical values of the supplied bounds make sense
    virtual void checkInputs() const;

    /// Print out problem statement for debugging
    void printProblem(std::ostream &stream=casadi::userOut()) const;

    // Creator function for internal class
    typedef SocpSolverInternal* (*Creator)(const std::map<std::string, Sparsity>& st);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "socp";}

  protected:

    /// Problem structure
    std::vector<Sparsity> st_;

    /// Sizes of each block
    std::vector<int> ni_;

    /// Total size of G
    int N_;

    /// Size of decision variable vector
    int n_;

    /// The number SOC constraints
    int m_;

    /// Number of linear constraints
    int nc_;

    /// Indicates whether problem is printed before solving
    bool print_problem_;

    /** \brief Convert Second Order Cone Programming (SOCP) problem in standard form to one in a conic form.
    \verbatim
    Dual:

    max          c' x
    x
    subject to
    || yi ||_2 <= ti  i = 1..m
    A x + b == 0
    lbx <= x
     
    Dimensions                       | Meaning in terms of primal variables
    ------------------------------------------------------------------------------------------------------------
    with x ( nx x 1)                 | [lag_ti' lag_yi' lag_lba' lag_uba' lag_lbx' lag_ubx']'
    c dense ( nx x 1 )               | [-fi' hi' LBA' -UBA' LBX' -UBX']'
    yi dense ( ni )                  | Lagrange multipliers for equality constraints: yi = Gi' x + hi
    ti dense ( m )                   | Lagrange multipliers for equality constraint: ti = ei' x + fi
    A  sparse ( n x nx )             | [-ei Gi -Alba' Auba' -Ilbx Iubx]
    b  dense ( n x 1 )               | [c]
    lbx dense ( nx x 1 )             | [-inf' 0']'
    nx = m + N + nlba + nuba + nlbx + nubx (alias: dual_n_)
    \endverbatim
    */
    void convertToDualSocp();

    /// Linear objective in dual SOCP
    std::vector<double>  dual_c_;

    /// Sparse representation linear inequality constraints in dual SOCP
    /// @{
    std::vector<double>  dual_A_data_;
    std::vector<int>     dual_A_row_;
    std::vector<int>     dual_A_colind_;
    /// @}

    /// Vector of affine terms in linear inequality constraint in dual SOCP
    std::vector<double>  dual_b_;

    /// Size of dual decision variable vector
    int dual_n_;

    /// Number of linear constraints in dual problem
    int dual_nc_;

    /// Transpose of A matrix in primal problem definition (fixed sparsity pattern)
    DMatrix           primal_A_T_;
    std::vector<int>  primal_A_T_temp_int_;

    /** Indices of lower bounded linear inequality constraints (LBA != -inf),
    used to set up dual SOCP variables */
    std::vector<int>  primal_idx_lba_;

    /** Indices of upper bounded linear inequality constraints (UBA != inf),
    used to set up dual SOCP variables */
    std::vector<int>  primal_idx_uba_;

    /// Indices of simple lower bounds  (LBX != -inf), used to set up dual SOCP variables
    std::vector<int>  primal_idx_lbx_;

    /// Indices of simple upper bounds (UBX != inf), used to set up dual SOCP variables
    std::vector<int>  primal_idx_ubx_;
};


} // namespace casadi

/// \endcond

#endif // CASADI_SOCP_SOLVER_INTERNAL_HPP

