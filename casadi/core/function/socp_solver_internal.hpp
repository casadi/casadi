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
  class CASADI_CORE_EXPORT
  SocpSolverInternal : public FunctionInternal,
                       public PluginInterface<SocpSolverInternal> {
  public:

    // Constructor
    SocpSolverInternal(const std::vector<Sparsity>& st);

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
    void printProblem(std::ostream &stream=std::cout) const;

    // Creator function for internal class
    typedef SocpSolverInternal* (*Creator)(const SOCPStructure& st);

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
};


} // namespace casadi

/// \endcond

#endif // CASADI_SOCP_SOLVER_INTERNAL_HPP

