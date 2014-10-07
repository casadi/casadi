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


#ifndef CASADI_LINEAR_SOLVER_INTERNAL_HPP
#define CASADI_LINEAR_SOLVER_INTERNAL_HPP

#include "linear_solver.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL

namespace casadi {

  /** Internal class
      @copydoc LinearSolver_doc
  */
  class CASADI_CORE_EXPORT
  LinearSolverInternal : public FunctionInternal,
                         public PluginInterface<LinearSolverInternal> {
  public:
    /// Constructor
    LinearSolverInternal(const Sparsity& sparsity, int nrhs);

    /// Destructor
    virtual ~LinearSolverInternal();

    /// Clone
    virtual LinearSolverInternal* clone() const { return new LinearSolverInternal(*this);}

    /// Initialize
    virtual void init();

    /// Solve the system of equations
    virtual void evaluate();

    /// Prepare the factorization
    virtual void prepare() {}

    /// Solve the system of equations, using internal vector
    virtual void solve(bool transpose);

    /// Solve the system of equations
    virtual void solve(double* x, int nrhs, bool transpose);

    /// Create a solve node
    MX solve(const MX& A, const MX& B, bool transpose);

    /// Evaluate numerically, possibly transposed
    virtual void evaluateDGen(const DMatrixPtrV& input, DMatrixPtrV& output, bool tr);

    /// Evaluate MX, possibly transposed
    virtual void evaluateSXGen(const SXPtrV& input, SXPtrV& output, bool tr);

    /// Evaluate MX, possibly transposed
    virtual void evaluateMXGen(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                               MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                               bool output_given, bool tr);

    /// Propagate sparsity, possibly transposed
    void propagateSparsityGen(DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                              std::vector<double>& rtmp, bool fwd, bool transpose);

    ///@{
    /// Propagate sparsity through a linear solve
    void spSolve(bvec_t* X, const bvec_t* B, bool transpose) const;
    void spSolve(DMatrix& X, const DMatrix& B, bool transpose) const;
    ///@}


    /// Solve the system of equations <tt>Lx = b</tt>
    virtual void solveL(double* x, int nrhs, bool transpose);

    /// Obtain a symbolic Cholesky factorization
    virtual Sparsity getFactorizationSparsity(bool transpose) const;

    /// Obtain a numeric Cholesky factorization
    virtual DMatrix getFactorization(bool transpose) const;

    /// Dulmage-Mendelsohn decomposition
    std::vector<int> rowperm_, colperm_, rowblock_, colblock_;

    /// Is prepared
    bool prepared_;

    /// Get sparsity pattern
    int nrow() const { return input(LINSOL_A).size1();}
    int ncol() const { return input(LINSOL_A).size2();}
    int nnz() const { return input(LINSOL_A).size();}
    const std::vector<int>& row() const { return input(LINSOL_A).row();}
    const std::vector<int>& colind() const { return input(LINSOL_A).colind();}

    // Creator function for internal class
    typedef LinearSolverInternal* (*Creator)(const Sparsity& sp, int nrhs);

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;
  };


} // namespace casadi
/// \endcond

#endif // CASADI_LINEAR_SOLVER_INTERNAL_HPP

