/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef SOLVE_HPP
#define SOLVE_HPP

#include "mx_node.hpp"
#include "../fx/linear_solver.hpp"

namespace CasADi{
  /** \brief An MX atomic for linear solver solution: x = r * A^-1 or x = r * A^-T
      
      Forward derivatives:
      x_dot = (r_dot - x * A_dot) * A^-1

      Adjoint derivatives:
      r_bar = x_bar * A^-T
      A_bar = -x^T * r_bar

      \author Joel Andersson
      \date 2013
  */
  template<bool Tr>
  class Solve : public MXNode{
  public:
    
    /** \brief  Constructor */
    Solve(const MX& r, const MX& A, const LinearSolver& linear_solver);

    /** \brief  Destructor */
    virtual ~Solve(){}

    /** \brief  Clone function */
    virtual Solve* clone() const{ return new Solve(*this);}

    /** \brief  Print expression (make sure number of calls is not exceeded) */
    virtual void print(std::ostream &stream, long& remaining_calls) const;

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /// Evaluate the function numerically
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp);

    /** \brief  Evaluate the function symbolically (MX) */
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp, bool fwd);

    /** \brief Get the operation */
    virtual int getOp() const{ return OP_SOLVE;}

    /// Can the operation be performed inplace (i.e. overwrite the result)
    virtual int numInplace() const{ return 1;}

    /** \brief  Get function reference */
    virtual FX& getFunction(){ return linear_solver_;}

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

    /// Get number of temporary variables needed
    virtual void nTmp(size_t& ni, size_t& nr){ ni=0; nr=sparsity().size1();}

    /// Linear Solver (may be shared between multiple nodes)
    LinearSolver linear_solver_;
  };


} // namespace CasADi


#endif // SOLVE_HPP
