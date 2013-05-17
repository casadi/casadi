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

#include "nonlinear_solve.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../stl_vector_tools.hpp"
#include "../fx/fx_internal.hpp"

using namespace std;

namespace CasADi{

  NonlinearSolve::NonlinearSolve(const std::vector<MX>& x, const ImplicitFunction& implicit_function) : implicit_function_(implicit_function){
    setDependencies(x);
    implicit_function_.init(false);
    setSparsity(implicit_function_.output().sparsity());
  }

  void NonlinearSolve::printPart(std::ostream &stream, int part) const{
    if(part==0){
      stream << implicit_function_ << ".solve(";
    } else if (part == ndep()) {
      stream << ",";
    } else {
      stream << ")";
    }
  }

  void NonlinearSolve::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    casadi_error("not implemented");
  }

  void NonlinearSolve::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
    casadi_error("not implemented");
  }
  
  //  void NonlinearSolve::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
  //    casadi_error("not implemented");
  //  }

  void NonlinearSolve::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    casadi_error("not implemented");
  }

  void NonlinearSolve::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    MXNode::deepCopyMembers(already_copied);
    implicit_function_ = deepcopy(implicit_function_, already_copied);
  }

} // namespace CasADi


