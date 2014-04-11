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

#include "matrix_tools.hpp"
#include "../std_vector_tools.hpp"
#include "../function/linear_solver.hpp"

using namespace std;

namespace casadi{

  Matrix<double> solve(const Matrix<double>& A, const Matrix<double>& b, linearSolverCreator lsolver, const Dictionary& dict) {
    LinearSolver mysolver = lsolver(A.sparsity(),b.size2());
    mysolver.setOption(dict);
    mysolver.init();
    mysolver.setInput(A,LINSOL_A);
    mysolver.setInput(b,LINSOL_B);
    mysolver.prepare();
    mysolver.solve(false);

    return mysolver.output(LINSOL_X);
  }


  Matrix<double> pinv(const Matrix<double>& A, linearSolverCreator lsolver, const Dictionary& dict) {
    if (A.size1()>=A.size2()) {
      return solve(mul(A.T(),A),A.T(),lsolver,dict);
    } else {
      return solve(mul(A,A.T()),A,lsolver,dict).T();
    }
  }

} // namespace casadi


