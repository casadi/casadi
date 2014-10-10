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


#include "cle_internal.hpp"
#include "../std_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../function/mx_function.hpp"
#include "../function/sx_function.hpp"

INPUTSCHEME(CLEInput)
OUTPUTSCHEME(CLEOutput)

using namespace std;
namespace casadi {

  CleInternal::CleInternal(const CleStructure& st,
                             int nrhs,
                             bool transp) :
      st_(st), nrhs_(nrhs), transp_(transp) {

    // set default options
    setOption("name", "unnamed_dple_solver"); // name of the function

    addOption("pos_def", OT_BOOLEAN, false, "Assume P positive definite");

    addOption("error_unstable", OT_BOOLEAN, false,
              "Throw an exception when it is detected that Product(A_i, i=N..1) "
              "has eigenvalues greater than 1-eps_unstable");
    addOption("eps_unstable", OT_REAL, 1e-4, "A margin for unstability detection");

    if (nrhs_==1) {
      input_.scheme = SCHEME_CLEInput;
      output_.scheme = SCHEME_CLEOutput;
    }

  }

  CleInternal::~CleInternal() {

  }

  void CleInternal::init() {

    pos_def_ = getOption("pos_def");
    error_unstable_ = getOption("error_unstable");
    eps_unstable_ = getOption("eps_unstable");

    A_ = st_[Cle_STRUCT_A];
    V_ = st_[Cle_STRUCT_V];

    casadi_assert_message(V_.isSymmetric(), "V must be symmetric but got "
                          << V_.dimString() << ".");

    casadi_assert_message(A_.size1()==V_.size1(), "First dimension of A ("
                          << A_.size1() << ") must match dimension of symmetric V ("
                          << V_.size1() << ")" << ".");

    int n = A_.size1();
    casadi_assert_message(A_.size1()==n, "You have set const_dim option, but found "
                           "an A with dimension ( " << A_.dimString()
                           << " ) deviating from n = " << n << ".");

    // Allocate inputs
    setNumInputs(1+nrhs_);

    input(0)  = DMatrix::zeros(A_);

    for (int i=0;i<nrhs_;++i) {
      input(1+i)  = DMatrix::zeros(V_);
    }

    // Allocate outputs
    Sparsity P = Sparsity::dense(V_.size1(), V_.size1());

    setNumOutputs(nrhs_);
    for (int i=0;i<nrhs_;++i) {
      output(i) = DMatrix::zeros(P);
    }

    FunctionInternal::init();

  }

  void CleInternal::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    FunctionInternal::deepCopyMembers(already_copied);
  }

  std::map<std::string, CleInternal::Plugin> CleInternal::solvers_;

  const std::string CleInternal::infix_ = "clesolver";


} // namespace casadi


