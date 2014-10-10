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


#include "dple_internal.hpp"
#include "dle_internal.hpp"
#include "../std_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../function/mx_function.hpp"
#include "../function/sx_function.hpp"

INPUTSCHEME(DPLEInput)
OUTPUTSCHEME(DPLEOutput)

using namespace std;
namespace casadi {

  DpleInternal::DpleInternal(const DpleStructure & st,
                             int nrhs,
                             bool transp) :
      st_(st), nrhs_(nrhs), transp_(transp) {

    // set default options
    setOption("name", "unnamed_dple_solver"); // name of the function

    addOption("const_dim", OT_BOOLEAN, true, "Assume constant dimension of P");
    addOption("pos_def", OT_BOOLEAN, false, "Assume P positive definite");

    addOption("error_unstable", OT_BOOLEAN, false,
              "Throw an exception when it is detected that Product(A_i, i=N..1) "
              "has eigenvalues greater than 1-eps_unstable");
    addOption("eps_unstable", OT_REAL, 1e-4, "A margin for unstability detection");

    if (nrhs_==1) {
      input_.scheme = SCHEME_DPLEInput;
      output_.scheme = SCHEME_DPLEOutput;
    }

  }

  DpleInternal::~DpleInternal() {

  }

  void DpleInternal::init() {

    const_dim_ = getOption("const_dim");
    pos_def_ = getOption("pos_def");
    error_unstable_ = getOption("error_unstable");
    eps_unstable_ = getOption("eps_unstable");
    A_ = st_[Dple_STRUCT_A];
    V_ = st_[Dple_STRUCT_V];

    // Dimension sanity checks
    casadi_assert_message(A_.size()==V_.size(), "A and V arguments must be of same length, but got "
                          << A_.size() << " and " << V_.size() << ".");
    K_ = A_.size();
    for (int k=0;k<K_;++k) {
      casadi_assert_message(V_[k].isSymmetric(), "V_i must be symmetric but got "
                            << V_[k].dimString() << " for i = " << k << ".");

      casadi_assert_message(A_[k].size1()==V_[k].size1(), "First dimension of A ("
                            << A_[k].size1() << ") must match dimension of symmetric V_i ("
                            << V_[k].size1() << ")" << " for i = " << k << ".");
    }

    if (const_dim_) {
      int n = A_[0].size1();
       for (int k=1;k<K_;++k) {
         casadi_assert_message(A_[k].size1()==n, "You have set const_dim option, but found "
                               "an A_i with dimension ( " << A_[k].dimString()
                               << " ) deviating from n = " << n << " at i = " << k << ".");
      }
    }

    // Allocate inputs
    setNumInputs(1+nrhs_);

    if (const_dim_) {
      input(0)  = DMatrix::zeros(horzcat(A_));
    } else {
      input(0)  = DMatrix::zeros(blkdiag(A_));
    }

    for (int i=0;i<nrhs_;++i) {
      if (const_dim_) {
        input(1+i)  = DMatrix::zeros(horzcat(V_));
      } else {
        input(1+i)  = DMatrix::zeros(blkdiag(V_));
      }
    }

    // Allocate outputs
    std::vector<Sparsity> P = LrDpleInternal::getSparsity(lrdpleStruct("a", A_, "v", V_));

    setNumOutputs(nrhs_);
    for (int i=0;i<nrhs_;++i) {
      if (const_dim_) {
        output(i) = DMatrix::zeros(horzcat(P));
      } else {
        output(i) = DMatrix::zeros(blkdiag(P));
      }
    }

    FunctionInternal::init();

  }

  void DpleInternal::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    FunctionInternal::deepCopyMembers(already_copied);
  }

  std::map<std::string, DpleInternal::Plugin> DpleInternal::solvers_;

  const std::string DpleInternal::infix_ = "dplesolver";

} // namespace casadi


