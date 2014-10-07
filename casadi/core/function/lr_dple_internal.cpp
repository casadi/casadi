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


#include "lr_dple_internal.hpp"
#include "lr_dle_internal.hpp"
#include "../std_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../function/mx_function.hpp"
#include "../function/sx_function.hpp"

INPUTSCHEME(LR_DPLEInput)
OUTPUTSCHEME(LR_DPLEOutput)

using namespace std;
namespace casadi {

  LrDpleInternal::LrDpleInternal(const LrDpleStructure & st,
                             const std::vector< std::vector<int> > &Hs,
                             int nrhs,
                             bool transp) :
      st_(st), Hs_(Hs), nrhs_(nrhs), transp_(transp) {

    // set default options
    setOption("name", "unnamed_dple_solver"); // name of the function

    addOption("const_dim", OT_BOOLEAN, true, "Assume constant dimension of P");
    addOption("pos_def", OT_BOOLEAN, false, "Assume P positive definite");

    addOption("error_unstable", OT_BOOLEAN, false,
              "Throw an exception when it is detected that Product(A_i, i=N..1) "
              "has eigenvalues greater than 1-eps_unstable");
    addOption("eps_unstable", OT_REAL, 1e-4, "A margin for unstability detection");

    if (nrhs_==1) {
      input_.scheme = SCHEME_LR_DPLEInput;
      output_.scheme = SCHEME_LR_DPLEOutput;
    }

  }

  LrDpleInternal::~LrDpleInternal() {

  }

  void LrDpleInternal::init() {

    const_dim_ = getOption("const_dim");
    pos_def_ = getOption("pos_def");
    error_unstable_ = getOption("error_unstable");
    eps_unstable_ = getOption("eps_unstable");

    casadi_assert(const_dim_);

    A_ = st_[LR_Dple_STRUCT_A];
    casadi_assert(A_.size()>0);
    int n = A_[0].size1();

    V_ = st_[LR_Dple_STRUCT_V];
    C_ = st_[LR_Dple_STRUCT_C];
    H_ = st_[LR_Dple_STRUCT_H];

    with_H_ = true;

    if (H_.empty()) {
      with_H_ = false;
    }



    with_C_ = true;
    if (C_.empty()) {
      with_C_ = false;
    }

    if (with_H_) {
      casadi_assert_message(A_.size()==H_.size(),
                          "A and H arguments must be of same length, but got "
                          << A_.size() << " and " << H_.size() << ".");
      casadi_assert_message(H_.size()==Hs_.size(),
                          "H and Hs arguments must be of same length, but got "
                          << H_.size() << " and " << Hs_.size() << ".");

      Hsi_.resize(1, 0);
      for (int k=0;k<H_.size();++k) {

        // Default hs: [H.size2()]
        if (Hs_[k].size()==0) {
          Hs_[k].push_back(H_[k].size2());
        }

        // Assert that sum of Hs entries match up to H.size2()
        double sum = 0;
        for (int i=0;i<Hs_[k].size();++i) {
          sum += Hs_[k][i];
        }

        Hsi_.push_back(Hsi_.back()+sum);

        casadi_assert_message(H_[k].size2()==sum,
                       "Number of columns in H @ k = " << k << "  (" << H_[k].size2() << "),"
                       << "must match sum(Hs[k]): " << sum << ".");
      }

    }


    // Dimension sanity checks
    casadi_assert_message(A_.size()==V_.size(), "A and V arguments must be of same length, but got "
                          << A_.size() << " and " << V_.size() << ".");
    K_ = A_.size();

    int m = with_C_? V_[0].size1(): n;

    for (int k=0;k<K_;++k) {
      casadi_assert_message(V_[k].isSymmetric(), "V_i must be symmetric but got "
                            << V_[k].dimString() << " for i = " << k << ".");
      if (with_C_) {
        casadi_assert_message(n==C_[k].size1(), "Number of rows in C ("
                              << C_[k].size1() << ") must match dimension of square A @ k = "
                              << k << " (" << n << ")" << ".");
        casadi_assert_message(m==C_[k].size2(), "Number of columns in C ("
                              << C_[k].size2() << ") must match dimension of symmetric V @ k = "
                              << k << " (" << m << ")" << ".");
      } else {
          casadi_assert_message(A_[k].size1()==V_[k].size1(), "First dimension of A ("
                            << A_[k].size1() << ") must match dimension of symmetric V @ k = "
                            << k << " (" << V_[k].size1() << ")" << ".");
      }
    }



    if (const_dim_) {
      int n = A_[0].size1();
       for (int k=1;k<K_;++k) {
         casadi_assert_message(A_[k].size1()==n, "You have set const_dim option, but found "
                               "an A_i with dimension ( " << A_[k].dimString()
                               << " ) deviating from n = " << n << " at i = " << k << ".");
      }
    }

    Hss_.resize(0);

    if (Hs_.size()>0) {
      for (int k=0;k<K_;++k) {
        Hss_.insert(Hss_.end(), Hs_[k].begin(), Hs_[k].end());
      }
    }

    // Allocate inputs
    setNumInputs(LR_DPLE_NUM_IN);

    if (const_dim_) {
      input(LR_DPLE_A)  = DMatrix::zeros(horzcat(A_));

      if (with_C_) {
        input(LR_DPLE_C)  = DMatrix::zeros(horzcat(C_));
      }
      input(LR_DPLE_V)  = DMatrix::zeros(horzcat(V_));
      if (with_H_) {
        input(LR_DPLE_H)  = DMatrix::zeros(horzcat(H_));
      }
    } else {
      input(LR_DPLE_A)  = DMatrix::zeros(blkdiag(A_));
    }

    /**for (int i=0;i<nrhs_;++i) {
      if (const_dim_) {
        input(1+i)  = DMatrix::zeros(horzcat(V_));
      } else {
        input(1+i)  = DMatrix::zeros(blkdiag(V_));
      }
    }*/

    // Allocate outputs
    std::vector<Sparsity> P = getSparsity(st_, Hs_);

    setNumOutputs(nrhs_*LR_DPLE_NUM_OUT);
    for (int i=0;i<nrhs_;++i) {
      if (const_dim_) {
        output(i) = DMatrix::zeros(horzcat(P));
      } else {
        output(i) = DMatrix::zeros(blkdiag(P));
      }
    }

    FunctionInternal::init();

  }

  void LrDpleInternal::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    FunctionInternal::deepCopyMembers(already_copied);
  }

  std::map<std::string, LrDpleInternal::Plugin> LrDpleInternal::solvers_;

  const std::string LrDpleInternal::infix_ = "lrdplesolver";


  std::vector<Sparsity> LrDpleInternal::getSparsity(
    const LrDpleStructure& st, const std::vector< std::vector<int> > &Hs_) {

    // Chop-up the arguments
    std::vector<Sparsity> As = st[LR_Dple_STRUCT_A];
    std::vector<Sparsity> Vs = st[LR_Dple_STRUCT_V];
    std::vector<Sparsity> Cs = st[LR_Dple_STRUCT_C];
    std::vector<Sparsity> Hs = st[LR_Dple_STRUCT_H];

    bool with_H = !st[LR_Dple_STRUCT_H].empty();

    int K = As.size();

    Sparsity A;
    if (K==1) {
      A = As[0];
    } else {
      Sparsity AL = blkdiag(vector_slice(As, range(As.size()-1)));

      Sparsity AL2 = horzcat(AL, Sparsity::sparse(AL.size1(), As[0].size2()));
      Sparsity AT = horzcat(Sparsity::sparse(As[0].size1(), AL.size2()), As.back());
      A = vertcat(AT, AL2);
    }

    Sparsity V = blkdiag(Vs.back(), blkdiag(vector_slice(Vs, range(Vs.size()-1))));
    Sparsity C;

    if (!st[LR_Dple_STRUCT_C].empty()) {
      C = blkdiag(Cs.back(), blkdiag(vector_slice(Cs, range(Cs.size()-1))));
    }
    Sparsity H;
    std::vector<int> Hs_agg;

    std::vector<int> Hi(1, 0);
    if (Hs_.size()>0 && with_H) {
      H = blkdiag(Hs.back(), blkdiag(vector_slice(Hs, range(Hs.size()-1))));
      casadi_assert(K==Hs_.size());
      for (int k=0;k<K;++k) {
        Hs_agg.insert(Hs_agg.end(), Hs_[k].begin(), Hs_[k].end());
        int sum=0;
        for (int i=0;i<Hs_[k].size();++i) {
          sum+= Hs_[k][i];
        }
        Hi.push_back(Hi.back()+sum);
      }
    }

    Sparsity res = LrDleInternal::getSparsity(lrdleStruct("a", A, "v", V, "c", C, "h", H), Hs_agg);

    if (with_H) {
      return diagsplit(res, Hi);
    } else {
      return diagsplit(res, As[0].size2());
    }
  }


} // namespace casadi


