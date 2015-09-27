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


#include "lr_dle_internal.hpp"
#include "../std_vector_tools.hpp"
#include "../function/mx_function.hpp"
#include "../function/sx_function.hpp"

INPUTSCHEME(LR_DLEInput)
OUTPUTSCHEME(LR_DLEOutput)

using namespace std;
namespace casadi {

  LrDleInternal::LrDleInternal(const std::map<std::string, Sparsity>& st,
                             int nrhs,
                             bool transp) :
    nrhs_(nrhs), transp_(transp) {

    // set default options
    setOption("name", "unnamed_dple_solver"); // name of the function

    addOption("pos_def", OT_BOOLEAN, false, "Assume P positive definite");
    addOption("Hs",      OT_INTEGERVECTOR, std::vector<int>(), "Column-sizes of H_i");
    addOption("error_unstable", OT_BOOLEAN, false,
              "Throw an exception when it is detected that Product(A_i, i=N..1) "
              "has eigenvalues greater than 1-eps_unstable");
    addOption("eps_unstable", OT_REAL, 1e-4, "A margin for unstability detection");

    st_.resize(LR_DLE_STRUCT_NUM);
    for (std::map<std::string, Sparsity>::const_iterator i=st.begin(); i!=st.end(); ++i) {
      if (i->first=="a") {
        st_[LR_DLE_STRUCT_A]=i->second;
      } else if (i->first=="v") {
        st_[LR_DLE_STRUCT_V]=i->second;
      } else if (i->first=="c") {
        st_[LR_DLE_STRUCT_C]=i->second;
      } else if (i->first=="h") {
        st_[LR_DLE_STRUCT_H]=i->second;
      } else {
        casadi_error("Unrecognized field in LR DLE structure: " << i->first);
      }
    }


    if (nrhs_==1) {
      ischeme_ = IOScheme(SCHEME_LR_DLEInput);
      oscheme_ = IOScheme(SCHEME_LR_DLEOutput);
    }
  }

  LrDleInternal::~LrDleInternal() {

  }

  void LrDleInternal::init() {

    pos_def_ = getOption("pos_def");
    error_unstable_ = getOption("error_unstable");
    eps_unstable_ = getOption("eps_unstable");
    Hs_ = getOption("Hs");

    A_ = st_[LR_DLE_STRUCT_A];
    V_ = st_[LR_DLE_STRUCT_V];
    C_ = st_[LR_DLE_STRUCT_C];
    H_ = st_[LR_DLE_STRUCT_H];

    int n = A_.size1();

    with_H_ = true;

    // Default H: unity
    if (H_.isNull() || H_.isempty()) {
      //H_ = Sparsity::diag(n);
      H_ = Sparsity(0, 0);
      //casadi_assert(Hs_.size()==0);
      with_H_ = false;
    }

    if (with_H_) {
      casadi_assert_message(H_.size1()==n, "Number of rows in H must be n (" << n << "), but got "
                          << H_.size1() << ".");

      // Default hs: [H.size2()]
      if (Hs_.size()==0) {
        Hs_.push_back(H_.size2());
      }

      // Assert that sum of Hs entries match up to H.size2()
      double sum = 0;
      for (int k=0;k<Hs_.size();++k) {
        sum += Hs_[k];
      }

      casadi_assert_message(H_.size2()==sum, "Number of columns in H (" << H_.size2() << "),"
                                             << "must match sum(Hs): " << sum << ".");

    }

    with_C_ = true;
    if (C_.isNull()  || C_.isempty()) {
      C_ = Sparsity(0, 0);
      //C_ = Sparsity::diag(n);
      //st_[Dle_STRUCT_C] = C_;
      with_C_ = false;
    }



    casadi_assert_message(V_.issymmetric(), "V must be symmetric but got "
                          << V_.dimString() << ".");

    casadi_assert_message(A_.size1()==A_.size2(), "A must be square but got "
                          << A_.dimString() << ".");

    int m = with_C_? V_.size1(): n;

    if (with_C_) {
      casadi_assert_message(n==C_.size1(), "Number of rows in C ("
                            << C_.size1() << ") must match dimension of square A ("
                            << n << ")" << ".");
      casadi_assert_message(m==C_.size2(), "Number of columns in C ("
                            << C_.size2() << ") must match dimension of symmetric V ("
                            << m << ")" << ".");
    } else {
      casadi_assert_message(A_.size1()==V_.size1(), "First dimension of A ("
                            << A_.size1() << ") must match dimension of symmetric V ("
                            << V_.size1() << ")" << ".");
    }

    casadi_assert(nrhs_==1);

    // Allocate inputs
    ibuf_.resize(LR_DLE_NUM_IN);
    input(LR_DLE_A)  = DMatrix::zeros(A_);


    if (with_C_) {
      input(LR_DLE_C)  = DMatrix::zeros(C_);
    }
    input(LR_DLE_V)  = DMatrix::zeros(V_);
    if (with_H_) {
      input(LR_DLE_H)  = DMatrix::zeros(H_);
    }

    obuf_.resize(nrhs_*LR_DLE_NUM_OUT);

    Sparsity Pnew = getSparsity(make_map("a", st_[LR_DLE_STRUCT_A],
                                         "v", st_[LR_DLE_STRUCT_V],
                                         "c", st_[LR_DLE_STRUCT_C],
                                         "h", st_[LR_DLE_STRUCT_H]),
                                Hs_);

    for (int i=0;i<nrhs_;++i) {
      output(i) = DMatrix::zeros(Pnew);
    }

    if (with_H_) {
      Hi_ = std::vector<int>(Hs_.size()+1, 0);
      for (int k=0;k<Hs_.size();++k) {
        Hi_[k+1] = Hi_[k] + Hs_[k];
      }
      Hv_ = horzsplit(input(LR_DLE_H), Hi_);

      Pv_.resize(Hs_.size());
      Pi_.resize(Hs_.size()+1, 0);
      std::vector<Sparsity> sp = diagsplit(Pnew, Hi_);
      for (int k=0;k<Hs_.size();++k) {
        Pv_[k] = DMatrix::zeros(sp[k]);
        Pi_[k+1] = Pi_[k] + sp[k].nnz();
      }
    }

    FunctionInternal::init();

  }

  Sparsity LrDleInternal::getSparsity(const std::map<std::string, Sparsity>& st,
                                      const std::vector<int> &Hs) {

    // Compute output sparsity by Smith iteration with frequency doubling
    Sparsity A, V, C, H;
    if (st.count("a")) A = st.at("a");

    int n = A.size1();
    if (st.count("v")) V = st.at("v");
    if (st.count("c")) C = st.at("c");

    if (C.isNull() || C.isempty()) C = Sparsity::diag(n);

    if (st.count("h")) H = st.at("h");
    bool with_H = !(H.isNull() || H.isempty());
    if (H.isNull() || H.isempty()) {
      H = Sparsity::diag(n);

    }

    Sparsity P = mul(mul(C, V), C.T());
    Sparsity Pprev(n, n);

    while (Pprev.nnz()!=P.nnz()) {
      // This can be much improved:
      //   * No need to make C and V grow
      //   * norm_0 instead of constructing P
      Pprev = P;
      C = horzcat(C, mul(A, C));
      V = diagcat(V, V);
      A = mul(A, A);
      P = mul(mul(C, V), C.T());
    }

    if (with_H) {
      int hn = Hs.size();

      std::vector<Sparsity> sp(hn);

      std::vector<int> Hi_(hn+1, 0);
      for (int k=0;k<hn;++k) {
        Hi_[k+1] = Hi_[k] + Hs[k];
      }
      std::vector<Sparsity> Hv_ = horzsplit(H, Hi_);

      for (int k=0;k<hn;++k) {
        sp[k] = mul(Hv_[k].T(), mul(P, Hv_[k]));
      }

      return diagcat(sp);
    } else {
      return P;
    }
  }

  std::map<std::string, LrDleInternal::Plugin> LrDleInternal::solvers_;

  const std::string LrDleInternal::infix_ = "lrdlesolver";

} // namespace casadi


