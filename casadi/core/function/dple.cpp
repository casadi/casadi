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


#include "dple_impl.hpp"
#include <typeinfo>

using namespace std;
namespace casadi {

  bool has_dple(const string& name) {
    return Dple::has_plugin(name);
  }

  void load_dple(const string& name) {
    Dple::load_plugin(name);
  }

  string doc_dple(const string& name) {
    return Dple::getPlugin(name).doc;
  }

  Function dple(const string& name, const string& solver,
                const SpDict& qp, const Dict& opts) {
    Function ret;
    ret.assignNode(Dple::instantiatePlugin(name, solver, qp));
    ret->construct(opts);
    return ret;
  }

  vector<string> dple_in() {
    vector<string> ret(dple_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=dple_in(i);
    return ret;
  }

  vector<string> dple_out() {
    vector<string> ret(dple_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=dple_out(i);
    return ret;
  }

  string dple_in(int ind) {
    switch (static_cast<DpleInput>(ind)) {
    case DPLE_A:      return "a";
    case DPLE_V:      return "v";
    case DPLE_NUM_IN: break;
    }
    return string();
  }

  string dple_out(int ind) {
    switch (static_cast<DpleOutput>(ind)) {
      case DPLE_P:      return "p";
      case DPLE_NUM_OUT: break;
    }
    return string();
  }

  int dple_n_in() {
    return DPLE_NUM_IN;
  }

  int dple_n_out() {
    return DPLE_NUM_OUT;
  }

  // Constructor
  Dple::Dple(const std::string& name, const std::map<std::string, Sparsity> &st, int nrhs, bool transp)
    : FunctionInternal(name), nrhs_(nrhs), trans_(transp) {
    for (auto i=st.begin(); i!=st.end(); ++i) {
      if (i->first=="a") {
        A_ = i->second;
      } else if (i->first=="v") {
        V_ = i->second;
      } else {
        casadi_error("Unrecognized field in Dple structure: " << i->first);
      }
    }

  }

  Sparsity Dple::get_sparsity_in(int i) {
    switch (static_cast<DpleInput>(i)) {
    case CONIC_X0:
    case CONIC_G:
    case CONIC_LBX:
    case CONIC_UBX:
    case CONIC_LAM_X0:
      return get_sparsity_out(CONIC_X);
    case CONIC_LBA:
    case CONIC_UBA:
      return get_sparsity_out(CONIC_LAM_A);
    case CONIC_A:
      return A_;
    case CONIC_H:
      return H_;
    case CONIC_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity Dple::get_sparsity_out(int i) {
    switch (static_cast<DpleOutput>(i)) {
    case CONIC_COST:
      return Sparsity::scalar();
    case CONIC_X:
    case CONIC_LAM_X:
      return Sparsity::dense(nx_, 1);
    case CONIC_LAM_A:
      return Sparsity::dense(na_, 1);
    case CONIC_NUM_OUT: break;
    }
    return Sparsity();
  }

  Options Dple::options_
  = {{&FunctionInternal::options_},
     {{"const_dim",
       {OT_BOOLEAN,
        "Assume constant dimension of P"}},
      {"pos_def",
        {OT_BOOLEAN,
         "Assume P positive definite"}},
      {"error_unstable",
        {OT_BOOLEAN,
        "Throw an exception when it is detected that Product(A_i, i=N..1)"
        "has eigenvalues greater than 1-eps_unstable"}},
      {"eps_unstable",
        {OT_REAL,
        "A margin for unstability detection"}}
     }
  };

  void Dple::init(const Dict& opts) {
    // Call the init method of the base class
    FunctionInternal::init(opts);

    // Default options
    const_dim_ = true;
    pos_def_ = false;
    error_unstable_ = false;
    eps_unstable = 1e-4;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="const_dim") {
        const_dim_ = op.second;
      } else if  (op.first=="pos_def") {
        pos_def_ = op.second;
      } else if  (op.first=="error_unstable") {
        error_unstable_ = op.second;
      } else if  (op.first=="eps_unstable") {
        eps_unstable_ = op.second;
      }
    }


    // Dimension sanity checks
    casadi_assert_message(A_.size()==V_.size(), "A and V arguments must be of same length, but got "
                          << A_.size() << " and " << V_.size() << ".");
    K_ = A_.size();
    for (int k=0;k<K_;++k) {
      casadi_assert_message(V_[k].issymmetric(), "V_i must be symmetric but got "
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
  }

  Dple::~Dple() {
  }

  std::map<std::string, Dple::Plugin> Dple::solvers_;

  const std::string Dple::infix_ = "dple";

  double Dple::default_in(int ind) const {
    return 0;
  }

} // namespace casadi
