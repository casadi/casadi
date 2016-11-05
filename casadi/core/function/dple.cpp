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

  MX dplesol(const MX& A, const MX& V, const std::string& solver, const Dict& opts) {
    SpDict sp;
    sp["a"] = A.sparsity();
    sp["v"] = V.sparsity();
    Function f = dplesol("dplesol", solver, sp, opts);
    MXDict f_in;
    f_in["a"] = A;
    f_in["v"] = V;
    MXDict f_out = f(f_in);
    return f_out["p"];
  }
  Function dplesol(const string& name, const string& solver,
                const SpDict& st, const Dict& opts) {
    Function ret;
    ret.assignNode(Dple::instantiatePlugin(name, solver, st));
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
  Dple::Dple(const std::string& name, const SpDict &st)
    : FunctionInternal(name) {
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
      case DPLE_A:
        return A_;
      case DPLE_V:
        return V_;
      case DPLE_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity Dple::get_sparsity_out(int i) {
    switch (static_cast<DpleOutput>(i)) {
      case DPLE_P:
        return V_;
      case DPLE_NUM_OUT: break;
    }
    return Sparsity();
  }

  Options Dple::options_
  = {{&FunctionInternal::options_},
     {{"const_dim",
       {OT_BOOL,
        "Assume constant dimension of P"}},
      {"pos_def",
        {OT_BOOL,
         "Assume P positive definite"}},
      {"error_unstable",
        {OT_BOOL,
        "Throw an exception when it is detected that Product(A_i, i=N..1)"
        "has eigenvalues greater than 1-eps_unstable"}},
      {"eps_unstable",
        {OT_DOUBLE,
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
    eps_unstable_ = 1e-4;

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

    casadi_assert(V_.size2() % V_.size1() == 0);
    nrhs_ = V_.size2() / V_.size1();
    casadi_assert(nrhs_>=1);

    std::vector<Sparsity> Vs = horzsplit(V_, V_.size1());
    Sparsity Vref = Vs[0];
    casadi_assert_message(Vref.is_symmetric(), "V must be symmetric but got " << Vref.dim() << ".");

    for (auto&& s : Vs)
      casadi_assert(s==Vref);

    casadi_assert_message(const_dim_, "Not implemented");

    int blocksize = Vref.colind()[1];
    int N = Vref.size1()/blocksize;
    Sparsity block = Sparsity::dense(blocksize, blocksize);

    std::vector<Sparsity> blocks(N, block);
    casadi_assert_message(Vref==diagcat(blocks), "Structure not recognised.");
    casadi_assert_message(A_==Vref, "Structure not recognised.");


  }

  Dple::~Dple() {
  }

  std::map<std::string, Dple::Plugin> Dple::solvers_;

  const std::string Dple::infix_ = "dple";

  double Dple::default_in(int ind) const {
    return 0;
  }

} // namespace casadi
