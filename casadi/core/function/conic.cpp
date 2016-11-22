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


#include "conic_impl.hpp"
#include "nlpsol_impl.hpp"
#include <typeinfo>

using namespace std;
namespace casadi {

  bool has_conic(const string& name) {
    return Conic::has_plugin(name);
  }

  void load_conic(const string& name) {
    Conic::load_plugin(name);
  }

  string doc_conic(const string& name) {
    return Conic::getPlugin(name).doc;
  }

  Function conic(const string& name, const string& solver,
                const SpDict& qp, const Dict& opts) {
    Function ret;
    ret.assignNode(Conic::instantiatePlugin(name, solver, qp));
    ret->construct(opts);
    return ret;
  }

  void Function::conic_debug(const string &filename) const {
    ofstream file;
    file.open(filename.c_str());
    conic_debug(file);
  }

  void Function::conic_debug(ostream &file) const {
    casadi_assert(!is_null());
    const Conic* n = dynamic_cast<const Conic*>(get());
    casadi_assert_message(n!=0, "Not a QP solver");
    return n->generateNativeCode(file);
  }

  vector<string> conic_in() {
    vector<string> ret(conic_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=conic_in(i);
    return ret;
  }

  vector<string> conic_out() {
    vector<string> ret(conic_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=conic_out(i);
    return ret;
  }

  string conic_in(int ind) {
    switch (static_cast<ConicInput>(ind)) {
    case CONIC_H:      return "h";
    case CONIC_G:      return "g";
    case CONIC_A:      return "a";
    case CONIC_LBA:    return "lba";
    case CONIC_UBA:    return "uba";
    case CONIC_LBX:    return "lbx";
    case CONIC_UBX:    return "ubx";
    case CONIC_X0:     return "x0";
    case CONIC_LAM_X0: return "lam_x0";
    case CONIC_LAM_A0: return "lam_a0";
    case CONIC_NUM_IN: break;
    }
    return string();
  }

  string conic_out(int ind) {
    switch (static_cast<ConicOutput>(ind)) {
    case CONIC_X:     return "x";
    case CONIC_COST:  return "cost";
    case CONIC_LAM_A: return "lam_a";
    case CONIC_LAM_X: return "lam_x";
    case CONIC_NUM_OUT: break;
    }
    return string();
  }

  int conic_n_in() {
    return CONIC_NUM_IN;
  }

  int conic_n_out() {
    return CONIC_NUM_OUT;
  }

  template<typename M>
  Function qpsol_nlp(const std::string& name, const std::string& solver,
                     const std::map<std::string, M>& qp, const Dict& opts) {
    // We have: minimize    f(x) = 1/2 * x' H x + c'x
    //          subject to  lbx <= x <= ubx
    //                      lbg <= g(x) = A x + b <= ubg
    M x, p, f, g;
    for (auto&& i : qp) {
      if (i.first=="x") {
        x = i.second;
      } else if (i.first=="p") {
        p = i.second;
      } else if (i.first=="f") {
        f = i.second;
      } else if (i.first=="g") {
        g = i.second;
      } else {
        casadi_error("No such field: " + i.first);
      }
    }
    if (g.is_empty(true)) g = M(0, 1); // workaround

    // Gradient of the objective: gf == Hx + g
    M gf = M::gradient(f, x);

    // Identify the linear term in the objective
    M c = substitute(gf, x, M::zeros(x.sparsity()));

    // Identify the quadratic term in the objective
    M H = M::jacobian(gf, x, {{"symmetric", true}});

    // Identify the constant term in the constraints
    M b = substitute(g, x, M::zeros(x.sparsity()));

    // Identify the linear term in the constraints
    M A = M::jacobian(g, x);

    // Create a function for calculating the required matrices vectors
    Function prob(name + "_qp", {x, p}, {H, c, A, b});

    // Create the QP solver
    Function conic_f = conic(name + "_qpsol", solver,
                             {{"h", H.sparsity()}, {"a", A.sparsity()}}, opts);

    // Create an MXFunction with the right signature
    vector<MX> ret_in(NLPSOL_NUM_IN);
    ret_in[NLPSOL_X0] = MX::sym("x0", x.sparsity());
    ret_in[NLPSOL_P] = MX::sym("p", p.sparsity());
    ret_in[NLPSOL_LBX] = MX::sym("lbx", x.sparsity());
    ret_in[NLPSOL_UBX] = MX::sym("ubx", x.sparsity());
    ret_in[NLPSOL_LBG] = MX::sym("lbg", g.sparsity());
    ret_in[NLPSOL_UBG] = MX::sym("ubg", g.sparsity());
    ret_in[NLPSOL_LAM_X0] = MX::sym("lam_x0", x.sparsity());
    ret_in[NLPSOL_LAM_G0] = MX::sym("lam_g0", g.sparsity());
    vector<MX> ret_out(NLPSOL_NUM_OUT);

    // Get expressions for the QP matrices and vectors
    vector<MX> v(NL_NUM_IN);
    v[NL_X] = ret_in[NLPSOL_X0];
    v[NL_P] = ret_in[NLPSOL_P];
    v = prob(v);

    // Call the QP solver
    vector<MX> w(CONIC_NUM_IN);
    w[CONIC_H] = v.at(0);
    w[CONIC_G] = v.at(1);
    w[CONIC_A] = v.at(2);
    w[CONIC_LBX] = ret_in[NLPSOL_LBX];
    w[CONIC_UBX] = ret_in[NLPSOL_UBX];
    w[CONIC_LBA] = ret_in[NLPSOL_LBG] - v.at(3);
    w[CONIC_UBA] = ret_in[NLPSOL_UBG] - v.at(3);
    w[CONIC_X0] = ret_in[NLPSOL_X0];
    w[CONIC_LAM_X0] = ret_in[NLPSOL_LAM_X0];
    w[CONIC_LAM_A0] = ret_in[NLPSOL_LAM_G0];
    w = conic_f(w);

    // Get expressions for the solution
    ret_out[NLPSOL_X] = w[CONIC_X];
    ret_out[NLPSOL_F] = w[CONIC_COST];
    ret_out[NLPSOL_G] = mtimes(v.at(2), w[CONIC_X]) + v.at(3);
    ret_out[NLPSOL_LAM_X] = w[CONIC_LAM_X];
    ret_out[NLPSOL_LAM_G] = w[CONIC_LAM_A];
    ret_out[NLPSOL_LAM_P] = MX::nan(p.sparsity());
    return Function(name, ret_in, ret_out, nlpsol_in(), nlpsol_out(),
                    {{"default_in", nlpsol_default_in()}});
  }

  Function qpsol(const std::string& name, const std::string& solver,
                 const SXDict& qp, const Dict& opts) {
    return qpsol_nlp(name, solver, qp, opts);
  }

  Function qpsol(const std::string& name, const std::string& solver,
                 const MXDict& qp, const Dict& opts) {
    return qpsol_nlp(name, solver, qp, opts);
  }

  // Constructor
  Conic::Conic(const std::string& name, const std::map<std::string, Sparsity> &st)
    : FunctionInternal(name) {
    for (auto i=st.begin(); i!=st.end(); ++i) {
      if (i->first=="a") {
        A_ = i->second;
      } else if (i->first=="h") {
        H_ = i->second;
      } else {
        casadi_error("Unrecognized field in QP structure: " << i->first);
      }
    }

    // We need either A or H
    casadi_assert_message(!A_.is_null() || !H_.is_null(),
      "Cannot determine dimension");

    // Generate A or H
    if (A_.is_null()) {
      A_ = Sparsity(0, H_.size2());
    } else if (H_.is_null()) {
      H_ = Sparsity(A_.size2(), A_.size2());
    } else {
      // Consistency check
      casadi_assert_message(A_.size2()==H_.size2(),
        "Got incompatible dimensions.   min          x'Hx + G'x s.t.   LBA <= Ax <= UBA :"
        << std::endl <<
        "H: " << H_.dim() << " - A: " << A_.dim() << std::endl <<
        "We need: H_.size2()==A_.size2()" << std::endl);
    }

    casadi_assert_message(H_.is_symmetric(),
      "Got incompatible dimensions.   min          x'Hx + G'x" << std::endl <<
      "H: " << H_.dim() <<
      "We need H square & symmetric" << std::endl);

    nx_ = A_.size2();
    na_ = A_.size1();
  }

  Sparsity Conic::get_sparsity_in(int i) {
    switch (static_cast<ConicInput>(i)) {
    case CONIC_X0:
    case CONIC_G:
    case CONIC_LBX:
    case CONIC_UBX:
    case CONIC_LAM_X0:
      return get_sparsity_out(CONIC_X);
    case CONIC_LBA:
    case CONIC_UBA:
    case CONIC_LAM_A0:
      return get_sparsity_out(CONIC_LAM_A);
    case CONIC_A:
      return A_;
    case CONIC_H:
      return H_;
    case CONIC_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity Conic::get_sparsity_out(int i) {
    switch (static_cast<ConicOutput>(i)) {
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

  Options Conic::options_
  = {{&FunctionInternal::options_},
     {{"discrete",
       {OT_BOOLVECTOR,
        "Indicates which of the variables are discrete, i.e. integer-valued"}}
     }
  };

  void Conic::init(const Dict& opts) {
    // Call the init method of the base class
    FunctionInternal::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="discrete") {
        discrete_ = op.second;
      }
    }

    // Check options
    if (!discrete_.empty()) {
      casadi_assert_message(discrete_.size()==nx_, "\"discrete\" option has wrong length");
      if (std::find(discrete_.begin(), discrete_.end(), true)!=discrete_.end()) {
        casadi_assert_message(integer_support(),
                              "Discrete variables require a solver with integer support");
      }
    }
  }

  Conic::~Conic() {
  }

  void Conic::checkInputs(const double* lbx, const double* ubx,
                          const double* lba, const double* uba) const {
    for (int i=0; i<nx_; ++i) {
      double lb = lbx ? lbx[i] : 0., ub = ubx ? ubx[i] : 0.;
      casadi_assert_message(lb <= ub,
                            "LBX[" << i << "] <= UBX[" << i << "] was violated. "
                            << "Got LBX[" << i << "]=" << lb <<
                            " and UBX[" << i << "] = " << ub << ".");
    }
    for (int i=0; i<na_; ++i) {
      double lb = lba ? lba[i] : 0., ub = uba ? uba[i] : 0.;
      casadi_assert_message(lb <= ub,
                            "LBA[" << i << "] <= UBA[" << i << "] was violated. "
                            << "Got LBA[" << i << "] = " << lb <<
                            " and UBA[" << i << "] = " << ub << ".");
    }
  }

  void Conic::generateNativeCode(std::ostream& file) const {
    casadi_error("Conic::generateNativeCode not defined for class "
                 << typeid(*this).name());
  }

  std::map<std::string, Conic::Plugin> Conic::solvers_;

  const std::string Conic::infix_ = "conic";

  double Conic::default_in(int ind) const {
    switch (ind) {
    case CONIC_LBX:
    case CONIC_LBA:
      return -std::numeric_limits<double>::infinity();
    case CONIC_UBX:
    case CONIC_UBA:
      return std::numeric_limits<double>::infinity();
    default:
      return 0;
    }
  }

  void Conic::print_fstats(const ConicMemory* m) const {

    size_t maxNameLen=0;

    // Retrieve all qp keys
    std::vector<std::string> keys;
    for (auto &&s : m->fstats) {
      maxNameLen = max(s.first.size(), maxNameLen);
      keys.push_back(s.first);
    }

    // Print header
    std::stringstream s;
    std::string blankName(maxNameLen, ' ');
    s
      << blankName
      << "      proc           wall      num           mean             mean"
      << endl << blankName
      << "      time           time     evals       proc time        wall time";
    userOut() << s.str() << endl;

    std::sort(keys.begin(), keys.end());
    for (auto k : keys) {
      const FStats& fs = m->fstats.at(k);
      print_stats_line(maxNameLen, k, fs.n_call, fs.t_proc, fs.t_wall);
    }

    // Sum the previously printed stats
    double t_wall_all_previous = 0;
    double t_proc_all_previous = 0;
    for (auto k : keys) {
      const FStats& fs = m->fstats.at(k);
      t_proc_all_previous += fs.t_proc;
      t_wall_all_previous += fs.t_wall;
    }
    print_stats_line(maxNameLen, "all previous", -1, t_proc_all_previous, t_wall_all_previous);

  }

} // namespace casadi
