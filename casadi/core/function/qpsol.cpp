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


#include "qpsol_impl.hpp"
#include "nlpsol_impl.hpp"
#include <typeinfo>

using namespace std;
namespace casadi {

  bool has_qpsol(const string& name) {
    return Qpsol::hasPlugin(name);
  }

  void load_qpsol(const string& name) {
    Qpsol::loadPlugin(name);
  }

  string doc_qpsol(const string& name) {
    return Qpsol::getPlugin(name).doc;
  }

  Function qpsol(const string& name, const string& solver,
                               const SpDict& qp, const Dict& opts) {
    Function ret;
    ret.assignNode(Qpsol::instantiatePlugin(name, solver, qp));
    ret->construct(opts);
    return ret;
  }

  void Function::qpsol_debug(const string &filename) const {
    ofstream file;
    file.open(filename.c_str());
    qpsol_debug(file);
  }

  void Function::qpsol_debug(ostream &file) const {
    casadi_assert(!is_null());
    const Qpsol* n = dynamic_cast<const Qpsol*>(get());
    casadi_assert_message(n!=0, "Not a QP solver");
    return n->generateNativeCode(file);
  }

  vector<string> qpsol_in() {
    vector<string> ret(qpsol_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=qpsol_in(i);
    return ret;
  }

  vector<string> qpsol_out() {
    vector<string> ret(qpsol_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=qpsol_out(i);
    return ret;
  }

  string qpsol_in(int ind) {
    switch (static_cast<QpsolInput>(ind)) {
    case QPSOL_H:      return "h";
    case QPSOL_G:      return "g";
    case QPSOL_A:      return "a";
    case QPSOL_LBA:    return "lba";
    case QPSOL_UBA:    return "uba";
    case QPSOL_LBX:    return "lbx";
    case QPSOL_UBX:    return "ubx";
    case QPSOL_X0:     return "x0";
    case QPSOL_LAM_X0: return "lam_x0";
    case QPSOL_NUM_IN: break;
    }
    return string();
  }

  string qpsol_out(int ind) {
    switch (static_cast<QpsolOutput>(ind)) {
    case QPSOL_X:     return "x";
    case QPSOL_COST:  return "cost";
    case QPSOL_LAM_A: return "lam_a";
    case QPSOL_LAM_X: return "lam_x";
    case QPSOL_NUM_OUT: break;
    }
    return string();
  }

  int qpsol_n_in() {
    return QPSOL_NUM_IN;
  }

  int qpsol_n_out() {
    return QPSOL_NUM_OUT;
  }

  template<typename M>
  Function qpsol_nlp(const std::string& name, const std::string& solver,
                     const Problem<M>& qp, const Dict& opts) {
    // We have: minimize    f(x) = 1/2 * x' H x + c'x
    //          subject to  lbx <= x <= ubx
    //                      lbg <= g(x) = A x + b <= ubg
    M x = qp.in[NL_X];
    M p = qp.in[NL_P];
    M f = qp.out[NL_F];
    M g = qp.out[NL_G];
    if (g.is_empty(true)) g = M(0, 1); // workaround

    // Gradient of the objective: gf == Hx + g
    M gf = M::gradient(f, x);

    // Identify the linear term in the objective
    M c = substitute(gf, x, M::zeros(x.sparsity()));

    // Identify the quadratic term in the objective
    M H = M::jacobian(gf, x, true);

    // Identify the constant term in the constraints
    M b = substitute(g, x, M::zeros(x.sparsity()));

    // Identify the linear term in the constraints
    M A = M::jacobian(g, x);

    // Create a function for calculating the required matrices vectors
    Function prob(name + "_qp", qp.in, {H, c, A, b});

    // Create the QP solver
    Function qpsol_f = qpsol(name + "_qpsol", solver,
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
    vector<MX> w(QPSOL_NUM_IN);
    w[QPSOL_H] = v.at(0);
    w[QPSOL_G] = v.at(1);
    w[QPSOL_A] = v.at(2);
    w[QPSOL_LBX] = ret_in[NLPSOL_LBX];
    w[QPSOL_UBX] = ret_in[NLPSOL_UBX];
    w[QPSOL_LBA] = ret_in[NLPSOL_LBG] - v.at(3);
    w[QPSOL_UBA] = ret_in[NLPSOL_UBG] - v.at(3);
    w[QPSOL_X0] = ret_in[NLPSOL_X0];
    w[QPSOL_LAM_X0] = ret_in[NLPSOL_LAM_X0];
    w = qpsol_f(w);

    // Get expressions for the solution
    ret_out[NLPSOL_X] = w[QPSOL_X];
    ret_out[NLPSOL_F] = w[QPSOL_COST];
    ret_out[NLPSOL_G] = mtimes(v.at(2), w[QPSOL_X]) + v.at(3);
    ret_out[NLPSOL_LAM_X] = w[QPSOL_LAM_X];
    ret_out[NLPSOL_LAM_G] = w[QPSOL_LAM_A];
    ret_out[NLPSOL_LAM_P] = MX::nan(p.sparsity());
    return Function(name, ret_in, ret_out, nlpsol_in(), nlpsol_out());
  }

  Function qpsol(const std::string& name, const std::string& solver,
                 const XProblem& qp, const Dict& opts) {
    if (qp.is_sx) {
      return qpsol_nlp<SX>(name, solver, qp, opts);
    } else {
      return qpsol_nlp<MX>(name, solver, qp, opts);
    }
  }

  Function qpsol(const std::string& name, const std::string& solver,
                 const SXDict& qp, const Dict& opts) {
    return qpsol(name, solver, Nlpsol::map2problem(qp), opts);
  }

  Function qpsol(const std::string& name, const std::string& solver,
                 const MXDict& qp, const Dict& opts) {
    return qpsol(name, solver, Nlpsol::map2problem(qp), opts);
  }

  // Constructor
  Qpsol::Qpsol(const std::string& name,
               const std::map<std::string, Sparsity> &st)
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

    n_ = H_.size2();
    nc_ = A_.is_null() ? 0 : A_.size1();

    if (!A_.is_null()) {
      casadi_assert_message(A_.size2()==n_,
        "Got incompatible dimensions.   min          x'Hx + G'x s.t.   LBA <= Ax <= UBA :"
        << std::endl <<
        "H: " << H_.dim() << " - A: " << A_.dim() << std::endl <<
        "We need: H_.size2()==A_.size2()" << std::endl);
    }

    casadi_assert_message(H_.is_symmetric(),
      "Got incompatible dimensions.   min          x'Hx + G'x" << std::endl <<
      "H: " << H_.dim() <<
      "We need H square & symmetric" << std::endl);

    // Sparsity
    Sparsity x_sparsity = Sparsity::dense(n_, 1);
    Sparsity bounds_sparsity = Sparsity::dense(nc_, 1);
  }

  Sparsity Qpsol::get_sparsity_in(int ind) const {
    switch (static_cast<QpsolInput>(ind)) {
    case QPSOL_X0:
    case QPSOL_G:
    case QPSOL_LBX:
    case QPSOL_UBX:
    case QPSOL_LAM_X0:
      return get_sparsity_out(QPSOL_X);
    case QPSOL_LBA:
    case QPSOL_UBA:
      return get_sparsity_out(QPSOL_LAM_A);
    case QPSOL_A:
      return A_;
    case QPSOL_H:
      return H_;
    case QPSOL_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity Qpsol::get_sparsity_out(int ind) const {
    switch (static_cast<QpsolOutput>(ind)) {
    case QPSOL_COST:
      return Sparsity::scalar();
    case QPSOL_X:
    case QPSOL_LAM_X:
      return Sparsity::dense(n_, 1);
    case QPSOL_LAM_A:
      return Sparsity::dense(nc_, 1);
    case QPSOL_NUM_OUT: break;
    }
    return Sparsity();
  }

  void Qpsol::init(const Dict& opts) {
    // Call the init method of the base class
    FunctionInternal::init(opts);
  }

  Qpsol::~Qpsol() {
  }

  void Qpsol::checkInputs(const double* lbx, const double* ubx,
                          const double* lba, const double* uba) const {
    for (int i=0; i<n_; ++i) {
      double lb = lbx ? lbx[i] : 0., ub = ubx ? ubx[i] : 0.;
      casadi_assert_message(lb <= ub,
                            "LBX[" << i << "] <= UBX[" << i << "] was violated. "
                            << "Got LBX[" << i << "]=" << lb <<
                            " and UBX[" << i << "] = " << ub << ".");
    }
    for (int i=0; i<nc_; ++i) {
      double lb = lba ? lba[i] : 0., ub = uba ? uba[i] : 0.;
      casadi_assert_message(lb <= ub,
                            "LBA[" << i << "] <= UBA[" << i << "] was violated. "
                            << "Got LBA[" << i << "] = " << lb <<
                            " and UBA[" << i << "] = " << ub << ".");
    }
  }

  void Qpsol::generateNativeCode(std::ostream& file) const {
    casadi_error("Qpsol::generateNativeCode not defined for class "
                 << typeid(*this).name());
  }

  std::map<std::string, Qpsol::Plugin> Qpsol::solvers_;

  const std::string Qpsol::infix_ = "qpsol";

  double Qpsol::default_in(int ind) const {
    switch (ind) {
    case QPSOL_LBX:
    case QPSOL_LBA:
      return -std::numeric_limits<double>::infinity();
    case QPSOL_UBX:
    case QPSOL_UBA:
      return std::numeric_limits<double>::infinity();
    default:
      return 0;
    }
  }

} // namespace casadi




