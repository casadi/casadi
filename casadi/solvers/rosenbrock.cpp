/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#include "rosenbrock.hpp"
#include "rosenbrock_tableau.hpp"

namespace casadi {

  extern "C"
  int CASADI_INTEGRATOR_ROSENBROCK_EXPORT
      casadi_register_integrator_rosenbrock(Integrator::Plugin* plugin) {
    plugin->creator = Rosenbrock::creator;
    plugin->name = "rosenbrock";
    plugin->doc = Rosenbrock::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Rosenbrock::options_;
    plugin->deserialize = &Rosenbrock::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_INTEGRATOR_ROSENBROCK_EXPORT casadi_load_integrator_rosenbrock() {
    Integrator::registerPlugin(casadi_register_integrator_rosenbrock);
  }

  Rosenbrock::Rosenbrock(const std::string& name, const Function& dae, double t0,
      const std::vector<double>& tout) : RuntimeIntegrator(name, dae, t0, tout) {
  }

  Rosenbrock::~Rosenbrock() {
    clear_mem();
  }

  const Options Rosenbrock::options_
  = {{&RuntimeIntegrator::options_},
     {{"init_tol",
       {OT_DOUBLE,
        "Tolerance on the scaled Newton step when making the algebraic states "
        "consistent [1e-3]"}},
      {"init_max_iter",
       {OT_INT,
        "Maximum number of Newton iterations for consistent algebraic states [20]"}}
     }
  };

  void Rosenbrock::init(const Dict& opts) {
    // Call the base class init
    RuntimeIntegrator::init(opts);

    // Default options
    init_tol_ = 1e-3;
    init_max_iter_ = 20;

    // Read options
    for (auto&& op : opts) {
      if (op.first == "init_tol") {
        init_tol_ = op.second;
      } else if (op.first == "init_max_iter") {
        init_max_iter_ = op.second;
      }
    }

    // Reverse communication problem
    set_iteration_matrix();
    set_coefficients();
    set_ivp_prob();
    set_nominal_prob();
    set_rosenbrock_prob();

    // Work vectors of the kernel, and of the nominal integration
    alloc_ivp_work();
    casadi_int sz_arg = 0, sz_res = 0, sz_iw = 0, sz_w = 0;
    casadi_rosenbrock_work(&pr_, &sz_arg, &sz_res, &sz_iw, &sz_w);
    if (replay_) casadi_rosenbrock_work(&pr_nom_, &sz_arg, &sz_res, &sz_iw, &sz_w);
    alloc_w(sz_w, true);
  }

  void Rosenbrock::set_iteration_matrix() {
    // Iteration matrix M/(h gamma) - J of the nominal [x; z]
    sp_iter_ = blockcat(sp_jac_[IVP_ODE_X] + Sparsity::diag(nx1_),
      sp_jac_[IVP_ODE_Z], sp_jac_[IVP_ALG_X], sp_jac_[IVP_ALG_Z]);
    // Positions of the nonzeros of the blocks ODE_X, ALG_X, ODE_Z, ALG_Z, then the x-diagonal
    iter_map_.clear();
    const casadi_int blk[4] = {IVP_ODE_X, IVP_ALG_X, IVP_ODE_Z, IVP_ALG_Z};
    for (casadi_int k = 0; k < 4; ++k) {
      std::vector<casadi_int> r = sp_jac_[blk[k]].get_row(), c = sp_jac_[blk[k]].get_col();
      casadi_int r_off = blk[k] == IVP_ALG_X || blk[k] == IVP_ALG_Z ? nx1_ : 0;
      casadi_int c_off = blk[k] == IVP_ODE_Z || blk[k] == IVP_ALG_Z ? nx1_ : 0;
      for (size_t i = 0; i < r.size(); ++i) {
        iter_map_.push_back(sp_iter_.get_nz(r[i] + r_off, c[i] + c_off));
      }
    }
    for (casadi_int i = 0; i < nx1_; ++i) iter_map_.push_back(sp_iter_.get_nz(i, i));
    // Symbolic sparse QR factorizations
    sp_iter_.qr_sparse(sp_iter_v_, sp_iter_r_, iter_prinv_, iter_pc_);
    if (nz1_ > 0) sp_jac_[IVP_ALG_Z].qr_sparse(sp_alg_v_, sp_alg_r_, alg_prinv_, alg_pc_);
  }

  void Rosenbrock::set_coefficients() {
    // Coefficients in original form (rosenbrock_tableau.hpp), Gamma lower triangular
    using namespace rosenbrock_tableau;
    const casadi_int s = n_stages;
    // Inverse of Gamma = (gamma_ij), lower triangular, in long double
    std::vector<long double> Gi(s * s, 0);
    for (casadi_int j = 0; j < s; ++j) {
      Gi[j * s + j] = 1 / gamma;
      for (casadi_int i = j + 1; i < s; ++i) {
        long double acc = 0;
        for (casadi_int k = j; k < i; ++k) acc += gamma_ij[i][k] * Gi[k * s + j];
        Gi[i * s + j] = -acc / gamma;
      }
    }
    // a = A Gamma^-1, c = diag(1/gamma) - Gamma^-1, m = Gamma^-T b, m_err = Gamma^-T (b - bhat)
    a_.assign(s * s, 0);
    c_.assign(s * s, 0);
    alpha_.assign(s, 0);
    gamma_sum_.assign(s, 0);
    m_.assign(s, 0);
    m_err_.assign(s, 0);
    for (casadi_int i = 0; i < s; ++i) {
      long double alpha_i = 0, gamma_i = gamma;
      for (casadi_int j = 0; j < i; ++j) {
        alpha_i += alpha[i][j];
        gamma_i += gamma_ij[i][j];
        long double aij = 0;
        for (casadi_int k = j; k < i; ++k) aij += alpha[i][k] * Gi[k * s + j];
        a_[i * s + j] = static_cast<double>(aij);
        c_[i * s + j] = static_cast<double>(-Gi[i * s + j]);
      }
      alpha_[i] = static_cast<double>(alpha_i);
      gamma_sum_[i] = static_cast<double>(gamma_i);
      long double mi = 0, mhi = 0;
      for (casadi_int k = i; k < s; ++k) {
        mi += b[k] * Gi[k * s + i];
        mhi += bhat[k] * Gi[k * s + i];
      }
      m_[i] = static_cast<double>(mi);
      m_err_[i] = static_cast<double>(mi - mhi);
    }
    order_err_ = order_embedded;
    // Dense output weights in the same form: m(theta) = Gamma^-T b(theta)
    n_dense_ = n_dense;
    m_dense_.assign(n_dense * s, 0);
    for (casadi_int q = 0; q < n_dense; ++q) {
      for (casadi_int i = 0; i < s; ++i) {
        long double mi = 0;
        for (casadi_int k = i; k < s; ++k) mi += b_dense[q][k] * Gi[k * s + i];
        m_dense_[q * s + i] = static_cast<double>(mi);
      }
    }
  }

  void Rosenbrock::set_rosenbrock_prob() {
    // Problem structure of the method
    pr_.ivp = &p_;
    pr_.gamma = static_cast<double>(rosenbrock_tableau::gamma);
    pr_.c = get_ptr(c_);
    pr_.gamma_sum = get_ptr(gamma_sum_);
    pr_.sp_iter = sp_iter_;
    pr_.iter_map = get_ptr(iter_map_);
    pr_.sp_iter_v = sp_iter_v_;
    pr_.sp_iter_r = sp_iter_r_;
    pr_.iter_prinv = get_ptr(iter_prinv_);
    pr_.iter_pc = get_ptr(iter_pc_);
    pr_.sp_alg_v = nz1_ > 0 ? static_cast<const casadi_int*>(sp_alg_v_) : nullptr;
    pr_.sp_alg_r = nz1_ > 0 ? static_cast<const casadi_int*>(sp_alg_r_) : nullptr;
    pr_.alg_prinv = get_ptr(alg_prinv_);
    pr_.alg_pc = get_ptr(alg_pc_);
    pr_.init_tol = init_tol_;
    pr_.init_max_iter = init_max_iter_;
    casadi_rosenbrock_setup(&pr_);
    // Method of the original integrator, if its steps are replayed
    if (replay_) {
      pr_nom_ = static_cast<const Rosenbrock*>(root())->pr_;
      pr_nom_.ivp = &p_nom_;
      casadi_rosenbrock_setup(&pr_nom_);
    }
  }

  void Rosenbrock::set_work(void* mem, const double**& arg, double**& res,
      casadi_int*& iw, double*& w) const {
    auto m = static_cast<RosenbrockMemory*>(mem);
    RuntimeIntegrator::set_work(mem, arg, res, iw, w);
    m->r.prob = &pr_;
    casadi_rosenbrock_set_work(&m->r, &arg, &res, &iw, &w);
    if (replay_) {
      m->r_nom.prob = &pr_nom_;
      casadi_rosenbrock_set_work(&m->r_nom, &arg, &res, &iw, &w);
      m->r.nom = &m->r_nom;
    }
  }

  int Rosenbrock::kernel(RuntimeIntegratorMemory* mem) const {
    auto m = static_cast<RosenbrockMemory*>(mem);
    return casadi_rosenbrock(&m->d, &m->r);
  }

  void Rosenbrock::codegen_method(CodeGenerator& g) const {
    g.add_auxiliary(CodeGenerator::AUX_ROSENBROCK);
    codegen_rosenbrock_prob(g, "pe", "p");
    g << "e.prob = &pe;\n";
    g << "casadi_rosenbrock_set_work(&e, &arg, &res, &iw, &w);\n";
    if (replay_) {
      static_cast<const Rosenbrock*>(root())->codegen_rosenbrock_prob(g, "pen", "pn");
      g << "en.prob = &pen;\n";
      g << "casadi_rosenbrock_set_work(&en, &arg, &res, &iw, &w);\n";
      g << "e.nom = &en;\n";
    }
  }

  void Rosenbrock::codegen_rosenbrock_prob(CodeGenerator& g, const std::string& pe,
      const std::string& p) const {
    // Problem structure of the method in variable pe, as set_rosenbrock_prob
    auto set = [&](const std::string& field, const std::string& value) -> void {
      g << pe << "." << field << " = " << value << ";\n";
    };
    set("ivp", "&" + p);
    set("gamma", g.constant(pr_.gamma));
    set("c", g.constant(c_));
    set("gamma_sum", g.constant(gamma_sum_));
    set("sp_iter", g.sparsity(sp_iter_));
    set("iter_map", g.constant(iter_map_));
    set("sp_iter_v", g.sparsity(sp_iter_v_));
    set("sp_iter_r", g.sparsity(sp_iter_r_));
    set("iter_prinv", g.constant(iter_prinv_));
    set("iter_pc", g.constant(iter_pc_));
    set("sp_alg_v", nz1_ > 0 ? g.sparsity(sp_alg_v_) : "0");
    set("sp_alg_r", nz1_ > 0 ? g.sparsity(sp_alg_r_) : "0");
    set("alg_prinv", nz1_ > 0 ? g.constant(alg_prinv_) : "0");
    set("alg_pc", nz1_ > 0 ? g.constant(alg_pc_) : "0");
    set("init_tol", g.constant(pr_.init_tol));
    set("init_max_iter", str(pr_.init_max_iter));
    g << "casadi_rosenbrock_setup(&" << pe << ");\n";
  }

  Rosenbrock::Rosenbrock(DeserializingStream& s) : RuntimeIntegrator(s) {
    s.version("Rosenbrock", 1);
    s.unpack("Rosenbrock::init_tol", init_tol_);
    s.unpack("Rosenbrock::init_max_iter", init_max_iter_);
    set_iteration_matrix();
    set_coefficients();
    set_ivp_prob();
    set_nominal_prob();
    set_rosenbrock_prob();
  }

  void Rosenbrock::serialize_body(SerializingStream &s) const {
    RuntimeIntegrator::serialize_body(s);
    s.version("Rosenbrock", 1);
    s.pack("Rosenbrock::init_tol", init_tol_);
    s.pack("Rosenbrock::init_max_iter", init_max_iter_);
  }

} // namespace casadi
