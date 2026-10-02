/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
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



#include "acados_integrator.hpp"
#include "casadi/core/code_generator.hpp"
#include "casadi/core/serializing_stream.hpp"
#include <acados_runtime_str.h>

#include <acados/sim/sim_collocation_utils.h>

#include <algorithm>
#include <cstdlib>

namespace casadi {


  extern "C"
  int CASADI_INTEGRATOR_ACADOS_EXPORT
  casadi_register_integrator_acados(Integrator::Plugin* plugin) {
    plugin->creator = AcadosInterface::creator;
    plugin->name = "acados";
    plugin->doc = AcadosInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &AcadosInterface::options_;
    plugin->deserialize = &AcadosFunction::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_INTEGRATOR_ACADOS_EXPORT casadi_load_integrator_acados() {
    Integrator::registerPlugin(casadi_register_integrator_acados);
  }

  template<typename M>
  static std::vector<Function> casadi_acados_model_functions(const Function& dae, bool irk,
      bool p_diff, const std::string& suffix) {
    casadi_int nx = dae.numel_in(DYN_X), nz = dae.numel_in(DYN_Z), np = dae.numel_in(DYN_P);
    casadi_int nu = dae.numel_in(DYN_U);
    // acados controls w = [u; p] if p is differentiable, else p is passed by CasADi only
    M x = M::sym("x", nx), u = M::sym("u", nu), pd = M::sym("p", p_diff ? np : 0);
    M p = M::sym("p", p_diff ? 0 : np), w = vertcat(u, pd), z = M::sym("z", nz);
    casadi_int nwu = w.numel();
    std::vector<M> dae_in(DYN_NUM_IN);
    dae_in[DYN_T] = M::zeros(dae.sparsity_in(DYN_T));
    dae_in[DYN_X] = x;
    dae_in[DYN_Z] = z;
    dae_in[DYN_P] = p_diff ? pd : p;
    dae_in[DYN_U] = u;
    std::vector<M> dae_out = dae(dae_in);
    M f = dae_out.at(DYN_ODE);
    M xw = vertcat(x, w);
    // Same signatures as the acados model functions, with the CasADi-only p as last input
    std::vector<std::pair<std::string, std::pair<std::vector<M>, std::vector<M>>>> d;
    if (irk) {
      M xdot = M::sym("xdot", nx), t = M::sym("t");
      M f_impl = vertcat(xdot - f, dae_out.at(DYN_ALG));
      M jac_x = M::jacobian(f_impl, x), jac_xdot = M::jacobian(f_impl, xdot);
      M jac_u = M::jacobian(f_impl, w), jac_z = M::jacobian(f_impl, z);
      std::vector<M> in = {x, xdot, w, z, t, p};
      d.push_back({"impl_dae_fun", {in, {f_impl}}});
      d.push_back({"impl_dae_fun_jac_x_xdot_z", {in, {f_impl, jac_x, jac_xdot, jac_z}}});
      d.push_back({"impl_dae_jac_x_xdot_u_z", {in, {jac_x, jac_xdot, jac_u, jac_z}}});
      M xxdotzu = vertcat(std::vector<M>{x, xdot, z, w});
      M mult = M::sym("multiplier", nx + nz);
      M adj = M::jtimes(f_impl, xxdotzu, mult, true);
      d.push_back({"impl_dae_hess", {{x, xdot, w, z, mult, t, p},
        {M::jacobian(adj, xxdotzu)}}});
    } else {
      M Sx = M::sym("Sx", nx, nx), Sp = M::sym("Sp", nx, nwu), lam = M::sym("lambdaX", nx);
      M vdeX = M::jtimes(f, x, Sx);
      M vdeP = M::jacobian(f, w) + M::jtimes(f, x, Sp);
      M adj = M::jtimes(f, xw, lam, true);
      M S_forw = vertcat(horzcat(Sx, Sp), horzcat(M::zeros(nwu, nx), M::eye(nwu)));
      M hess = mtimes(S_forw.T(), M::jtimes(adj, xw, S_forw));
      std::vector<M> hess2;
      for (casadi_int j = 0; j < nx + nwu; ++j) {
        for (casadi_int i = j; i < nx + nwu; ++i) hess2.push_back(hess(i, j));
      }
      d.push_back({"expl_ode_fun", {{x, w, p}, {f}}});
      d.push_back({"expl_vde_forw", {{x, Sx, Sp, w, p}, {f, vdeX, vdeP}}});
      d.push_back({"expl_vde_adj", {{x, lam, w, p}, {adj}}});
      d.push_back({"expl_ode_hess", {{x, Sx, Sp, lam, w, p}, {adj, vertcat(hess2)}}});
    }
    std::vector<Function> ret;
    for (auto&& e : d) {
      std::vector<M> out = e.second.second;
      for (auto&& o : out) o = densify(o);
      ret.push_back(Function("acados_" + e.first + suffix, e.second.first, out));
    }
    return ret;
  }

  const std::vector<std::string> AcadosModel::option_names = {
    "scheme",            // erk|irk (default irk)
    "num_stages",        // default 4
    "num_steps",         // per integration interval, default 1
    "newton_iter",       // irk, default 3
    "collocation_type"}; // gauss_legendre|gauss_radau_iia

  void AcadosModel::init(const Function& dae, const Dict& opts, bool p_diff) {
    this->p_diff = p_diff;
    irk = true;
    num_stages = 4;
    num_steps = 1;
    newton_iter = 3;
    collocation_type = "gauss_legendre";
    for (auto&& op : opts) {
      if (op.first=="scheme") {
        irk = op.second.to_string()=="irk";
      } else if (op.first=="num_stages") {
        num_stages = op.second;
      } else if (op.first=="num_steps") {
        num_steps = op.second;
      } else if (op.first=="newton_iter") {
        newton_iter = op.second;
      } else if (op.first=="collocation_type") {
        collocation_type = op.second.to_string();
      } else {
        casadi_error("acados: unknown option " + op.first);
      }
    }
    nx = dae.numel_in(DYN_X);
    np = dae.numel_in(DYN_P);
    nu = dae.numel_in(DYN_U);
    nz = dae.numel_in(DYN_Z);
    npd = p_diff ? np : 0;
    casadi_assert(irk || nz == 0, "acados: algebraic states require scheme irk");
    casadi_assert(dae.numel_out(DYN_QUAD) == 0, "acados: quadratures not supported");
    casadi_assert(dae.numel_out(DYN_ZERO) == 0, "acados: events not supported");
    casadi_assert(dae.numel_in(DYN_T) == 0 || dae.sparsity_jac(DYN_T, DYN_ODE).nnz() == 0,
      "acados: time-dependent dynamics not supported");
    this->dae = dae;
    if (irk) {
      fields = {"impl_dae_fun", "impl_dae_fun_jac_x_xdot_z", "impl_dae_jac_x_xdot_u_z",
                "impl_dae_hess"};
    } else {
      fields = {"expl_ode_fun", "expl_vde_forw", "expl_vde_adj", "expl_ode_hess"};
    }
  }

  const Function& AcadosModel::fcn(bool v, casadi_int k) const {
    if (fcns[v][k].is_null()) {
      std::string suffix = v ? "_p" : "";
      std::vector<Function> f = dae.is_a("SXFunction") ?
        casadi_acados_model_functions<SX>(dae, irk, v, suffix) :
        casadi_acados_model_functions<MX>(dae, irk, v, suffix);
      for (casadi_int i = 0; i < 4; ++i) if (fcns[v][i].is_null()) fcns[v][i] = f[i];
    }
    return fcns[v][k];
  }

  std::string AcadosModel::fcn_name(bool v, casadi_int k) const {
    return "acados_" + fields.at(k) + (v ? "_p" : "");
  }

  void AcadosModel::register_functions(OracleFunction* owner, bool v, casadi_int n) const {
    for (casadi_int k = 0; k < n; ++k) owner->set_function(fcn(v, k), fcn_name(v, k), true);
  }

  int AcadosModel::collocation_type_enum() const {
    if (!irk) return EXPLICIT_RUNGE_KUTTA;
    return collocation_type=="gauss_radau_iia" ? GAUSS_RADAU_IIA : GAUSS_LEGENDRE;
  }

  AcadosInterface::AcadosInterface(const std::string& name, const Function& dae,
      double t0, const std::vector<double>& tout) : Integrator(name, dae, t0, tout) {
  }

  AcadosInterface::~AcadosInterface() {
    clear_mem();
  }

  const Options AcadosInterface::options_
  = {{&Integrator::options_},
     {{"scheme",
       {OT_STRING, "erk|irk (default irk)"}},
      {"num_stages",
       {OT_INT, "Number of Runge-Kutta stages (default 4)"}},
      {"num_steps",
       {OT_INT, "Number of integration steps (default 1)"}},
      {"newton_iter",
       {OT_INT, "Number of Newton iterations, irk (default 3)"}},
      {"collocation_type",
       {OT_STRING, "gauss_legendre|gauss_radau_iia (default gauss_legendre)"}}
     }
  };

  Function AcadosInterface::create_advanced(const Dict& opts) {
    // Initialize as a regular Integrator: checks the options, fixes the I/O signature
    Function temp = Function::create(this, opts);
    Dict model_opts, fun_opts, nom_opts;
    const auto& on = AcadosModel::option_names;
    for (auto&& op : opts) {
      if (std::find(on.begin(), on.end(), op.first) != on.end()) {
        model_opts[op.first] = op.second;
      } else if (op.first == "is_diff_in" || op.first == "is_diff_out") {
        // Integrator signature only, not that of the derivative functions
        nom_opts[op.first] = op.second;
      } else if (OracleFunction::options_.find(op.first)) {
        fun_opts[op.first] = op.second;
      }
    }
    // p not differentiable unless declared otherwise: sensitivities w.r.t. p cost
    // O((nx+nu+np)^2) memory in the hessian
    bool p_diff = false;
    if (nom_opts.find("is_diff_in") == nom_opts.end()) {
      std::vector<bool> is_diff_in(INTEGRATOR_NUM_IN, true);
      is_diff_in[INTEGRATOR_P] = false;
      nom_opts["is_diff_in"] = is_diff_in;
    } else {
      p_diff = temp.is_diff_in(INTEGRATOR_P);
    }
    auto model = std::make_shared<AcadosModel>();
    model->init(oracle_, model_opts, p_diff);
    std::vector<Sparsity> sp_in(INTEGRATOR_NUM_IN), sp_out(INTEGRATOR_NUM_OUT);
    for (casadi_int i = 0; i < INTEGRATOR_NUM_IN; ++i) sp_in[i] = sparsity_in(i);
    for (casadi_int i = 0; i < INTEGRATOR_NUM_OUT; ++i) sp_out[i] = sparsity_out(i);
    // Interval lengths
    std::vector<double> T;
    double t = t0_;
    for (double tk : tout_) {
      T.push_back(tk - t);
      t = tk;
    }
    return Function::create(new AcadosFunction(name_, oracle_, model, model_opts,
      fun_opts, T, AcadosFunction::NOM, 0, false, sp_in, sp_out,
      integrator_in(), integrator_out()), combine(nom_opts, fun_opts));
  }

  int AcadosInterface::advance_noevent(IntegratorMemory* mem) const {
    casadi_error("acados: not evaluated as an Integrator");
    return 1;
  }

  void AcadosInterface::resetB(IntegratorMemory* mem) const {
    casadi_error("acados: not evaluated as an Integrator");
  }

  void AcadosInterface::impulseB(IntegratorMemory* mem,
      const double* adj_x, const double* adj_z, const double* adj_q) const {
    casadi_error("acados: not evaluated as an Integrator");
  }

  void AcadosInterface::retreat(IntegratorMemory* mem, const double* u,
      double* adj_x, double* adj_p, double* adj_u) const {
    casadi_error("acados: not evaluated as an Integrator");
  }

  AcadosFunction::AcadosFunction(const std::string& name, const Function& dae,
      const std::shared_ptr<AcadosModel>& model, const Dict& model_opts, const Dict& fun_opts,
      const std::vector<double>& T, Mode mode, casadi_int nadj, bool pv,
      const std::vector<Sparsity>& sp_in_nom, const std::vector<Sparsity>& sp_out_nom,
      const std::vector<std::string>& names_in, const std::vector<std::string>& names_out)
      : OracleFunction(name, dae), model_(model), model_opts_(model_opts), fun_opts_(fun_opts),
        T_(T), mode_(mode), nadj_(nadj), pv_(pv), sp_in_nom_(sp_in_nom), sp_out_nom_(sp_out_nom),
        names_in_(names_in), names_out_(names_out) {
  }

  AcadosFunction::~AcadosFunction() {
    clear_mem();
  }

  casadi_int AcadosFunction::acados_offset(casadi_int i) const {
    if (i == INTEGRATOR_X0) return 0;
    if (i == INTEGRATOR_U) return model_->nx;
    if (i == INTEGRATOR_P && pv_) return model_->nx + model_->nu * nt();
    return -1;
  }

  // Sparsity of input j of the ADJ function with nadj directions
  static Sparsity casadi_acados_adj_in(const std::vector<Sparsity>& sp_in,
      const std::vector<Sparsity>& sp_out, casadi_int j, casadi_int nadj) {
    casadi_int ni = sp_in.size(), no = sp_out.size();
    if (j < ni) return sp_in[j];
    if (j < ni + no) return Sparsity(sp_out[j - ni].size());  // nominal outputs not needed
    return repmat(sp_out[j - ni - no], 1, nadj);
  }

  Sparsity AcadosFunction::get_sparsity_in(casadi_int i) {
    casadi_int ni = INTEGRATOR_NUM_IN, no = INTEGRATOR_NUM_OUT;
    switch (mode_) {
      case NOM: return sp_in_nom_.at(i);
      case JAC: return i < ni ? sp_in_nom_.at(i) : Sparsity(sp_out_nom_.at(i - ni).size());
      case ADJ: return casadi_acados_adj_in(sp_in_nom_, sp_out_nom_, i, nadj_);
      case HESS:
        if (i < ni + 2 * no) return casadi_acados_adj_in(sp_in_nom_, sp_out_nom_, i, 1);
        return Sparsity(sp_in_nom_.at(i - ni - 2 * no).size());
    }
    return Sparsity();
  }

  Sparsity AcadosFunction::get_sparsity_out(casadi_int i) {
    casadi_int ni = INTEGRATOR_NUM_IN, no = INTEGRATOR_NUM_OUT, nx = model_->nx;
    switch (mode_) {
      case NOM: return sp_out_nom_.at(i);
      case JAC:
      {
        casadi_int o = i / ni, j = i % ni;
        casadi_int nr = sp_out_nom_[o].numel(), nc = sp_in_nom_[j].numel();
        if (o == INTEGRATOR_XF && acados_offset(j) >= 0) return Sparsity::dense(nr, nc);
        return Sparsity(nr, nc);
      }
      case ADJ:
      {
        casadi_int nr = sp_in_nom_.at(i).size1(), nc = sp_in_nom_.at(i).size2() * nadj_;
        return acados_offset(i) >= 0 ? Sparsity::dense(nr, nc) : Sparsity(nr, nc);
      }
      case HESS:
      {
        casadi_int nj = ni + 2 * no;
        casadi_int o = i / nj, j = i % nj;
        casadi_int nr = sp_in_nom_[o].numel();
        casadi_int nc = casadi_acados_adj_in(sp_in_nom_, sp_out_nom_, j, 1).numel();
        if (acados_offset(o) >= 0) {
          if (j < ni && acados_offset(j) >= 0) return Sparsity::dense(nr, nc);
          if (j == adj_seed_xf()) return Sparsity::dense(nr, nc);
        }
        return Sparsity(nr, nc);
      }
    }
    return Sparsity();
  }

  void AcadosFunction::set_acados_prob() {
    p_.sim.irk = model_->irk;
    p_.sim.num_stages = model_->num_stages;
    p_.sim.num_steps = model_->num_steps;
    p_.sim.newton_iter = model_->newton_iter;
    p_.sim.collocation_type = model_->collocation_type_enum();
    p_.sim.nz = model_->nz;
    p_.sim.sim_bytes = 0;
    for (casadi_int k = 0; k < nfun(); ++k) {
      const Function& f = model_->fcn(pv_, k);
      p_.sim.cb[k] = OracleCallback(f.name(), this);
      p_.sim.field[k] = model_->fields[k].c_str();
      dims_[k] = {f.n_in(), f.n_out()};
      for (casadi_int i = 0; i < f.n_in(); ++i) {
        dims_[k].push_back(f.size1_in(i));
        dims_[k].push_back(f.size2_in(i));
      }
      for (casadi_int i = 0; i < f.n_out(); ++i) {
        dims_[k].push_back(f.size1_out(i));
        dims_[k].push_back(f.size2_out(i));
      }
      p_.sim.dims[k] = get_ptr(dims_[k]);
    }
    p_.chain.nx = model_->nx;
    p_.chain.nu = model_->nu;
    p_.chain.np = model_->np;
    p_.chain.npd = pv_ ? model_->npd : 0;
    p_.chain.nt = nt();
    p_.chain.T = get_ptr(T_);
    p_.mode = mode_;
    p_.nadj = nadj_;
    p_.nan = nan;
    casadi_acados_setup(&p_);
  }

  casadi_int AcadosFunction::measure_sim() const {
    // Model functions are not evaluated while laying out
    casadi_acados_sim_data<double> d;
    d.prob = &p_.sim;
    d.oracle = nullptr;
    d.z0 = nullptr;
    d.measure = 1;
    for (casadi_int n = 1 << 16; ; n *= 4) {
      std::vector<char> buf(n);
      d.raw = buf.data();
      d.nraw = n;
      casadi_int used = casadi_acados_sim_init(&d);
      if (used >= 0) return used;
      casadi_assert(n < (1LL << 34), "acados: cannot lay out the sim objects");
    }
  }

  void AcadosFunction::init(const Dict& opts) {
    OracleFunction::init(opts);
    model_->register_functions(this, pv_, nfun());
    set_acados_prob();
    p_.sim.sim_bytes = measure_sim();
    casadi_int sz_w = 0;
    casadi_acados_work(&p_, &sz_w);
    alloc_w(sz_w, true);
  }

  void AcadosFunction::set_work(void* mem, const double**& arg, double**& res,
      casadi_int*& iw, double*& w) const {
    auto* m = static_cast<AcadosFunctionMemory*>(mem);
    OracleFunction::set_work(mem, arg, res, iw, w);
    m->d.prob = &p_;
    m->d.sim.oracle = &m->d_oracle;
    casadi_acados_set_work(&m->d, &arg, &res, &iw, &w);
    // Model functions are evaluated with this memory's oracle work
    m->d_oracle.m = static_cast<void*>(m);
  }

  int AcadosFunction::init_mem(void* mem) const {
    return OracleFunction::init_mem(mem);
  }

  void AcadosFunction::free_mem(void* mem) const {
    delete static_cast<AcadosFunctionMemory*>(mem);
  }

  int AcadosFunction::eval(const double** arg, double** res, casadi_int* iw, double* w,
      void* mem) const {
    auto* m = static_cast<AcadosFunctionMemory*>(mem);
    setup(m, arg + n_in_, res + n_out_, iw, w);
    int flag = casadi_acados_eval(&m->d, arg, res);
    casadi_assert(flag != 2, "acados: work memory too small for the sim objects");
    join_results(m);
    return flag;
  }

  Function AcadosFunction::derivative(const std::string& name, Mode mode, casadi_int nadj,
      bool pv, const std::vector<std::string>& inames, const std::vector<std::string>& onames,
      const Dict& opts) const {
    return Function::create(new AcadosFunction(name, oracle_, model_, model_opts_, fun_opts_,
      T_, mode, nadj, pv, sp_in_nom_, sp_out_nom_, inames, onames), combine(opts, fun_opts_));
  }

  Function AcadosFunction::split_p(const std::string& name, const Function& f0,
      const Function& f1, const std::vector<bool>& from_f1,
      const std::vector<std::string>& inames, const std::vector<std::string>& onames,
      const Dict& opts) {
    std::vector<MX> arg = f0.mx_in(), r0 = f0(arg), r1 = f1(arg);
    for (casadi_int i = 0; i < r0.size(); ++i) if (from_f1[i]) r0[i] = r1[i];
    // Inlined: a caller not using the blocks w.r.t. p does not evaluate (or size) f1
    Dict wopts = opts;
    wopts["always_inline"] = true;
    return Function(name, arg, r0, inames, onames, wopts);
  }

  Function AcadosFunction::get_jacobian(const std::string& name,
      const std::vector<std::string>& inames, const std::vector<std::string>& onames,
      const Dict& opts) const {
    Mode m = mode_==NOM ? JAC : HESS;
    if (pv_ || model_->npd == 0) return derivative(name, m, nadj_, pv_, inames, onames, opts);
    // Columns w.r.t. p (of nonzero rows) with p in the acados controls
    casadi_int ni = INTEGRATOR_NUM_IN, nj = m==JAC ? ni : ni + 2 * INTEGRATOR_NUM_OUT;
    std::vector<bool> from_f1(onames.size());
    for (casadi_int i = 0; i < from_f1.size(); ++i) {
      from_f1[i] = i % nj == INTEGRATOR_P && (m==JAC || i / nj != INTEGRATOR_P);
    }
    return split_p(name, derivative(name + "_0", m, nadj_, false, inames, onames, Dict()),
      derivative(name + "_p", m, nadj_, true, inames, onames, Dict()), from_f1, inames, onames,
      opts);
  }

  Function AcadosFunction::get_reverse(casadi_int nadj, const std::string& name,
      const std::vector<std::string>& inames, const std::vector<std::string>& onames,
      const Dict& opts) const {
    if (model_->npd == 0) return derivative(name, ADJ, nadj, false, inames, onames, opts);
    // Adjoint of p with p in the acados controls
    std::vector<bool> from_f1(onames.size());
    from_f1.at(INTEGRATOR_P) = true;
    return split_p(name, derivative(name + "_0", ADJ, nadj, false, inames, onames, Dict()),
      derivative(name + "_p", ADJ, nadj, true, inames, onames, Dict()), from_f1, inames, onames,
      opts);
  }

  void AcadosFunction::codegen_declarations(CodeGenerator& g) const {
    // Model functions, called by acados through casadi_acados_ext_eval and calc_function
    g.add_auxiliary(CodeGenerator::AUX_ORACLE);
    g.add_auxiliary(CodeGenerator::AUX_ORACLE_CALLBACK);
    for (casadi_int k = 0; k < nfun(); ++k) g.add_dependency(model_->fcn(pv_, k));
  }

  void AcadosFunction::codegen_body(CodeGenerator& g) const {
    g.add_include("string.h");
    g.add_include("acados_c/sim_interface.h");
    g.add_include("acados/sim/sim_erk_integrator.h");
    g.add_include("acados/sim/sim_irk_integrator.h");
    g.add_include("blasfeo_d_aux.h");
    g.add_auxiliary(CodeGenerator::AUX_COPY);
    g.add_auxiliary(CodeGenerator::AUX_CLEAR);
    g.add_auxiliary(CodeGenerator::AUX_AXPY);
    g.add_auxiliary(CodeGenerator::AUX_FILL);
    g.add_auxiliary(CodeGenerator::AUX_NAN);
    if (g.auxiliaries.str().find("struct casadi_acados_sim_data {") == std::string::npos) {
      g.auxiliaries << g.sanitize_source(acados_sim_str, {"casadi_real"});
      g.auxiliaries << g.sanitize_source(acados_chain_str, {"casadi_real"});
      g.auxiliaries << g.sanitize_source(acados_runtime_str, {"casadi_real"});
    }
    codegen_body_enter(g);
    g.local("p", "struct casadi_acados_prob");
    g.local("d", "struct casadi_acados_data");
    g.local("oarg", "const casadi_real", "**");
    g.local("ores", "casadi_real", "**");
    // Problem: same values as p_ in the VM (set_acados_prob)
    g << "p.sim.irk = " << p_.sim.irk << ";\n";
    g << "p.sim.num_stages = " << p_.sim.num_stages << ";\n";
    g << "p.sim.num_steps = " << p_.sim.num_steps << ";\n";
    g << "p.sim.newton_iter = " << p_.sim.newton_iter << ";\n";
    g << "p.sim.collocation_type = " << p_.sim.collocation_type << ";\n";
    g << "p.sim.nz = " << p_.sim.nz << ";\n";
    g << "p.sim.sim_bytes = " << p_.sim.sim_bytes << ";\n";
    for (casadi_int k = 0; k < nfun(); ++k) {
      g.setup_callback("p.sim.cb[" + str(k) + "]", model_->fcn(pv_, k));
      g << "p.sim.field[" << k << "] = \"" << p_.sim.field[k] << "\";\n";
      g << "p.sim.dims[" << k << "] = " << g.constant(dims_[k]) << ";\n";
    }
    g << "p.chain.nx = " << p_.chain.nx << ";\n";
    g << "p.chain.nu = " << p_.chain.nu << ";\n";
    g << "p.chain.np = " << p_.chain.np << ";\n";
    g << "p.chain.npd = " << p_.chain.npd << ";\n";
    g << "p.chain.nt = " << p_.chain.nt << ";\n";
    g << "p.chain.T = " << g.constant(T_) << ";\n";
    g << "p.mode = " << p_.mode << ";\n";
    g << "p.nadj = " << p_.nadj << ";\n";
    g << "p.nan = casadi_nan;\n";
    g << "casadi_acados_setup(&p);\n";
    // Work: acados glue and chain kernel, then the model functions (oracle)
    g << "d.prob = &p;\n";
    g << "d.sim.oracle = &d_oracle;\n";
    g << "casadi_acados_set_work(&d, &arg, &res, &iw, &w);\n";
    g << "oarg = arg + " << n_in_ << ";\n";
    g << "ores = res + " << n_out_ << ";\n";
    g << "casadi_oracle_set_work(&d_oracle, &oarg, &ores, &iw, &w);\n";
    g << "if (casadi_acados_eval(&d, arg, res)) return 1;\n";
    codegen_body_exit(g);
    g << "return 0;\n";
  }

  void AcadosFunction::serialize_type(SerializingStream &s) const {
    OracleFunction::serialize_type(s);
    PluginInterface<Integrator>::serialize_type(s);
  }

  void AcadosFunction::serialize_body(SerializingStream &s) const {
    OracleFunction::serialize_body(s);
    s.version("AcadosFunction", 1);
    s.pack("AcadosFunction::model_opts", model_opts_);
    s.pack("AcadosFunction::fun_opts", fun_opts_);
    s.pack("AcadosFunction::p_diff", model_->p_diff);
    s.pack("AcadosFunction::T", T_);
    s.pack("AcadosFunction::mode", static_cast<casadi_int>(mode_));
    s.pack("AcadosFunction::nadj", nadj_);
    s.pack("AcadosFunction::sp_in_nom", sp_in_nom_);
    s.pack("AcadosFunction::sp_out_nom", sp_out_nom_);
    s.pack("AcadosFunction::names_in", names_in_);
    s.pack("AcadosFunction::names_out", names_out_);
    s.pack("AcadosFunction::pv", pv_);
    s.pack("AcadosFunction::sim_bytes", static_cast<casadi_int>(p_.sim.sim_bytes));
  }

  AcadosFunction::AcadosFunction(DeserializingStream& s) : OracleFunction(s) {
    s.version("AcadosFunction", 1);
    s.unpack("AcadosFunction::model_opts", model_opts_);
    s.unpack("AcadosFunction::fun_opts", fun_opts_);
    bool p_diff;
    s.unpack("AcadosFunction::p_diff", p_diff);
    s.unpack("AcadosFunction::T", T_);
    casadi_int mode;
    s.unpack("AcadosFunction::mode", mode);
    mode_ = static_cast<Mode>(mode);
    s.unpack("AcadosFunction::nadj", nadj_);
    s.unpack("AcadosFunction::sp_in_nom", sp_in_nom_);
    s.unpack("AcadosFunction::sp_out_nom", sp_out_nom_);
    s.unpack("AcadosFunction::names_in", names_in_);
    s.unpack("AcadosFunction::names_out", names_out_);
    s.unpack("AcadosFunction::pv", pv_);
    casadi_int sim_bytes;
    s.unpack("AcadosFunction::sim_bytes", sim_bytes);
    // Rebuild the acados glue; the model functions themselves were deserialized
    model_ = std::make_shared<AcadosModel>();
    model_->init(oracle_, model_opts_, p_diff);
    for (casadi_int k = 0; k < nfun(); ++k) {
      model_->fcns[pv_][k] = get_function(model_->fcn_name(pv_, k));
    }
    set_acados_prob();
    // Layout size consistent with the serialized work sizes
    p_.sim.sim_bytes = sim_bytes;
  }

  ProtoFunction* AcadosFunction::deserialize(DeserializingStream& s) {
    return new AcadosFunction(s);
  }

} // namespace casadi
