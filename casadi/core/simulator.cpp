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


#include "simulator_impl.hpp"
#include "casadi_misc.hpp"

using namespace std;
namespace casadi {

  bool has_simulator(const string& name) {
    return Simulator::has_plugin(name);
  }

  void load_simulator(const string& name) {
    Simulator::load_plugin(name);
  }

  string doc_simulator(const string& name) {
    return Simulator::getPlugin(name).doc;
  }

  Function simulator(const string& name, const string& solver,
                      const SXDict& dae, const Dict& opts) {
    return simulator(name, solver, Simulator::map2oracle("dae", dae), opts);
  }

  Function simulator(const string& name, const string& solver,
                      const MXDict& dae, const Dict& opts) {
    return simulator(name, solver, Simulator::map2oracle("dae", dae), opts);
  }

  Function simulator(const string& name, const string& solver,
                      const Function& dae, const Dict& opts) {
    // Make sure that dae is sound
    if (dae.has_free()) {
      casadi_error("Cannot create '" + name + "' since " + str(dae.get_free()) + " are free.");
    }
    Simulator* intg = Simulator::getPlugin(solver).creator(name, dae);
    return intg->create_advanced(opts);
  }

  vector<string> simulator_in() {
    vector<string> ret(simulator_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=simulator_in(i);
    return ret;
  }

  vector<string> simulator_out() {
    vector<string> ret(simulator_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=simulator_out(i);
    return ret;
  }

  string simulator_in(casadi_int ind) {
    switch (static_cast<SimulatorInput>(ind)) {
    case SIMULATOR_X0:  return "x0";
    case SIMULATOR_P:   return "p";
    case SIMULATOR_Z0:  return "z0";
    case SIMULATOR_RX0: return "rx0";
    case SIMULATOR_RP:  return "rp";
    case SIMULATOR_RZ0: return "rz0";
    case SIMULATOR_NUM_IN: break;
    }
    return string();
  }

  string simulator_out(casadi_int ind) {
    switch (static_cast<SimulatorOutput>(ind)) {
    case SIMULATOR_XF:  return "xf";
    case SIMULATOR_QF:  return "qf";
    case SIMULATOR_ZF:  return "zf";
    case SIMULATOR_RXF: return "rxf";
    case SIMULATOR_RQF: return "rqf";
    case SIMULATOR_RZF: return "rzf";
    case SIMULATOR_NUM_OUT: break;
    }
    return string();
  }

  casadi_int simulator_n_in() {
    return SIMULATOR_NUM_IN;
  }

  casadi_int simulator_n_out() {
    return SIMULATOR_NUM_OUT;
  }

  Simulator::Simulator(const std::string& name, const Function& oracle)
    : OracleFunction(name, oracle) {

    // Negative number of parameters for consistancy checking
    np_ = -1;

    // Default options
    print_stats_ = false;
    output_t0_ = false;
  }

  Simulator::~Simulator() {
  }

  Sparsity Simulator::get_sparsity_in(casadi_int i) {
    switch (static_cast<SimulatorInput>(i)) {
    case SIMULATOR_X0: return x();
    case SIMULATOR_P: return p();
    case SIMULATOR_Z0: return z();
    case SIMULATOR_RX0: return repmat(rx(), 1, ntout_);
    case SIMULATOR_RP: return repmat(rp(), 1, ntout_);
    case SIMULATOR_RZ0: return repmat(rz(), 1, ntout_);
    case SIMULATOR_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity Simulator::get_sparsity_out(casadi_int i) {
    switch (static_cast<SimulatorOutput>(i)) {
    case SIMULATOR_XF: return repmat(x(), 1, ntout_);
    case SIMULATOR_QF: return repmat(q(), 1, ntout_);
    case SIMULATOR_ZF: return repmat(z(), 1, ntout_);
    case SIMULATOR_RXF: return rx();
    case SIMULATOR_RQF: return rq();
    case SIMULATOR_RZF: return rz();
    case SIMULATOR_NUM_OUT: break;
    }
    return Sparsity();
  }

  Function Simulator::create_advanced(const Dict& opts) {
    return Function::create(this, opts);
  }

  int Simulator::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<SimulatorMemory*>(mem);

    // Read inputs
    const double* x0 = arg[SIMULATOR_X0];
    const double* z0 = arg[SIMULATOR_Z0];
    const double* p = arg[SIMULATOR_P];
    const double* rx0 = arg[SIMULATOR_RX0];
    const double* rz0 = arg[SIMULATOR_RZ0];
    const double* rp = arg[SIMULATOR_RP];
    arg += SIMULATOR_NUM_IN;

    // Read outputs
    double* x = res[SIMULATOR_XF];
    double* z = res[SIMULATOR_ZF];
    double* q = res[SIMULATOR_QF];
    double* rx = res[SIMULATOR_RXF];
    double* rz = res[SIMULATOR_RZF];
    double* rq = res[SIMULATOR_RQF];
    res += SIMULATOR_NUM_OUT;

    // Setup memory object
    setup(m, arg, res, iw, w);

    // Reset solver, take time to t0
    reset(m, grid_.front(), x0, z0, p);

    // Integrate forward
    for (casadi_int k=0; k<grid_.size(); ++k) {
      // Skip t0?
      if (k==0 && !output_t0_) continue;

      // Integrate forward
      advance(m, grid_[k], x, z, q);
      if (x) x += nx_;
      if (z) z += nz_;
      if (q) q += nq_;
    }

    // If backwards integration is needed
    if (nrx_>0) {
      // Integrate backward
      resetB(m, grid_.back(), rx0, rz0, rp);

      // Proceed to t0
      retreat(m, grid_.front(), rx, rz, rq);
    }

    if (print_stats_) print_stats(m);

    return 0;
  }

  const Options Simulator::options_
  = {{&OracleFunction::options_},
     {{"expand",
       {OT_BOOL,
        "Replace MX with SX expressions in problem formulation [false]"}},
      {"print_stats",
       {OT_BOOL,
        "Print out statistics after integration"}},
      {"t0",
       {OT_DOUBLE,
        "Beginning of the time horizon"}},
      {"tf",
       {OT_DOUBLE,
        "End of the time horizon"}},
      {"grid",
       {OT_DOUBLEVECTOR,
        "Time grid"}},
      {"augmented_options",
       {OT_DICT,
        "Options to be passed down to the augmented simulator, if one is constructed."}},
      {"output_t0",
       {OT_BOOL,
        "Output the state at the initial time"}}
     }
  };

  void Simulator::init(const Dict& opts) {
    // Default (temporary) options
    double t0=0, tf=1;
    bool expand = false;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="expand") {
        expand = op.second;
      } else if (op.first=="output_t0") {
        output_t0_ = op.second;
      } else if (op.first=="print_stats") {
        print_stats_ = op.second;
      } else if (op.first=="grid") {
        grid_ = op.second;
      } else if (op.first=="augmented_options") {
        augmented_options_ = op.second;
      } else if (op.first=="t0") {
        t0 = op.second;
      } else if (op.first=="tf") {
        tf = op.second;
      }
    }

    // Replace MX oracle with SX oracle?
    if (expand) this->expand();

    // Store a copy of the options, for creating augmented simulators
    opts_ = opts;

    // If grid unset, default to [t0, tf]
    if (grid_.empty()) {
      grid_ = {t0, tf};
    }

    ngrid_ = grid_.size();
    ntout_ = output_t0_ ? ngrid_ : ngrid_-1;

    // Call the base class method
    OracleFunction::init(opts);

    // For sparsity pattern propagation
    alloc(oracle_);

    // Error if sparse input
    casadi_assert(x().is_dense(), "Sparse DAE not supported");
    casadi_assert(z().is_dense(), "Sparse DAE not supported");
    casadi_assert(p().is_dense(), "Sparse DAE not supported");
    casadi_assert(rx().is_dense(), "Sparse DAE not supported");
    casadi_assert(rz().is_dense(), "Sparse DAE not supported");
    casadi_assert(rp().is_dense(), "Sparse DAE not supported");

    // Get dimensions (excluding sensitivity equations)
    nx1_ = x().size1();
    nz1_ = z().size1();
    nq1_ = q().size1();
    np1_  = p().size1();
    nrx1_ = rx().size1();
    nrz1_ = rz().size1();
    nrp1_ = rp().size1();
    nrq1_ = rq().size1();

    // Get dimensions (including sensitivity equations)
    nx_ = x().nnz();
    nz_ = z().nnz();
    nq_ = q().nnz();
    np_  = p().nnz();
    nrx_ = rx().nnz();
    nrz_ = rz().nnz();
    nrp_ = rp().nnz();
    nrq_ = rq().nnz();

    // Number of sensitivities
    ns_ = x().size2()-1;

    // Get the sparsities of the forward and reverse DAE
    sp_jac_dae_ = sp_jac_dae();
    casadi_assert(!sp_jac_dae_.is_singular(),
                          "Jacobian of the forward problem is structurally rank-deficient. "
                          "sprank(J)=" + str(sprank(sp_jac_dae_)) + "<"
                          + str(nx_+nz_));
    if (nrx_>0) {
      sp_jac_rdae_ = sp_jac_rdae();
      casadi_assert(!sp_jac_rdae_.is_singular(),
                            "Jacobian of the backward problem is structurally rank-deficient. "
                            "sprank(J)=" + str(sprank(sp_jac_rdae_)) + "<"
                            + str(nrx_+nrz_));
    }

    // Consistency check

    // Allocate sufficiently large work vectors
    alloc_w(nx_+nz_);
    alloc_w(nrx_+nrz_);
    alloc_w(nx_ + nz_ + nrx_ + nrz_, true);
  }

  int Simulator::init_mem(void* mem) const {
    if (OracleFunction::init_mem(mem)) return 1;

    //auto m = static_cast<SimulatorMemory*>(mem);
    return 0;
  }

  template<typename MatType>
  std::map<string, MatType> Simulator::aug_fwd(casadi_int nfwd) const {
    if (verbose_) casadi_message(name_ + "::aug_fwd");

    // Get input expressions
    vector<MatType> arg = MatType::get_input(oracle_);
    vector<MatType> aug_x, aug_z, aug_p, aug_rx, aug_rz, aug_rp;
    MatType aug_t = arg.at(DE2_T);
    aug_x.push_back(vec(arg.at(DE2_X)));
    aug_z.push_back(vec(arg.at(DE2_Z)));
    aug_p.push_back(vec(arg.at(DE2_P)));
    aug_rx.push_back(vec(arg.at(DE2_RX)));
    aug_rz.push_back(vec(arg.at(DE2_RZ)));
    aug_rp.push_back(vec(arg.at(DE2_RP)));

    // Get output expressions
    vector<MatType> res = oracle_(arg);
    vector<MatType> aug_ode, aug_alg, aug_quad, aug_rode, aug_ralg, aug_rquad;
    aug_ode.push_back(vec(res.at(DE2_ODE)));
    aug_alg.push_back(vec(res.at(DE2_ALG)));
    aug_quad.push_back(vec(res.at(DE2_QUAD)));
    aug_rode.push_back(vec(res.at(DE2_RODE)));
    aug_ralg.push_back(vec(res.at(DE2_RALG)));
    aug_rquad.push_back(vec(res.at(DE2_RQUAD)));

    // Zero of time dimension
    MatType zero_t = MatType::zeros(t());

    // Forward directional derivatives
    vector<vector<MatType>> seed(nfwd, vector<MatType>(DE2_NUM_IN));
    for (casadi_int d=0; d<nfwd; ++d) {
      seed[d][DE2_T] = zero_t;
      string pref = "aug" + str(d) + "_";
      aug_x.push_back(vec(seed[d][DE2_X] = MatType::sym(pref + "x", x())));
      aug_z.push_back(vec(seed[d][DE2_Z] = MatType::sym(pref + "z", z())));
      aug_p.push_back(vec(seed[d][DE2_P] = MatType::sym(pref + "p", p())));
      aug_rx.push_back(vec(seed[d][DE2_RX] = MatType::sym(pref + "rx", rx())));
      aug_rz.push_back(vec(seed[d][DE2_RZ] = MatType::sym(pref + "rz", rz())));
      aug_rp.push_back(vec(seed[d][DE2_RP] = MatType::sym(pref + "rp", rp())));
    }

    // Calculate directional derivatives
    vector<vector<MatType>> sens;
    oracle_->call_forward(arg, res, seed, sens, true, false);

    // Collect sensitivity equations
    casadi_assert_dev(sens.size()==nfwd);
    for (casadi_int d=0; d<nfwd; ++d) {
      casadi_assert_dev(sens[d].size()==DE2_NUM_OUT);
      aug_ode.push_back(vec(project(sens[d][DE2_ODE], x())));
      aug_alg.push_back(vec(project(sens[d][DE2_ALG], z())));
      aug_quad.push_back(vec(project(sens[d][DE2_QUAD], q())));
      aug_rode.push_back(vec(project(sens[d][DE2_RODE], rx())));
      aug_ralg.push_back(vec(project(sens[d][DE2_RALG], rz())));
      aug_rquad.push_back(vec(project(sens[d][DE2_RQUAD], rq())));
    }

    // Construct return object
    std::map<string, MatType> ret;
    ret["t"] = aug_t;
    ret["x"] = horzcat(aug_x);
    ret["z"] = horzcat(aug_z);
    ret["p"] = horzcat(aug_p);
    ret["ode"] = horzcat(aug_ode);
    ret["alg"] = horzcat(aug_alg);
    ret["quad"] = horzcat(aug_quad);
    ret["rx"] = horzcat(aug_rx);
    ret["rz"] = horzcat(aug_rz);
    ret["rp"] = horzcat(aug_rp);
    ret["rode"] = horzcat(aug_rode);
    ret["ralg"] = horzcat(aug_ralg);
    ret["rquad"] = horzcat(aug_rquad);
    return ret;
  }

  template<typename MatType>
  std::map<string, MatType> Simulator::aug_adj(casadi_int nadj) const {
    if (verbose_) casadi_message(name_ + "::aug_adj");

    // Get input expressions
    vector<MatType> arg = MatType::get_input(oracle_);
    vector<MatType> aug_x, aug_z, aug_p, aug_rx, aug_rz, aug_rp;
    MatType aug_t = arg.at(DE2_T);
    aug_x.push_back(vec(arg.at(DE2_X)));
    aug_z.push_back(vec(arg.at(DE2_Z)));
    aug_p.push_back(vec(arg.at(DE2_P)));
    aug_rx.push_back(vec(arg.at(DE2_RX)));
    aug_rz.push_back(vec(arg.at(DE2_RZ)));
    aug_rp.push_back(vec(arg.at(DE2_RP)));

    // Get output expressions
    vector<MatType> res = oracle_(arg);
    vector<MatType> aug_ode, aug_alg, aug_quad, aug_rode, aug_ralg, aug_rquad;
    aug_ode.push_back(vec(res.at(DE2_ODE)));
    aug_alg.push_back(vec(res.at(DE2_ALG)));
    aug_quad.push_back(vec(res.at(DE2_QUAD)));
    aug_rode.push_back(vec(res.at(DE2_RODE)));
    aug_ralg.push_back(vec(res.at(DE2_RALG)));
    aug_rquad.push_back(vec(res.at(DE2_RQUAD)));

    // Zero of time dimension
    MatType zero_t = MatType::zeros(t());

    // Reverse mode directional derivatives
    vector<vector<MatType>> seed(nadj, vector<MatType>(DE2_NUM_OUT));
    for (casadi_int d=0; d<nadj; ++d) {
      string pref = "aug" + str(d) + "_";
      aug_rx.push_back(vec(seed[d][DE2_ODE] = MatType::sym(pref + "ode", x())));
      aug_rz.push_back(vec(seed[d][DE2_ALG] = MatType::sym(pref + "alg", z())));
      aug_rp.push_back(vec(seed[d][DE2_QUAD] = MatType::sym(pref + "quad", q())));
      aug_x.push_back(vec(seed[d][DE2_RODE] = MatType::sym(pref + "rode", rx())));
      aug_z.push_back(vec(seed[d][DE2_RALG] = MatType::sym(pref + "ralg", rz())));
      aug_p.push_back(vec(seed[d][DE2_RQUAD] = MatType::sym(pref + "rquad", rq())));
    }

    // Calculate directional derivatives
    vector<vector<MatType>> sens;
    oracle_->call_reverse(arg, res, seed, sens, true, false);

    // Collect sensitivity equations
    casadi_assert_dev(sens.size()==nadj);
    for (casadi_int d=0; d<nadj; ++d) {
      casadi_assert_dev(sens[d].size()==DE2_NUM_IN);
      aug_rode.push_back(vec(project(sens[d][DE2_X], x())));
      aug_ralg.push_back(vec(project(sens[d][DE2_Z], z())));
      aug_rquad.push_back(vec(project(sens[d][DE2_P], p())));
      aug_ode.push_back(vec(project(sens[d][DE2_RX], rx())));
      aug_alg.push_back(vec(project(sens[d][DE2_RZ], rz())));
      aug_quad.push_back(vec(project(sens[d][DE2_RP], rp())));
    }

    // Construct return object
    std::map<string, MatType> ret;
    ret["t"] = aug_t;
    ret["x"] = vertcat(aug_x);
    ret["z"] = vertcat(aug_z);
    ret["p"] = vertcat(aug_p);
    ret["ode"] = vertcat(aug_ode);
    ret["alg"] = vertcat(aug_alg);
    ret["quad"] = vertcat(aug_quad);
    ret["rx"] = vertcat(aug_rx);
    ret["rz"] = vertcat(aug_rz);
    ret["rp"] = vertcat(aug_rp);
    ret["rode"] = vertcat(aug_rode);
    ret["ralg"] = vertcat(aug_ralg);
    ret["rquad"] = vertcat(aug_rquad);

    // Make sure that forward problem does not depend on backward states
    Function f("f", {ret["t"], ret["x"], ret["z"], ret["p"]},
                    {ret["ode"], ret["alg"], ret["quad"]});
    if (f.has_free()) {
      // Replace dependencies of rx, rz and rp with zeros
      f = Function("f", {ret["t"], ret["x"], ret["z"], ret["p"],
                         ret["rx"], ret["rz"], ret["rp"]},
                        {ret["ode"], ret["alg"], ret["quad"]});
      vector<MatType> v = {ret["t"], ret["x"], ret["z"], ret["p"], 0, 0, 0};
      v = f(v);
      ret["ode"] = v.at(0);
      ret["alg"] = v.at(1);
      ret["quad"] = v.at(2);
    }

    return ret;
  }

  int Simulator::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
    if (verbose_) casadi_message(name_ + "::sp_forward");

    // Work vectors
    bvec_t *tmp_x = w; w += nx_;
    bvec_t *tmp_z = w; w += nz_;
    bvec_t *tmp_rx = w; w += nrx_;
    bvec_t *tmp_rz = w; w += nrz_;

    // Propagate forward
    const bvec_t** arg1 = arg+n_in_;
    fill_n(arg1, static_cast<size_t>(DE2_NUM_IN), nullptr);
    arg1[DE2_X] = arg[SIMULATOR_X0];
    arg1[DE2_P] = arg[SIMULATOR_P];
    bvec_t** res1 = res+n_out_;
    fill_n(res1, static_cast<size_t>(DE2_NUM_OUT), nullptr);
    res1[DE2_ODE] = tmp_x;
    res1[DE2_ALG] = tmp_z;
    oracle_(arg1, res1, iw, w, 0);
    if (arg[SIMULATOR_X0]) {
      const bvec_t *tmp = arg[SIMULATOR_X0];
      for (casadi_int i=0; i<nx_; ++i) tmp_x[i] |= *tmp++;
    }

    // "Solve" in order to resolve interdependencies (cf. Rootfinder)
    copy_n(tmp_x, nx_+nz_, w);
    fill_n(tmp_x, nx_+nz_, 0);
    sp_jac_dae_.spsolve(tmp_x, w, false);

    // Get xf and zf
    if (res[SIMULATOR_XF]) copy_n(tmp_x, nx_, res[SIMULATOR_XF]);
    if (res[SIMULATOR_ZF]) copy_n(tmp_z, nz_, res[SIMULATOR_ZF]);

    // Propagate to quadratures
    if (nq_>0 && res[SIMULATOR_QF]) {
      arg1[DE2_X] = tmp_x;
      arg1[DE2_Z] = tmp_z;
      res1[DE2_ODE] = res1[DE2_ALG] = nullptr;
      res1[DE2_QUAD] = res[SIMULATOR_QF];
      if (oracle_(arg1, res1, iw, w, 0)) return 1;
    }

    if (nrx_>0) {
      // Propagate through g
      fill_n(arg1, static_cast<size_t>(DE2_NUM_IN), nullptr);
      arg1[DE2_X] = tmp_x;
      arg1[DE2_P] = arg[SIMULATOR_P];
      arg1[DE2_Z] = tmp_z;
      arg1[DE2_RX] = arg[SIMULATOR_X0];
      arg1[DE2_RX] = arg[SIMULATOR_RX0];
      arg1[DE2_RP] = arg[SIMULATOR_RP];
      fill_n(res1, static_cast<size_t>(DE2_NUM_OUT), nullptr);
      res1[DE2_RODE] = tmp_rx;
      res1[DE2_RALG] = tmp_rz;
      oracle_(arg1, res1, iw, w, 0);
      if (arg[SIMULATOR_RX0]) {
        const bvec_t *tmp = arg[SIMULATOR_RX0];
        for (casadi_int i=0; i<nrx_; ++i) tmp_rx[i] |= *tmp++;
      }

      // "Solve" in order to resolve interdependencies (cf. Rootfinder)
      copy_n(tmp_rx, nrx_+nrz_, w);
      fill_n(tmp_rx, nrx_+nrz_, 0);
      sp_jac_rdae_.spsolve(tmp_rx, w, false);

      // Get rxf and rzf
      if (res[SIMULATOR_RXF]) copy_n(tmp_rx, nrx_, res[SIMULATOR_RXF]);
      if (res[SIMULATOR_RZF]) copy_n(tmp_rz, nrz_, res[SIMULATOR_RZF]);

      // Propagate to quadratures
      if (nrq_>0 && res[SIMULATOR_RQF]) {
        arg1[DE2_RX] = tmp_rx;
        arg1[DE2_RZ] = tmp_rz;
        res1[DE2_RODE] = res1[DE2_RALG] = nullptr;
        res1[DE2_RQUAD] = res[SIMULATOR_RQF];
        if (oracle_(arg1, res1, iw, w, 0)) return 1;
      }
    }
    return 0;
  }

  int Simulator::sp_reverse(bvec_t** arg, bvec_t** res,
      casadi_int* iw, bvec_t* w, void* mem) const {
    if (verbose_) casadi_message(name_ + "::sp_reverse");

    // Work vectors
    bvec_t** arg1 = arg+n_in_;
    bvec_t** res1 = res+n_out_;
    bvec_t *tmp_x = w; w += nx_;
    bvec_t *tmp_z = w; w += nz_;

    // Shorthands
    bvec_t* x0 = arg[SIMULATOR_X0];
    bvec_t* p = arg[SIMULATOR_P];
    bvec_t* xf = res[SIMULATOR_XF];
    bvec_t* zf = res[SIMULATOR_ZF];
    bvec_t* qf = res[SIMULATOR_QF];

    // Propagate from outputs to state vectors
    if (xf) {
      copy_n(xf, nx_, tmp_x);
      fill_n(xf, nx_, 0);
    } else {
      fill_n(tmp_x, nx_, 0);
    }
    if (zf) {
      copy_n(zf, nz_, tmp_z);
      fill_n(zf, nz_, 0);
    } else {
      fill_n(tmp_z, nz_, 0);
    }

    if (nrx_>0) {
      // Work vectors
      bvec_t *tmp_rx = w; w += nrx_;
      bvec_t *tmp_rz = w; w += nrz_;

      // Shorthands
      bvec_t* rx0 = arg[SIMULATOR_RX0];
      bvec_t* rp = arg[SIMULATOR_RP];
      bvec_t* rxf = res[SIMULATOR_RXF];
      bvec_t* rzf = res[SIMULATOR_RZF];
      bvec_t* rqf = res[SIMULATOR_RQF];

      // Propagate from outputs to state vectors
      if (rxf) {
        copy_n(rxf, nrx_, tmp_rx);
        fill_n(rxf, nrx_, 0);
      } else {
        fill_n(tmp_rx, nrx_, 0);
      }
      if (rzf) {
        copy_n(rzf, nrz_, tmp_rz);
        fill_n(rzf, nrz_, 0);
      } else {
        fill_n(tmp_rz, nrz_, 0);
      }

      // Get dependencies from backward quadratures
      fill_n(res1, static_cast<size_t>(DE2_NUM_OUT), nullptr);
      fill_n(arg1, static_cast<size_t>(DE2_NUM_IN), nullptr);
      res1[DE2_RQUAD] = rqf;
      arg1[DE2_X] = tmp_x;
      arg1[DE2_Z] = tmp_z;
      arg1[DE2_P] = p;
      arg1[DE2_RX] = tmp_rx;
      arg1[DE2_RZ] = tmp_rz;
      arg1[DE2_RP] = rp;
      if (oracle_.rev(arg1, res1, iw, w, 0)) return 1;

      // Propagate interdependencies
      fill_n(w, nrx_+nrz_, 0);
      sp_jac_rdae_.spsolve(w, tmp_rx, true);
      copy_n(w, nrx_+nrz_, tmp_rx);

      // Direct dependency rx0 -> rxf
      if (rx0) for (casadi_int i=0; i<nrx_; ++i) rx0[i] |= tmp_rx[i];

      // Indirect dependency via g
      res1[DE2_RODE] = tmp_rx;
      res1[DE2_RALG] = tmp_rz;
      res1[DE2_RQUAD] = nullptr;
      arg1[DE2_RX] = rx0;
      arg1[DE2_RZ] = nullptr; // arg[SIMULATOR_RZ0] is a guess, no dependency
      if (oracle_.rev(arg1, res1, iw, w, 0)) return 1;
    }

    // Get dependencies from forward quadratures
    fill_n(res1, static_cast<size_t>(DE2_NUM_OUT), nullptr);
    fill_n(arg1, static_cast<size_t>(DE2_NUM_IN), nullptr);
    res1[DE2_QUAD] = qf;
    arg1[DE2_X] = tmp_x;
    arg1[DE2_Z] = tmp_z;
    arg1[DE2_P] = p;
    if (qf && nq_>0) {
      if (oracle_.rev(arg1, res1, iw, w, 0)) return 1;
    }

    // Propagate interdependencies
    fill_n(w, nx_+nz_, 0);
    sp_jac_dae_.spsolve(w, tmp_x, true);
    copy_n(w, nx_+nz_, tmp_x);

    // Direct dependency x0 -> xf
    if (x0) for (casadi_int i=0; i<nx_; ++i) x0[i] |= tmp_x[i];

    // Indirect dependency through f
    res1[DE2_ODE] = tmp_x;
    res1[DE2_ALG] = tmp_z;
    res1[DE2_QUAD] = nullptr;
    arg1[DE2_X] = x0;
    arg1[DE2_Z] = nullptr; // arg[SIMULATOR_Z0] is a guess, no dependency
    if (oracle_.rev(arg1, res1, iw, w, 0)) return 1;
    return 0;
  }

  Function Simulator::
  get_forward(casadi_int nfwd, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {
    if (verbose_) casadi_message(name_ + "::get_forward");

    // Simulator options
    Dict aug_opts = getDerivativeOptions(true);
    for (auto&& i : augmented_options_) {
      aug_opts[i.first] = i.second;
    }

    // Create simulator for augmented DAE
    Function aug_dae;
    string aug_prefix = "fsens" + str(nfwd) + "_";
    string dae_name = aug_prefix + oracle_.name();
    Dict dae_opts = {{"derivative_of", oracle_}};
    if (oracle_.is_a("SXFunction")) {
      aug_dae = map2oracle(dae_name, aug_fwd<SX>(nfwd));
    } else {
      aug_dae = map2oracle(dae_name, aug_fwd<MX>(nfwd));
    }
    aug_opts["derivative_of"] = self();
    Function aug_int = simulator(aug_prefix + name_, plugin_name(),
      aug_dae, aug_opts);

    // All inputs of the return function
    vector<MX> ret_in;
    ret_in.reserve(SIMULATOR_NUM_IN*(1+nfwd) + SIMULATOR_NUM_OUT);

    // Augmented state
    vector<MX> x0_aug, p_aug, z0_aug, rx0_aug, rp_aug, rz0_aug;

    // Add nondifferentiated inputs and forward seeds
    for (casadi_int dir=-1; dir<nfwd; ++dir) {
      // Suffix
      string suff;
      if (dir>=0) suff = "_" + str(dir);

      // Augmented problem
      vector<MX> din(SIMULATOR_NUM_IN);
      x0_aug.push_back(vec(din[SIMULATOR_X0] = MX::sym("x0" + suff, x())));
      p_aug.push_back(vec(din[SIMULATOR_P] = MX::sym("p" + suff, p())));
      z0_aug.push_back(vec(din[SIMULATOR_Z0] = MX::sym("z0" + suff, z())));
      rx0_aug.push_back(vec(din[SIMULATOR_RX0] = MX::sym("rx0" + suff, rx())));
      rp_aug.push_back(vec(din[SIMULATOR_RP] = MX::sym("rp" + suff, rp())));
      rz0_aug.push_back(vec(din[SIMULATOR_RZ0] = MX::sym("rz0" + suff, rz())));
      ret_in.insert(ret_in.end(), din.begin(), din.end());

      // Dummy outputs
      if (dir==-1) {
        vector<MX> dout(SIMULATOR_NUM_OUT);
        dout[SIMULATOR_XF]  = MX::sym("xf_dummy", Sparsity(size_out(SIMULATOR_XF)));
        dout[SIMULATOR_QF]  = MX::sym("qf_dummy", Sparsity(q().size()));
        dout[SIMULATOR_ZF]  = MX::sym("zf_dummy", Sparsity(z().size()));
        dout[SIMULATOR_RXF]  = MX::sym("rxf_dummy", Sparsity(rx().size()));
        dout[SIMULATOR_RQF]  = MX::sym("rqf_dummy", Sparsity(rq().size()));
        dout[SIMULATOR_RZF]  = MX::sym("rzf_dummy", Sparsity(rz().size()));
        ret_in.insert(ret_in.end(), dout.begin(), dout.end());
      }
    }

    // Call the simulator
    vector<MX> simulator_in(SIMULATOR_NUM_IN);
    simulator_in[SIMULATOR_X0] = horzcat(x0_aug);
    simulator_in[SIMULATOR_P] = horzcat(p_aug);
    simulator_in[SIMULATOR_Z0] = horzcat(z0_aug);
    simulator_in[SIMULATOR_RX0] = horzcat(rx0_aug);
    simulator_in[SIMULATOR_RP] = horzcat(rp_aug);
    simulator_in[SIMULATOR_RZ0] = horzcat(rz0_aug);
    vector<MX> simulator_out = aug_int(simulator_in);
    for (auto&& e : simulator_out) {
      // Workaround
      if (e.size2()!=1+nfwd) e = reshape(e, -1, 1+nfwd);
    }

    // Augmented results
    vector<casadi_int> offset = range(1+nfwd+1);
    vector<MX> xf_aug = horzsplit(simulator_out[SIMULATOR_XF], offset);
    vector<MX> qf_aug = horzsplit(simulator_out[SIMULATOR_QF], offset);
    vector<MX> zf_aug = horzsplit(simulator_out[SIMULATOR_ZF], offset);
    vector<MX> rxf_aug = horzsplit(simulator_out[SIMULATOR_RXF], offset);
    vector<MX> rqf_aug = horzsplit(simulator_out[SIMULATOR_RQF], offset);
    vector<MX> rzf_aug = horzsplit(simulator_out[SIMULATOR_RZF], offset);

    // All outputs of the return function
    vector<MX> ret_out;
    ret_out.reserve(SIMULATOR_NUM_OUT*nfwd);

    // Collect the forward sensitivities
    vector<MX> dd(SIMULATOR_NUM_IN);
    for (casadi_int dir=0; dir<nfwd; ++dir) {
      dd[SIMULATOR_XF]  = reshape(xf_aug.at(dir+1), x().size());
      dd[SIMULATOR_QF]  = reshape(qf_aug.at(dir+1), q().size());
      dd[SIMULATOR_ZF]  = reshape(zf_aug.at(dir+1), z().size());
      dd[SIMULATOR_RXF] = reshape(rxf_aug.at(dir+1), rx().size());
      dd[SIMULATOR_RQF] = reshape(rqf_aug.at(dir+1), rq().size());
      dd[SIMULATOR_RZF] = reshape(rzf_aug.at(dir+1), rz().size());
      ret_out.insert(ret_out.end(), dd.begin(), dd.end());
    }

    // Concatenate forward seeds
    vector<MX> v(nfwd);
    auto r_it = ret_in.begin() + n_in_ + n_out_;
    for (casadi_int i=0; i<n_in_; ++i) {
      for (casadi_int d=0; d<nfwd; ++d) v[d] = *(r_it + d*n_in_);
      *r_it++ = horzcat(v);
    }
    ret_in.resize(n_in_ + n_out_ + n_in_);

    // Concatenate forward sensitivites
    r_it = ret_out.begin();
    for (casadi_int i=0; i<n_out_; ++i) {
      for (casadi_int d=0; d<nfwd; ++d) v[d] = *(r_it + d*n_out_);
      *r_it++ = horzcat(v);
    }
    ret_out.resize(n_out_);

    // Create derivative function and return
    return Function(name, ret_in, ret_out, inames, onames, opts);
  }

  Function Simulator::
  get_reverse(casadi_int nadj, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {
    if (verbose_) casadi_message(name_ + "::get_reverse");

    // Simulator options
    Dict aug_opts = getDerivativeOptions(false);
    for (auto&& i : augmented_options_) {
      aug_opts[i.first] = i.second;
    }

    // Create simulator for augmented DAE
    Function aug_dae;
    string aug_prefix = "asens" + str(nadj) + "_";
    string dae_name = aug_prefix + oracle_.name();
    if (oracle_.is_a("SXFunction")) {
      aug_dae = map2oracle(dae_name, aug_adj<SX>(nadj));
    } else {
      aug_dae = map2oracle(dae_name, aug_adj<MX>(nadj));
    }
    aug_opts["derivative_of"] = self();
    Function aug_int = simulator(aug_prefix + name_, plugin_name(),
      aug_dae, aug_opts);

    // All inputs of the return function
    vector<MX> ret_in;
    ret_in.reserve(SIMULATOR_NUM_IN + SIMULATOR_NUM_OUT*(1+nadj));

    // Augmented state
    vector<MX> x0_aug, p_aug, z0_aug, rx0_aug, rp_aug, rz0_aug;

    // Inputs or forward/adjoint seeds in one direction
    vector<MX> dd(SIMULATOR_NUM_IN);
    x0_aug.push_back(vec(dd[SIMULATOR_X0] = MX::sym("x0", x())));
    p_aug.push_back(vec(dd[SIMULATOR_P] = MX::sym("p", p())));
    z0_aug.push_back(vec(dd[SIMULATOR_Z0] = MX::sym("r0", z())));
    rx0_aug.push_back(vec(dd[SIMULATOR_RX0] = MX::sym("rx0", rx())));
    rp_aug.push_back(vec(dd[SIMULATOR_RP] = MX::sym("rp", rp())));
    rz0_aug.push_back(vec(dd[SIMULATOR_RZ0] = MX::sym("rz0", rz())));
    ret_in.insert(ret_in.end(), dd.begin(), dd.end());

    // Add dummy inputs (outputs of the nondifferentiated funciton)
    dd.resize(SIMULATOR_NUM_OUT);
    fill(dd.begin(), dd.end(), MX());
    dd[SIMULATOR_XF]  = MX::sym("xf_dummy", Sparsity(x().size()));
    dd[SIMULATOR_QF]  = MX::sym("qf_dummy", Sparsity(q().size()));
    dd[SIMULATOR_ZF]  = MX::sym("zf_dummy", Sparsity(z().size()));
    dd[SIMULATOR_RXF]  = MX::sym("rxf_dummy", Sparsity(rx().size()));
    dd[SIMULATOR_RQF]  = MX::sym("rqf_dummy", Sparsity(rq().size()));
    dd[SIMULATOR_RZF]  = MX::sym("rzf_dummy", Sparsity(rz().size()));
    ret_in.insert(ret_in.end(), dd.begin(), dd.end());

    // Add adjoint seeds
    dd.resize(SIMULATOR_NUM_OUT);
    fill(dd.begin(), dd.end(), MX());
    for (casadi_int dir=0; dir<nadj; ++dir) {
      // Suffix
      string suff;
      if (dir>=0) suff = "_" + str(dir);

      // Augmented problem
      rx0_aug.push_back(vec(dd[SIMULATOR_XF] = MX::sym("xf" + suff, x())));
      rp_aug.push_back(vec(dd[SIMULATOR_QF] = MX::sym("qf" + suff, q())));
      rz0_aug.push_back(vec(dd[SIMULATOR_ZF] = MX::sym("zf" + suff, z())));
      x0_aug.push_back(vec(dd[SIMULATOR_RXF] = MX::sym("rxf" + suff, rx())));
      p_aug.push_back(vec(dd[SIMULATOR_RQF] = MX::sym("rqf" + suff, rq())));
      z0_aug.push_back(vec(dd[SIMULATOR_RZF] = MX::sym("rzf" + suff, rz())));
      ret_in.insert(ret_in.end(), dd.begin(), dd.end());
    }

    // Call the simulator
    vector<MX> simulator_in(SIMULATOR_NUM_IN);
    simulator_in[SIMULATOR_X0] = vertcat(x0_aug);
    simulator_in[SIMULATOR_P] = vertcat(p_aug);
    simulator_in[SIMULATOR_Z0] = vertcat(z0_aug);
    simulator_in[SIMULATOR_RX0] = vertcat(rx0_aug);
    simulator_in[SIMULATOR_RP] = vertcat(rp_aug);
    simulator_in[SIMULATOR_RZ0] = vertcat(rz0_aug);
    vector<MX> simulator_out = aug_int(simulator_in);

    // Get offset in the splitted problem
    vector<casadi_int> off_x = {0, x().numel()};
    vector<casadi_int> off_z = {0, z().numel()};
    vector<casadi_int> off_q = {0, q().numel()};
    vector<casadi_int> off_p = {0, p().numel()};
    vector<casadi_int> off_rx = {0, rx().numel()};
    vector<casadi_int> off_rz = {0, rz().numel()};
    vector<casadi_int> off_rq = {0, rq().numel()};
    vector<casadi_int> off_rp = {0, rp().numel()};
    for (casadi_int dir=0; dir<nadj; ++dir) {
      off_x.push_back(off_x.back() + rx().numel());
      off_z.push_back(off_z.back() + rz().numel());
      off_q.push_back(off_q.back() + rp().numel());
      off_p.push_back(off_p.back() + rq().numel());
      off_rx.push_back(off_rx.back() + x().numel());
      off_rz.push_back(off_rz.back() + z().numel());
      off_rq.push_back(off_rq.back() + p().numel());
      off_rp.push_back(off_rp.back() + q().numel());
    }

    // Augmented results
    vector<MX> xf_aug = vertsplit(simulator_out[SIMULATOR_XF], off_x);
    vector<MX> qf_aug = vertsplit(simulator_out[SIMULATOR_QF], off_q);
    vector<MX> zf_aug = vertsplit(simulator_out[SIMULATOR_ZF], off_z);
    vector<MX> rxf_aug = vertsplit(simulator_out[SIMULATOR_RXF], off_rx);
    vector<MX> rqf_aug = vertsplit(simulator_out[SIMULATOR_RQF], off_rq);
    vector<MX> rzf_aug = vertsplit(simulator_out[SIMULATOR_RZF], off_rz);

    // All outputs of the return function
    vector<MX> ret_out;
    ret_out.reserve(SIMULATOR_NUM_IN*nadj);

    // Collect the adjoint sensitivities
    dd.resize(SIMULATOR_NUM_IN);
    fill(dd.begin(), dd.end(), MX());
    for (casadi_int dir=0; dir<nadj; ++dir) {
      dd[SIMULATOR_X0]  = reshape(rxf_aug.at(dir+1), x().size());
      dd[SIMULATOR_P]   = reshape(rqf_aug.at(dir+1), p().size());
      dd[SIMULATOR_Z0]  = reshape(rzf_aug.at(dir+1), z().size());
      dd[SIMULATOR_RX0] = reshape(xf_aug.at(dir+1), rx().size());
      dd[SIMULATOR_RP]  = reshape(qf_aug.at(dir+1), rp().size());
      dd[SIMULATOR_RZ0] = reshape(zf_aug.at(dir+1), rz().size());
      ret_out.insert(ret_out.end(), dd.begin(), dd.end());
    }

    // Concatenate forward seeds
    vector<MX> v(nadj);
    auto r_it = ret_in.begin() + n_in_ + n_out_;
    for (casadi_int i=0; i<n_out_; ++i) {
      for (casadi_int d=0; d<nadj; ++d) v[d] = *(r_it + d*n_out_);
      *r_it++ = horzcat(v);
    }
    ret_in.resize(n_in_ + n_out_ + n_out_);

    // Concatenate forward sensitivites
    r_it = ret_out.begin();
    for (casadi_int i=0; i<n_in_; ++i) {
      for (casadi_int d=0; d<nadj; ++d) v[d] = *(r_it + d*n_in_);
      *r_it++ = horzcat(v);
    }
    ret_out.resize(n_in_);

    // Create derivative function and return
    return Function(name, ret_in, ret_out, inames, onames, opts);
  }

  Dict Simulator::getDerivativeOptions(bool fwd) const {
    // Copy all options
    return opts_;
  }

  Sparsity Simulator::sp_jac_dae() {
    // Start with the sparsity pattern of the ODE part
    Sparsity jac_ode_x = oracle_.jac_sparsity(DE2_ODE, DE2_X);

    // Add diagonal to get interdependencies
    jac_ode_x = jac_ode_x + Sparsity::diag(nx_);

    // Quick return if no algebraic variables
    if (nz_==0) return jac_ode_x;

    // Add contribution from algebraic variables and equations
    Sparsity jac_ode_z = oracle_.jac_sparsity(DE2_ODE, DE2_Z);
    Sparsity jac_alg_x = oracle_.jac_sparsity(DE2_ALG, DE2_X);
    Sparsity jac_alg_z = oracle_.jac_sparsity(DE2_ALG, DE2_Z);
    return blockcat(jac_ode_x, jac_ode_z,
                    jac_alg_x, jac_alg_z);
  }

  Sparsity Simulator::sp_jac_rdae() {
    // Start with the sparsity pattern of the ODE part
    Sparsity jac_ode_x = oracle_.jac_sparsity(DE2_RODE, DE2_RX);

    // Add diagonal to get interdependencies
    jac_ode_x = jac_ode_x + Sparsity::diag(nrx_);

    // Quick return if no algebraic variables
    if (nrz_==0) return jac_ode_x;

    // Add contribution from algebraic variables and equations
    Sparsity jac_ode_z = oracle_.jac_sparsity(DE2_RODE, DE2_RZ);
    Sparsity jac_alg_x = oracle_.jac_sparsity(DE2_RALG, DE2_RX);
    Sparsity jac_alg_z = oracle_.jac_sparsity(DE2_RALG, DE2_RZ);
    return blockcat(jac_ode_x, jac_ode_z,
                    jac_alg_x, jac_alg_z);
  }

  std::map<std::string, Simulator::Plugin> Simulator::solvers_;

  const std::string Simulator::infix_ = "simulator";

  void Simulator::setStopTime(SimulatorMemory* mem, double tf) const {
    casadi_error("setStopTime not defined for class " + class_name());
  }

  template<typename XType>
  Function Simulator::map2oracle(const std::string& name,
    const std::map<std::string, XType>& d, const Dict& opts) {
    std::vector<XType> de_in(DE2_NUM_IN), de_out(DE2_NUM_OUT);

    for (auto&& i : d) {
      if (i.first=="t") {
        de_in[DE2_T]=i.second;
      } else if (i.first=="x") {
        de_in[DE2_X]=i.second;
      } else if (i.first=="z") {
        de_in[DE2_Z]=i.second;
      } else if (i.first=="p") {
        de_in[DE2_P]=i.second;
      } else if (i.first=="rx") {
        de_in[DE2_RX]=i.second;
      } else if (i.first=="rz") {
        de_in[DE2_RZ]=i.second;
      } else if (i.first=="rp") {
        de_in[DE2_RP]=i.second;
      } else if (i.first=="ode") {
        de_out[DE2_ODE]=i.second;
      } else if (i.first=="alg") {
        de_out[DE2_ALG]=i.second;
      } else if (i.first=="quad") {
        de_out[DE2_QUAD]=i.second;
      } else if (i.first=="rode") {
        de_out[DE2_RODE]=i.second;
      } else if (i.first=="ralg") {
        de_out[DE2_RALG]=i.second;
      } else if (i.first=="rquad") {
        de_out[DE2_RQUAD]=i.second;
      } else {
        casadi_error("No such field: " + i.first);
      }
    }

    // Make sure x and ode exist
    casadi_assert(!de_in[DE2_X].is_empty(), "Ill-posed ODE - no state");

    // Number of right-hand-sides
    casadi_int nrhs = de_in[DE2_X].size2();

    // Make sure consistent number of right-hand-sides
    for (bool b : {true, false}) {
      for (auto&& e : b ? de_in : de_out) {
        // Skip time
        if (&e == &de_in[DE2_T]) continue;
        // Number of rows
        casadi_int nr = e.size1();
        // Make sure no change in number of elements
        casadi_assert(e.numel()==nr*nrhs, "Inconsistent number of rhs");
        e = reshape(e, nr, nrhs);
      }
    }

    // Consistent sparsity for x
    casadi_assert(de_in[DE2_X].size()==de_out[DE2_ODE].size(),
      "Dimension mismatch for 'ode'");
    de_out[DE2_ODE] = project(de_out[DE2_ODE], de_in[DE2_X].sparsity());

    // Consistent sparsity for z
    casadi_assert(de_in[DE2_Z].size()==de_out[DE2_ALG].size(),
      "Dimension mismatch for 'alg'");
    de_out[DE2_ALG] = project(de_out[DE2_ALG], de_in[DE2_Z].sparsity());

    // Consistent sparsity for rx
    casadi_assert(de_in[DE2_RX].size()==de_out[DE2_RODE].size(),
      "Dimension mismatch for 'rode'");
    de_out[DE2_RODE] = project(de_out[DE2_RODE], de_in[DE2_RX].sparsity());

    // Consistent sparsity for rz
    casadi_assert(de_in[DE2_RZ].size()==de_out[DE2_RALG].size(),
      "Dimension mismatch for 'ralg'");
    de_out[DE2_RALG] = project(de_out[DE2_RALG], de_in[DE2_RZ].sparsity());

    // Construct
    return Function(name, de_in, de_out, DE2_INPUTS, DE2_OUTPUTS, opts);
  }

  void Simulator::serialize_body(SerializingStream &s) const {
    OracleFunction::serialize_body(s);

    s.version("Simulator", 1);
    s.pack("Simulator::sp_jac_dae", sp_jac_dae_);
    s.pack("Simulator::sp_jac_rdae", sp_jac_rdae_);
    s.pack("Simulator::nx", nx_);
    s.pack("Simulator::nz", nz_);
    s.pack("Simulator::nq", nq_);
    s.pack("Simulator::nx1", nx1_);
    s.pack("Simulator::nz1", nz1_);
    s.pack("Simulator::nq1", nq1_);
    s.pack("Simulator::nrx", nrx_);
    s.pack("Simulator::nrz", nrz_);
    s.pack("Simulator::nrq", nrq_);
    s.pack("Simulator::nrx1", nrx1_);
    s.pack("Simulator::nrz1", nrz1_);
    s.pack("Simulator::nrq1", nrq1_);
    s.pack("Simulator::np", np_);
    s.pack("Simulator::nrp", nrp_);
    s.pack("Simulator::np1", np1_);
    s.pack("Simulator::nrp1", nrp1_);
    s.pack("Simulator::ns", ns_);
    s.pack("Simulator::grid", grid_);
    s.pack("Simulator::ngrid", ngrid_);
    s.pack("Simulator::augmented_options", augmented_options_);
    s.pack("Simulator::opts", opts_);
    s.pack("Simulator::onestep", onestep_);
    s.pack("Simulator::print_stats", print_stats_);
    s.pack("Simulator::output_t0", output_t0_);
    s.pack("Simulator::ntout", ntout_);
  }

  void Simulator::serialize_type(SerializingStream &s) const {
    OracleFunction::serialize_type(s);
    PluginInterface<Simulator>::serialize_type(s);
  }

  ProtoFunction* Simulator::deserialize(DeserializingStream& s) {
    return PluginInterface<Simulator>::deserialize(s);
  }

  Simulator::Simulator(DeserializingStream & s) : OracleFunction(s) {
    s.version("Simulator", 1);
    s.unpack("Simulator::sp_jac_dae", sp_jac_dae_);
    s.unpack("Simulator::sp_jac_rdae", sp_jac_rdae_);
    s.unpack("Simulator::nx", nx_);
    s.unpack("Simulator::nz", nz_);
    s.unpack("Simulator::nq", nq_);
    s.unpack("Simulator::nx1", nx1_);
    s.unpack("Simulator::nz1", nz1_);
    s.unpack("Simulator::nq1", nq1_);
    s.unpack("Simulator::nrx", nrx_);
    s.unpack("Simulator::nrz", nrz_);
    s.unpack("Simulator::nrq", nrq_);
    s.unpack("Simulator::nrx1", nrx1_);
    s.unpack("Simulator::nrz1", nrz1_);
    s.unpack("Simulator::nrq1", nrq1_);
    s.unpack("Simulator::np", np_);
    s.unpack("Simulator::nrp", nrp_);
    s.unpack("Simulator::np1", np1_);
    s.unpack("Simulator::nrp1", nrp1_);
    s.unpack("Simulator::ns", ns_);
    s.unpack("Simulator::grid", grid_);
    s.unpack("Simulator::ngrid", ngrid_);
    s.unpack("Simulator::augmented_options", augmented_options_);
    s.unpack("Simulator::opts", opts_);
    s.unpack("Simulator::onestep", onestep_);
    s.unpack("Simulator::print_stats", print_stats_);
    s.unpack("Simulator::output_t0", output_t0_);
    s.unpack("Simulator::ntout", ntout_);
  }

} // namespace casadi
