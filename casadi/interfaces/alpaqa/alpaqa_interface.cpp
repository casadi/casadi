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


#include "alpaqa_interface.hpp"
#include <vector>

namespace casadi {

  extern "C"
  int CASADI_NLPSOL_ALPAQA_EXPORT
  casadi_register_nlpsol_alpaqa(Nlpsol::Plugin* plugin) {
    plugin->creator = AlpaqaInterface::creator;
    plugin->name = "alpaqa";
    plugin->doc = AlpaqaInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &AlpaqaInterface::options_;
    plugin->deserialize = &AlpaqaInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_ALPAQA_EXPORT casadi_load_nlpsol_alpaqa() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_alpaqa);
  }

  AlpaqaInterface::AlpaqaInterface(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }

  AlpaqaInterface::~AlpaqaInterface() {
    clear_mem();
  }

  const Options AlpaqaInterface::options_
  = {{&Nlpsol::options_},
      {{"alpaqa",
        {OT_DICT,
        "Options to be passed to Alpaqa"}},
       {"print_level",
        {OT_INT,
        "Print level of SLEQP (default: 3)"}},
       {"max_iter",
        {OT_INT,
        "Maximum number of iterations"}},
       {"max_wall_time",
        {OT_DOUBLE,
        "maximum wall time allowed"}}
      }
  };

  const std::string AlpaqaInterface::meta_doc = "";


  void AlpaqaInterface::init(const Dict& opts) {
    Nlpsol::init(opts);
    MX x = MX::sym("x", nx_);
    MX p = MX::sym("p", np_);
    MX y = MX::sym("y", ng_);
    MX s = MX::sym("s");
    MX v = MX::sym("v", nx_);
    MX sigma = MX::sym("sigma", ng_);
    MX zl = MX::sym("zl", ng_);
    MX zu = MX::sym("zu", ng_);

    std::vector<MX> ret = oracle_(std::vector<MX>{x, p});
    MX f = ret[0];
    MX g = ret[1];

    MX sL = s * f + dot(y, g);
    MX L = f + dot(y, g);
    MX zeta = g + (y / sigma);
    MX z_hat = fmax(zl, fmin(zeta, zu));
    MX d = zeta - z_hat;
    MX y_hat = sigma * d;
    MX s_psi = s * f + 0.5 * dot(y_hat, d);
    MX psi = f + 0.5 * dot(y_hat, d);

    // Setup NLP functions
    create_function("nlp_f", {"x", "p"}, {"f"});
    create_function("nlp_grad_f", {"x", "p"}, {"f", "grad:f:x"});
    create_function("nlp_g", {"x", "p"}, {"g"});
    create_function("nlp_jac_g", {"x", "p"}, {"jac:g:x"});

    create_function("nlp_psi", {x, p, y, sigma, zl, zu}, {psi, y_hat},
      {"x", "p", "y", "sigma", "zl", "zu"}, {"psi", "y_hat"});

    MX grad_psi = gradient(psi, x);
    create_function("nlp_grad_psi", {x, p, y, sigma, zl, zu}, {psi, grad_psi},
      {"x", "p", "y", "sigma", "zl", "zu"}, {"psi", "grad_psi"});

    MX grad_L = gradient(L, x);
    create_function("nlp_grad_L", {x, p, y}, {grad_L},
      {"x", "p", "y"}, {"grad_L"});

    MX hess_L = hessian(sL, x)(0);
    if (!hess_L.is_dense()) hess_L = triu(hess_L);
    create_function("nlp_hess_L", {x, p, y, s}, {hess_L},
      {"x", "p", "y", "s"}, {"hess_L"});

    MX L_prod = gradient(jtimes(sL, x, v, false), x);
    create_function("nlp_hess_L_prod", {x, p, y, s, v}, {L_prod},
      {"x", "p", "y", "s", "v"}, {"L_prod"});

    MX hess_psi = hessian(sL, x)(0);
    if (!hess_psi.is_dense()) hess_psi = triu(hess_psi);
    create_function("nlp_hess_psi", {x, p, y, sigma, s, zl, zu}, {hess_psi},
      {"x", "p", "y", "sigma", "s", "zl", "zu"}, {"hess_psi"});

    MX hess_psi_prod = gradient(jtimes(s_psi, x, v, false), x);
    create_function("nlp_hess_psi_prod", {x, p, y, sigma, s, zl, zu, v}, {hess_psi},
      {"x", "p", "y", "sigma", "s", "zl", "zu", "v"}, {"hess_psi_prod"});
  }

  /** \brief Initalize memory block */
  int AlpaqaInterface::init_mem(void* mem) const {
    if (Nlpsol::init_mem(mem)) return 1;

    AlpaqaMemory* m = static_cast<AlpaqaMemory*>(mem);


    return 0;
  }

  void AlpaqaInterface::clear_mem_at(AlpaqaMemory* m) const {

  }

  void AlpaqaInterface::free_mem(void *mem) const {
    AlpaqaMemory* m = static_cast<AlpaqaMemory*>(mem);


    delete m;
  }


  /// Get all statistics
  Dict AlpaqaInterface::get_stats(void* mem) const {
    Dict stats = Nlpsol::get_stats(mem);

    AlpaqaMemory* m = static_cast<AlpaqaMemory*>(mem);

    return stats;
  }

  /** \brief Set the (persistent) work vectors */
  void AlpaqaInterface::set_work(void* mem, const double**& arg, double**& res,
                                casadi_int*& iw, double*& w) const {

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    AlpaqaMemory* m = static_cast<AlpaqaMemory*>(mem);

    check_inputs(m);

    // clear_mem(m);

    casadi_nlpsol_data<double>& d_nlp = m->d_nlp;


    return;
  }

  // Solve the NLP
  int AlpaqaInterface::solve(void* mem) const {
    AlpaqaMemory* m = static_cast<AlpaqaMemory*>(mem);
    auto d_nlp = &m->d_nlp;
    AlpaqaProblem problem(*this, m);

    //problem.l1_reg = 0;
    //problem.penalty_alm_split = 0;
    casadi_copy(d_nlp->lbz, nx_, problem.C.lowerbound.data());
    casadi_copy(d_nlp->ubz, nx_, problem.C.upperbound.data());
    casadi_copy(d_nlp->lbz+nx_, ng_, problem.D.lowerbound.data());
    casadi_copy(d_nlp->ubz+nx_, ng_, problem.D.upperbound.data());

    // Wrap the problem to count the function evaluations
    auto counted_problem = alpaqa::problem_with_counters_ref(problem);

    // Define the solvers to use
    using Direction   = alpaqa::LBFGSDirection<alpaqa::DefaultConfig>;
    using InnerSolver = alpaqa::PANOCSolver<Direction>;
    using OuterSolver = alpaqa::ALMSolver<InnerSolver>;

    // Settings for the outer augmented Lagrangian method
    OuterSolver::Params almparam;
    almparam.tolerance             = 1e-8; // tolerance
    almparam.dual_tolerance        = 1e-8;
    almparam.penalty_update_factor = 10; // penalty update factor
    almparam.max_iter              = 20;
    almparam.print_interval        = 1;

    // Settings for the inner PANOC solver
    InnerSolver::Params panocparam;
    panocparam.max_iter       = 500;
    panocparam.print_interval = 1;
    // Settings for the L-BFGS algorithm used by PANOC
    Direction::AcceleratorParams lbfgsparam;
    lbfgsparam.memory = 2;

    // Create an ALM solver using PANOC as inner solver
    OuterSolver solver{
        almparam,                 // params for outer solver
        {panocparam, lbfgsparam}, // inner solver
    };

    // Initial guess
    alpaqa::DefaultConfig::vec x(nx_);
    alpaqa::DefaultConfig::vec y(ng_);
    casadi_copy(d_nlp->x0, nx_, x.data());
    casadi_copy(d_nlp->lam_g0, ng_, y.data());

    // Solve the problem
    auto stats = solver(counted_problem, x, y);
    // y and x have been overwritten by the solution

    casadi_copy(x.data(), nx_, d_nlp->x);
    casadi_copy(y.data(), ng_, d_nlp->lam_g);
    *d_nlp->f = problem.eval_f(x);
    alpaqa::DefaultConfig::vec g(ng_);
    problem.eval_g(x, g);
    casadi_copy(g.data(), ng_, d_nlp->g);

    // Print the results
    std::cout << '\n' << *counted_problem.evaluations << '\n';
    std::cout << "status: " << stats.status << '\n'
              << "f = " << problem.eval_f(x) << '\n'
              << "inner iterations: " << stats.inner.iterations << '\n'
              << "outer iterations: " << stats.outer_iterations << '\n'
              << "ε = " << stats.ε << '\n'
              << "δ = " << stats.δ << '\n'
              << "elapsed time:     "
              << std::chrono::duration<double>{stats.elapsed_time}.count()
              << " s" << '\n'
              << "x = " << x.transpose() << '\n'
              << "y = " << y.transpose() << '\n'
              << "avg τ = " << (stats.inner.sum_τ / stats.inner.count_τ) << '\n'
              << "L-BFGS rejected = " << stats.inner.lbfgs_rejected << '\n'
              << "L-BFGS failures = " << stats.inner.lbfgs_failures << '\n'
              << "Line search failures = " << stats.inner.linesearch_failures
              << '\n'
              << std::endl;

    return 0;
  }


  AlpaqaInterface::AlpaqaInterface(DeserializingStream& s) : Nlpsol(s) {
    s.version("AlpaqaInterface", 1);
  }

  void AlpaqaInterface::serialize_body(SerializingStream &s) const {
    Nlpsol::serialize_body(s);
    s.version("AlpaqaInterface", 1);
  }

} // namespace casadi
