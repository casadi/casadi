#include "panoc-alm/inner/decl/second-order-panoc-lbfgs.hpp"
#include "panoc-alm/util/problem.hpp"
#include <panoc-alm/alm.hpp>
#include <panoc-alm/inner/decl/panoc.hpp>
#include <panoc-alm/inner/directions/decl/lbfgs.hpp>
#include <panoc-alm/inner/second-order-panoc-lbfgs.hpp>

#include <panoc-alm/interop/casadi/CasADiLoader.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace fs = std::filesystem;

template <typename M>
M load_csv(const std::string &path) {
    std::ifstream indata(path);
    assert(indata);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    return Eigen::Map<const Eigen::Matrix<
        typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime,
        M::ColsAtCompileTime == 1 ? 0 : Eigen::RowMajor>>(values.data(), rows,
                                                          values.size() / rows);
}

using pa::crvec;
using pa::real_t;
using pa::rvec;
using pa::vec;

int main(int argc, char *argv[]) {
    fs::path path = "examples/mpc/mmehrez";
    // if (argc > 1)
    //     path = argv[1];
    fs::path so_name = path / "libmmehrez_functions.so";

    pa::Problem p;
    p.C.lowerbound = load_csv<vec>(path / "mmehrez_functions-lbx.csv");
    p.C.upperbound = load_csv<vec>(path / "mmehrez_functions-ubx.csv");
    p.D.lowerbound = load_csv<vec>(path / "mmehrez_functions-lbg.csv");
    p.D.upperbound = load_csv<vec>(path / "mmehrez_functions-ubg.csv");
    assert(p.C.upperbound.size() == p.C.lowerbound.size());
    assert(p.D.upperbound.size() == p.D.lowerbound.size());
    p.n           = p.C.lowerbound.size();
    p.m           = p.D.lowerbound.size();
    p.f           = load_CasADi_objective(so_name.c_str());
    p.grad_f      = load_CasADi_gradient_objective(so_name.c_str());
    p.g           = load_CasADi_constraints(so_name.c_str());
    p.grad_g_prod = load_CasADi_gradient_constraints_prod(so_name.c_str());
    vec w         = vec::Zero(p.m);
    p.grad_gi     = [&p, w](crvec x, unsigned i, rvec g) mutable {
        w(i) = 1;
        p.grad_g_prod(x, w, g);
        w(i) = 0;
    };
    p.hess_L      = load_CasADi_hessian_lagrangian(so_name.c_str());
    p.hess_L_prod = load_CasADi_hessian_lagrangian_prod(so_name.c_str());

    pa::ALMParams almparam;
    almparam.ε               = 1e-5;
    almparam.δ               = 1e-5;
    almparam.Δ               = 10;
    almparam.max_iter        = 200;
    almparam.Σ₀              = 5e5;
    almparam.print_interval  = 1;
    almparam.preconditioning = false;
    almparam.max_num_retries = 5;

    pa::PANOCParams panocparam;
    panocparam.max_iter       = 1000;
    panocparam.print_interval = 50;

    pa::SecondOrderPANOCLBFGSParams secondpanocparam;
    secondpanocparam.max_iter       = panocparam.max_iter;
    secondpanocparam.print_interval = panocparam.print_interval;

    pa::LBFGSParams lbfgsparam;
    lbfgsparam.memory = 20;

    using Stats =
        std::variant<pa::ALMSolver<pa::PANOCSolver<pa::LBFGS>>::Stats,
                     pa::ALMSolver<pa::SecondOrderPANOCLBFGSSolver>::Stats>;
    using Solver_sig =
        std::function<Stats(const pa::Problem &p, rvec x, rvec y)>;

    std::vector<Solver_sig> solvers{
        pa::ALMSolver<pa::PANOCSolver<pa::LBFGS>>{
            almparam,
            {panocparam, lbfgsparam},
        },
        pa::ALMSolver<pa::SecondOrderPANOCLBFGSSolver>{
            almparam,
            {secondpanocparam, lbfgsparam},
        },
    };

    size_t solv_idx = 0;
    if (argc > 1)
        solv_idx = std::stoul(argv[1]);
    if (solv_idx >= solvers.size())
        throw std::out_of_range("invalid solver index");

    vec x = vec::Zero(p.n);
    vec y = vec::Zero(p.m);

    pa::ProblemWithCounters pc{p};
    Stats stats = solvers[solv_idx](pc, y, x);

    std::cout << "x = " << x.transpose() << std::endl;
    std::cout << "y = " << y.transpose() << std::endl;
    vec g(p.m);
    p.g(x, g);
    std::cout << "g = " << g.transpose() << std::endl;
    std::cout << "f = " << p.f(x) << std::endl;

    std::cout << "f:           " << pc.evaluations.f << std::endl;
    std::cout << "grad_f:      " << pc.evaluations.grad_f << std::endl;
    std::cout << "g:           " << pc.evaluations.g << std::endl;
    std::cout << "grad_g_prod: " << pc.evaluations.grad_g_prod << std::endl;
    std::cout << "grad_gi:     " << pc.evaluations.grad_gi << std::endl;
    std::cout << "hess_L_prod: " << pc.evaluations.hess_L_prod << std::endl;
    std::cout << "hess_L:      " << pc.evaluations.hess_L << std::endl;

    std::cout << "Solver: " << solv_idx << std::endl;
    std::visit(
        [](auto &&s) {
            auto color = s.status == pa::SolverStatus::Converged ? "\x1b[0;32m"
                                                                 : "\x1b[0;31m";
            auto color_end = "\x1b[0m";
            std::cout << "initial_penalty_reduced:    "
                      << s.initial_penalty_reduced << std::endl;
            std::cout << "penalty_reduced:            " << s.penalty_reduced
                      << std::endl;
            std::cout << "inner_convergence_failures: "
                      << s.inner_convergence_failures << "/"
                      << s.inner.iterations << std::endl;
            std::cout << "ε:   " << s.ε << std::endl;
            std::cout << "δ:   " << s.δ << std::endl;
            std::cout << "‖Σ‖: " << s.norm_penalty << std::endl;
            std::cout << "avg τ:        " << s.inner.sum_τ / s.inner.count_τ
                      << std::endl;
            std::cout << "τ=1 accepted: " << s.inner.τ_1_accepted << "/"
                      << s.inner.count_τ << " = "
                      << 100. * s.inner.τ_1_accepted / s.inner.count_τ << "%"
                      << std::endl;
            std::cout << "inner: " << s.inner.iterations << std::endl;
            std::cout << "outer: " << s.outer_iterations << std::endl;
            std::cout << "status: " << color << s.status << color_end
                      << std::endl;
        },
        stats);
}