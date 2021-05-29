#include <panoc-alm/alm.hpp>
#include <panoc-alm/inner/directions/lbfgs.hpp>
#include <panoc-alm/inner/panoc.hpp>
#include <panoc-alm/inner/second-order-panoc-lbfgs.hpp>
#include <panoc-alm/inner/second-order-panoc.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>

#include <yaml-cpp/emitter.h>

using namespace pa;
using namespace std::string_literals;
using namespace std::chrono_literals;

Problem build_test_problem2() {
    pa::Problem p;
    p.n            = 2;
    p.m            = 2;
    p.C.upperbound = vec::Constant(2, 10);
    p.C.lowerbound = vec::Constant(2, -10);
    p.D.upperbound.resize(2);
    p.D.upperbound << -1, 0;
    p.D.lowerbound.resize(2);
    p.D.lowerbound << -inf, -inf;
    p.f = [](crvec x) { // f(x) = 1/48 x₁⁴ - 2x₁x₂ + 1/24 x₁²x₂⁴ + 10
        return 1. / 48 * std::pow(x(0), 4) - 2 * x(0) * x(1) +
               1. / 24 * std::pow(x(0), 2) * std::pow(x(1), 4) + 10;
    };
    p.grad_f = [](crvec x, rvec grad) {
        grad(0) = 1. / 12 * std::pow(x(0), 3) - 2 * x(1) +
                  1. / 12 * x(0) * std::pow(x(1), 4);
        grad(1) = -2 * x(0) + 1. / 6 * std::pow(x(0), 2) * std::pow(x(1), 3);
    };
    p.g = [](crvec x, rvec g) {
        g(0) = -4 * std::pow(x(0), 2) +
               0.25 * std::pow(x(0), 2) * std::pow(x(1), 2);
        g(1) = 0.125 * std::pow(x(0), 4) - x(0) * x(1);
    };
    p.grad_g_prod = [](crvec x, crvec y, rvec grad) {
        pa::mat gradmat(2, 2);
        gradmat <<                                      //
            -8 * x(0) + 0.5 * x(0) * std::pow(x(1), 2), //
            0.5 * std::pow(x(0), 3) - x(1),             //
            0.5 * std::pow(x(0), 2) * x(1),             //
            -x(0);
        grad = gradmat * y;
    };
    p.grad_gi = [](crvec x, unsigned i, rvec grad) {
        pa::mat gradmat(2, 2);
        gradmat <<                                      //
            -8 * x(0) + 0.5 * x(0) * std::pow(x(1), 2), //
            0.5 * std::pow(x(0), 3) - x(1),             //
            0.5 * std::pow(x(0), 2) * x(1),             //
            -x(0);
        grad = gradmat.col(i);
    };
    p.hess_L = [](crvec x, crvec y, rmat H) {
        // Hessian of f
        H(0, 0) = 1. / 4 * std::pow(x(0), 2) + 1. / 12 * std::pow(x(1), 4);
        H(0, 1) = -2 + 1. / 3 * x(0) * std::pow(x(1), 3);
        H(1, 0) = H(0, 1);
        H(1, 1) = 1. / 2 * std::pow(x(0), 2) * std::pow(x(1), 2);

        mat Hg1(2, 2);
        Hg1(0, 0) = 0.5 * std::pow(x(1), 2) - 8;
        Hg1(0, 1) = x(0) * x(1);
        Hg1(1, 0) = x(0) * x(1);
        Hg1(1, 1) = 0.5 * std::pow(x(0), 2);

        mat Hg2(2, 2);
        Hg2(0, 0) = 1.5 * std::pow(x(0), 2);
        Hg2(0, 1) = -1;
        Hg2(1, 0) = -1;
        Hg2(1, 1) = 0;

        H += y(0) * Hg1;
        H += y(1) * Hg2;
    };
    return p;
}
Problem build_test_problem3() {
    pa::Problem p;
    p.n            = 2;
    p.m            = 0;
    p.C.upperbound = vec::Constant(2, inf);
    p.C.lowerbound = vec::Constant(2, -inf);
    mat A(2, 2);
    A << 4, -2, -2, 3;
    p.f           = [A](crvec x) { return 0.5 * x.transpose() * A * x; };
    p.grad_f      = [A](crvec x, rvec g) { g = A * x; };
    p.g           = [](crvec, rvec) {};
    p.grad_g_prod = [](crvec, crvec, rvec grad) { grad.setZero(); };
    p.hess_L      = [A](crvec, crvec, rmat H) { H = A; };
    return p;
}

inline YAML::Emitter &operator<<(YAML::Emitter &out, crvec v) {
    out << YAML::Flow;
    out << YAML::BeginSeq;
    for (pa::vec::Index i = 0; i < v.size(); ++i)
        out << v[i];
    out << YAML::EndSeq;
    return out;
}

template <class InnersolverT>
auto logger(YAML::Emitter &ymlout) {
    return [&ymlout, alm_it{0u}](
               const typename InnersolverT::ProgressInfo &info) mutable {
        if (info.k == 0) {
            if (alm_it == 0) {
                ymlout << YAML::BeginMap;
            } else {
                ymlout << YAML::EndSeq << YAML::EndMap;
            }
            ymlout << YAML::Key << alm_it << YAML::Value << YAML::BeginMap
                   << YAML::Key << "Σ" << YAML::Value << info.Σ << YAML::Key
                   << "y" << YAML::Value << info.y << YAML::Key << "PANOC"
                   << YAML::Value << YAML::BeginSeq;
            ++alm_it;
        }
        ymlout << YAML::BeginMap;
        ymlout << YAML::Key << "x" << YAML::Value << info.x;
        ymlout << YAML::Key << "p" << YAML::Value << info.p;
        ymlout << YAML::Key << "‖p‖" << YAML::Value
               << std::sqrt(info.norm_sq_p);
        ymlout << YAML::Key << "x_hat" << YAML::Value << info.x_hat;
        ymlout << YAML::Key << "ψ" << YAML::Value << info.ψ;
        ymlout << YAML::Key << "∇ψ" << YAML::Value << info.grad_ψ;
        ymlout << YAML::Key << "ψ_hat" << YAML::Value << info.ψ_hat;
        ymlout << YAML::Key << "∇ψ_hat" << YAML::Value << info.grad_ψ_hat;
        ymlout << YAML::Key << "γ" << YAML::Value << info.γ;
        ymlout << YAML::Key << "ε" << YAML::Value << info.ε;
        ymlout << YAML::EndMap;
    };
}

int main(int argc, char *argv[]) {
    using namespace std::string_literals;
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <output folder>" << std::endl;
        return 1;
    }

    auto p = build_test_problem2();

    vec x1(p.n), y1(p.m);
    x1(0) = 0.05;
    x1(1) = 4.5;
    y1.fill(0);
    vec x2 = x1, y2 = y1;
    vec x3 = x1, y3 = y1;

    std::ofstream outf1(argv[1] + "/1.yaml"s);
    YAML::Emitter ymlout1(outf1);
    std::ofstream outf2(argv[1] + "/2.yaml"s);
    YAML::Emitter ymlout2(outf2);
    std::ofstream outf3(argv[1] + "/3.yaml"s);
    YAML::Emitter ymlout3(outf3);

    pa::ALMParams almparams;
    almparams.max_iter        = 20;
    almparams.max_time        = 1min + 30s;
    almparams.preconditioning = false;
    almparams.print_interval  = 1;
    almparams.Σ₀              = 1e-1;
    almparams.Δ               = 10;
    pa::PANOCParams panocparams;
    panocparams.max_iter                       = 1000;
    panocparams.update_lipschitz_in_linesearch = true;
    panocparams.alternative_linesearch_cond    = false;
    panocparams.print_interval                 = 1;
    pa::LBFGSParams lbfgsparams;
    lbfgsparams.rescale_when_γ_changes = false;
    lbfgsparams.memory                 = 5;

    pa::SecondOrderPANOCParams panocparams2;
    panocparams2.max_iter                       = 1000;
    panocparams2.update_lipschitz_in_linesearch = true;
    panocparams2.alternative_linesearch_cond    = false;
    panocparams2.print_interval                 = 1;

    pa::SecondOrderPANOCLBFGSParams panocparams3;
    panocparams3.max_iter                       = 1000;
    panocparams3.update_lipschitz_in_linesearch = true;
    panocparams3.alternative_linesearch_cond    = false;
    panocparams3.print_interval                 = 1;
    pa::LBFGSParams lbfgsparams3;
    lbfgsparams3.memory = 5;

    ALMSolver<PANOCSolver<>> solver1{almparams, {panocparams, lbfgsparams}};
    ALMSolver<SecondOrderPANOCSolver> solver2{almparams, panocparams2};
    ALMSolver<SecondOrderPANOCLBFGSSolver> solver3{
        almparams, {panocparams3, lbfgsparams3}};
    solver1.inner_solver.set_progress_callback(
        logger<decltype(solver1.inner_solver)>(ymlout1));
    solver2.inner_solver.set_progress_callback(
        logger<decltype(solver2.inner_solver)>(ymlout2));
    solver3.inner_solver.set_progress_callback(
        logger<decltype(solver3.inner_solver)>(ymlout3));

    solver1(p, y1, x1);
    std::cout << std::endl << "=======" << std::endl << std::endl;
    solver2(p, y2, x2);
    std::cout << std::endl << "=======" << std::endl << std::endl;
    solver3(p, y3, x3);
}