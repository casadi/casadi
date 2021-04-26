#include "panoc-alm/inner/directions/decl/lbfgs-fwd.hpp"
#include "panoc-alm/util/vec.hpp"
#include <limits>
#include <panoc-alm/inner/directions/lbfgs.hpp>
#include <panoc-alm/inner/panoc.hpp>
#include <panoc-alm/inner/pga.hpp>
#include <panoc-alm/inner/second-order-panoc-lbfgs.hpp>
#include <panoc-alm/inner/second-order-panoc.hpp>
#include <panoc-alm/reference-problems/himmelblau.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>

#include <yaml-cpp/emitter.h>

using namespace pa;
using namespace std::string_literals;
using namespace std::chrono_literals;

Problem build_test_problem() {
    pa::Problem p     = pa::problems::himmelblau_problem();
    p.C.lowerbound(0) = -inf;
    p.C.lowerbound(1) = 2.1;
    p.C.upperbound(0) = +inf;
    p.C.upperbound(1) = +inf;
    return p;
}

template <>
inline void
YAML::Emitter::SetStreamablePrecision<long double>(std::stringstream &stream) {
    stream.precision(static_cast<std::streamsize>(
        std::numeric_limits<long double>::max_digits10));
}

inline YAML::Emitter &operator<<(YAML::Emitter &emitter, long double v) {
    return emitter.WriteStreamable(v);
}

inline YAML::Emitter &operator<<(YAML::Emitter &out, const pa::vec &v) {
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

    auto p = build_test_problem();

    vec x0(p.n), y0(p.m);
    x0(0) = 3.3;
    x0(1) = 3;
    y0.fill(0);
    vec x1 = x0, y1 = y0;
    vec x2 = x0, y2 = y0;
    vec x3 = x0, y3 = y0;

    std::string pf = argc == 2 ? (argv[1] + "/"s) : (argv[1] + "/long-"s);
    std::string sf = ".yaml";

    std::ofstream outf0(pf + "pga" + sf);
    YAML::Emitter ymlout0(outf0);
    std::ofstream outf1(pf + "panoc" + sf);
    YAML::Emitter ymlout1(outf1);
    std::ofstream outf2(pf + "panoc-2lbfgs" + sf);
    YAML::Emitter ymlout2(outf2);
    std::ofstream outf3(pf + "panoc-2newton" + sf);
    YAML::Emitter ymlout3(outf3);

    auto L₀        = 0.95 / 5e-3;
    auto lbfgs_mem = 5;

    pa::PGAParams pgaparams0;
    pgaparams0.Lipschitz.L₀   = L₀;
    pgaparams0.max_iter       = 1000;
    pgaparams0.print_interval = 1;

    pa::PANOCParams panocparams1;
    panocparams1.Lipschitz.L₀                   = L₀;
    panocparams1.max_iter                       = 1000;
    panocparams1.update_lipschitz_in_linesearch = true;
    panocparams1.alternative_linesearch_cond    = false;
    panocparams1.print_interval                 = 1;
    panocparams1.lbfgs_stepsize = LBFGSStepSize::BasedOnGradientStepSize;
    pa::LBFGSParams lbfgsparams;
    lbfgsparams.rescale_when_γ_changes = false;
    lbfgsparams.memory                 = lbfgs_mem;

    pa::SecondOrderPANOCParams panocparams3;
    panocparams3.Lipschitz.L₀                   = L₀;
    panocparams3.max_iter                       = 1000;
    panocparams3.update_lipschitz_in_linesearch = true;
    panocparams3.alternative_linesearch_cond    = false;
    panocparams3.print_interval                 = 1;

    pa::SecondOrderPANOCLBFGSParams panocparams2;
    panocparams2.Lipschitz.L₀                   = L₀;
    panocparams2.max_iter                       = 1000;
    panocparams2.update_lipschitz_in_linesearch = true;
    panocparams2.alternative_linesearch_cond    = false;
    panocparams2.lbfgs_mem                      = lbfgs_mem;
    panocparams2.print_interval                 = 1;
    panocparams2.lbfgs_stepsize = LBFGSStepSize::BasedOnGradientStepSize;
    pa::LBFGSParams lbfgsparams2;

    PGASolver solver0{pgaparams0};
    PANOCSolver<> solver1{panocparams1, lbfgsparams};
    SecondOrderPANOCLBFGSSolver solver2{panocparams2, lbfgsparams2};
    SecondOrderPANOCSolver solver3{panocparams3};
    solver0.set_progress_callback(logger<decltype(solver0)>(ymlout0));
    solver1.set_progress_callback(logger<decltype(solver1)>(ymlout1));
    solver2.set_progress_callback(logger<decltype(solver2)>(ymlout2));
    solver3.set_progress_callback(logger<decltype(solver3)>(ymlout3));

    const pa::vec Σ = pa::vec::Zero(p.m);
    pa::vec err(p.m);
    const pa::real_t ε =
        argc == 2 ? 1e-9 : std::numeric_limits<real_t>::epsilon();

    unsigned iters[4];

    iters[0] = solver0(p, Σ, ε, true, x0, y0, err).iterations;
    std::cout << std::endl << "=======" << std::endl << std::endl;
    iters[1] = solver1(p, Σ, ε, true, x1, y1, err).iterations;
    std::cout << std::endl << "=======" << std::endl << std::endl;
    iters[2] = solver2(p, Σ, ε, true, x2, y2, err).iterations;
    std::cout << std::endl << "=======" << std::endl << std::endl;
    iters[3] = solver3(p, Σ, ε, true, x3, y3, err).iterations;

    for (auto it : iters)
        std::cout << it << ", ";
    std::cout << std::endl;
    auto digits = std::numeric_limits<real_t>::max_digits10;
    std::cout << std::setprecision(digits) << x0.transpose() << std::endl;
}