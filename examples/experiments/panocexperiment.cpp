#include <panoc-alm/inner/decl/panoc.hpp>
#include <panoc-alm/inner/directions/decl/lbfgs.hpp>
#include <panoc-alm/reference-problems/riskaverse-mpc.hpp>

#include <atomic>
#include <chrono>
#include <csignal>
#include <fstream>
#include <iostream>
#include <sstream>
#include <yaml-cpp/emitter.h>
#include <yaml-cpp/emittermanip.h>

using namespace std::chrono_literals;

using Solver = pa::PANOCSolver<>;

std::atomic<Solver *> acitve_solver{nullptr};
void signal_callback_handler(int signum) {
    if (signum == SIGINT) {
        if (auto *s = acitve_solver.load(std::memory_order_relaxed)) {
            std::atomic_signal_fence(std::memory_order_acquire);
            s->stop();
        }
    }
}

inline YAML::Emitter &operator<<(YAML::Emitter &out, const pa::vec &v) {
    out << YAML::Flow;
    out << YAML::BeginSeq;
    for (pa::vec::Index i = 0; i < v.size(); ++i)
        out << v[i];
    out << YAML::EndSeq;
    return out;
}

int main(int argc, char *argv[]) {
    using namespace std::string_literals;
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <output folder>" << std::endl;
        return 1;
    }
    unsigned alm_it = 0;
    std::ofstream outf(argv[1] + "/"s + "panoc-experiment" + ".yaml");
    YAML::Emitter ymlout(outf);
    auto logger = [&ymlout,
                   &alm_it](const pa::PANOCSolver<>::ProgressInfo &info) {
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
        ymlout << YAML::Key << "L" << YAML::Value << info.L;
        ymlout << YAML::Key << "ε" << YAML::Value << info.ε;
        ymlout << YAML::EndMap;
    };

    pa::PANOCParams panocparams;
    panocparams.max_iter                       = 1000;
    panocparams.update_lipschitz_in_linesearch = true;
    panocparams.alternative_linesearch_cond    = false;
    panocparams.lbfgs_mem                      = 20;
    // panocparams.print_interval = 500;
    pa::LBFGSParams lbfgsparams;

    Solver solver{panocparams, lbfgsparams};
    solver.set_progress_callback(logger);
    std::atomic_signal_fence(std::memory_order_release);
    acitve_solver.store(&solver, std::memory_order_relaxed);
    signal(SIGINT, signal_callback_handler);

    std::cout << solver.get_name() << std::endl;

    auto riskaverseproblem = pa::problems::riskaverse_mpc_problem();

    auto problem = pa::ProblemWithCounters(riskaverseproblem);

    pa::vec Σ(problem.m);
    Σ << 3.0797607617330747, 3.9569193214493699, 4.9537113927653671;
    pa::real_t ε = 1.0000000000000004e-11;
    pa::vec x(problem.n);
    x << -4.8780397699768994, -4.8780397699768994, 8156.220204170053,
        4108.6295384399791, 8156.792300872844, 0.57213825479611202, 0;
    pa::vec y(problem.m);
    y << 1.0000034126079587, 1.0000097609953398, 1.0000109987521664;
    pa::vec err_z(problem.m);
    auto status = solver(problem, Σ, ε, true, x, y, err_z);
    // ??? TODO: fence
    acitve_solver.store(nullptr, std::memory_order_relaxed);
    // ??? TODO: fence

    ymlout << YAML::EndMap << YAML::EndSeq << YAML::EndMap;
    std::cout << status.status << ": (" << status.iterations << ")"
              << std::endl;
}