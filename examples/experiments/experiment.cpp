#include <panoc-alm/decl/alm.hpp>
#include <panoc-alm/inner/decl/panoc.hpp>
#include <panoc-alm/inner/directions/decl/lbfgs.hpp>

#include <panoc-alm/interop/cutest/CUTEstLoader.hpp>

#include <atomic>
#include <chrono>
#include <csignal>
#include <fstream>
#include <iostream>
#include <sstream>
#include <yaml-cpp/emitter.h>
#include <yaml-cpp/emittermanip.h>

using namespace std::chrono_literals;

using Solver = pa::ALMSolver<>;

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
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <problem name> <output folder>"
                  << std::endl;
        return 1;
    }
    std::string prob_dir = "CUTEst/"s + argv[1];
    CUTEstProblem cp(prob_dir + "/libcutest-problem-" + argv[1] + ".so",
                     prob_dir + "/OUTSDIF.d");

    unsigned alm_it = 0;
    std::ofstream outf(argv[2] + "/"s + argv[1] + ".yaml");
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
        ymlout << YAML::Key << "ε" << YAML::Value << info.ε;
        ymlout << YAML::EndMap;
    };

    pa::ALMParams almparams;
    almparams.max_iter        = 20;
    almparams.max_time        = 1min + 30s;
    almparams.preconditioning = false;
    // almparams.print_interval  = 1;
    almparams.Σ₀ = 1;
    // almparams.ε₀ = 1e-5;
    // almparams.Δ = 1.1;
    pa::PANOCParams panocparams;
    panocparams.max_iter                       = 1000;
    panocparams.update_lipschitz_in_linesearch = true;
    panocparams.alternative_linesearch_cond    = false;
    panocparams.lbfgs_mem                      = 20;
    // panocparams.print_interval = 500;
    pa::LBFGSParams lbfgsparams;

    Solver solver{almparams, {panocparams, lbfgsparams}};
    solver.inner_solver.set_progress_callback(logger);
    std::atomic_signal_fence(std::memory_order_release);
    acitve_solver.store(&solver, std::memory_order_relaxed);
    signal(SIGINT, signal_callback_handler);

    std::cout << solver.get_name() << std::endl;

    auto problem = pa::ProblemWithCounters(cp.problem);

    auto status = solver(problem, cp.y0, cp.x0);
    // ??? TODO: fence
    acitve_solver.store(nullptr, std::memory_order_relaxed);
    // ??? TODO: fence
    // auto report = cp.get_report();

    ymlout << YAML::EndMap << YAML::EndSeq << YAML::EndMap;
    std::cout << status.status << ": (" << status.outer_iterations << ", "
              << status.inner_iterations << ")" << std::endl;
}