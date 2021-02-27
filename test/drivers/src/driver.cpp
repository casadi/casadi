#include "panoc-alm/alm.hpp"
#include "panoc-alm/panoc.hpp"
#include <chrono>
#include <drivers/YAMLEncoder.hpp>
#include <fstream>
#include <iostream>
#include <panoc-alm/interop/cutest/CUTEstLoader.hpp>
#include <sstream>
#include <yaml-cpp/emittermanip.h>

using namespace std::chrono_literals;

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

    pa::ALMParams almparams;
    almparams.max_iter        = 1000;
    almparams.max_time        = 2min + 30s;
    almparams.preconditioning = false;
    pa::PANOCParams panocparams;
    panocparams.max_iter                                    = 1000;
    panocparams.experimental.specialized_lbfgs              = false;
    panocparams.experimental.update_lipschitz_in_linesearch = false;

    pa::ALMSolver solver(almparams, panocparams);

    auto problem = pa::ProblemWithCounters(cp.problem);

    auto status = solver(problem, cp.y0, cp.x0);
    auto report = cp.get_report();

    auto f_star = cp.problem.f(cp.x0);

    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << "status" << YAML::Value << status.status;
    out << YAML::Key << "outer iterations" << YAML::Value
        << status.outer_iterations;
    out << YAML::Key << "inner iterations" << YAML::Value
        << status.inner_iterations;
    out << YAML::Key << "inner convergence failures" << YAML::Value
        << status.inner_convergence_failures;
    out << YAML::Key << "elapsed time" << YAML::Value
        << std::chrono::duration<double>(status.elapsed_time).count();
    out << YAML::Key << "ε" << YAML::Value << status.ε;
    out << YAML::Key << "δ" << YAML::Value << status.δ;
    out << YAML::Key << "f" << YAML::Value << f_star;
    out << YAML::Key << "counters" << YAML::Value << problem.evaluations;
    out << YAML::Key << "linesearch failures" << YAML::Value
        << status.inner_linesearch_failures;
    out << YAML::Key << "L-BFGS failures" << YAML::Value
        << status.inner_lbfgs_failures;
    out << YAML::Key << "L-BFGS rejected" << YAML::Value
        << status.inner_lbfgs_rejected;
    out << YAML::Key << "x" << YAML::Value << cp.x0;
    out << YAML::Key << "y" << YAML::Value << cp.y0;
    out << YAML::EndMap;
    out << report;

    std::ofstream(argv[2] + "/"s + argv[1] + ".yaml")
        << out.c_str() << std::endl;

    std::cout << out.c_str() << std::endl;
}