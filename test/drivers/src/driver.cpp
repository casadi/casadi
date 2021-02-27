#include "panoc-alm/alm.hpp"
#include "panoc-alm/panoc.hpp"
#include <drivers/YAMLEncoder.hpp>
#include <fstream>
#include <iostream>
#include <panoc-alm/interop/cutest/CUTEstLoader.hpp>
#include <yaml-cpp/emittermanip.h>

int main(int argc, char *argv[]) {
    using namespace std::string_literals;
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <problem name>" << std::endl;
        return 1;
    }
    std::string prob_dir = "CUTEst/"s + argv[1];
    CUTEstProblem cp(prob_dir + "/libcutest-problem-" + argv[1] + ".so",
                     prob_dir + "/OUTSDIF.d");

    pa::ALMParams almparams;
    almparams.max_iter        = 1000;
    almparams.preconditioning = false;
    pa::PANOCParams panocparams;
    panocparams.max_iter = 1000;
    panocparams.experimental.specialized_lbfgs = false;
    panocparams.experimental.update_lipschitz_in_linesearch = false;

    pa::ALMSolver solver(almparams, panocparams);

    auto status = solver(cp.problem, cp.y0, cp.x0);
    auto report = cp.get_report();

    auto f_star = cp.problem.f(cp.x0);

    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << "x" << YAML::Value << cp.x0;
    out << YAML::Key << "y" << YAML::Value << cp.y0;
    out << YAML::Key << "outer iterations" << YAML::Value
        << status.outer_iterations;
    out << YAML::Key << "inner iterations" << YAML::Value
        << status.inner_iterations;
    out << YAML::Key << "inner convergence failures" << YAML::Value
        << status.inner_convergence_failures;
    out << YAML::Key << "ε" << YAML::Value << status.ε;
    out << YAML::Key << "δ" << YAML::Value << status.δ;
    out << YAML::Key << "f" << YAML::Value << f_star;
    out << YAML::Key << "failed" << YAML::Value << status.failed;
    out << YAML::Key << "converged" << YAML::Value << status.converged;
    out << YAML::EndMap;
    out << report;

    std::ofstream("testresults/CUTEst/"s + argv[1] + ".yaml")
        << out.c_str() << std::endl;

    std::cout << out.c_str() << std::endl;
}