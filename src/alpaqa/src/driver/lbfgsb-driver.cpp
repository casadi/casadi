#include <alpaqa/implementation/outer/alm.tpp>
#include <alpaqa/lbfgsb/lbfgsb-adapter.hpp>

#include "alm-driver.hpp"
#include "cancel.hpp"
#include "lbfgsb-driver.hpp"
#include "solver-driver.hpp"
#include "util.hpp"

namespace {

using InnerLBFGSBSolver = alpaqa::lbfgsb::LBFGSBSolver;

auto make_inner_lbfgsb_solver(Options &opts) {
    // Settings for the solver
    InnerLBFGSBSolver::Params solver_param;
    solver_param.max_iter       = 50'000;
    solver_param.print_interval = 0;
    solver_param.stop_crit      = alpaqa::PANOCStopCrit::ProjGradUnitNorm;
    set_params(solver_param, "solver", opts);
    return InnerLBFGSBSolver{solver_param};
}

} // namespace

solver_func_t make_lbfgsb_driver(std::string_view direction, Options &opts) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    if (!direction.empty())
        throw std::invalid_argument(
            "L-BFGS-B solver does not support any directions");
    auto inner_solver = make_inner_lbfgsb_solver(opts);
    auto solver       = make_alm_solver(std::move(inner_solver), opts);
    unsigned N_exp    = 0;
    set_params(N_exp, "num_exp", opts);
    return [solver{std::move(solver)},
            N_exp](LoadedProblem &problem,
                   std::ostream &os) mutable -> SolverResults {
        static std::atomic<decltype(solver) *> solver_to_stop;
        auto cancel = attach_cancellation<solver_to_stop>(solver);
        return run_alm_solver(problem, solver, os, N_exp);
    };
}

template class alpaqa::ALMSolver<InnerLBFGSBSolver>;
