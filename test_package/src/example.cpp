#include <alpaqa/panoc-alm.hpp>
#include <alpaqa/problem/unconstr-problem.hpp>

int main() {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    struct Problem : alpaqa::UnconstrProblem<config_t> {
        length_t get_n() const { return 1; }
        real_t eval_f(crvec) const { return 0; }
        void eval_grad_f(crvec, rvec gr) const { gr.setZero(); }
    };

    Problem problem;

    using Direction   = alpaqa::LBFGSDirection<config_t>;
    using InnerSolver = alpaqa::PANOCSolver<Direction>;
    using OuterSolver = alpaqa::ALMSolver<InnerSolver>;

    vec x = vec::Zero(1), y;
    OuterSolver solver{{}, InnerSolver{{}}};
    auto stats = solver(problem, x, y);
    return stats.status == alpaqa::SolverStatus::Converged ? 0 : 1;
}
