#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/print.hpp>

#include <iomanip>
#include <iostream>

struct Problem {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    using Box = alpaqa::Box<config_t>;

    length_t N = 32,      // Horizon length
        nu     = 1,       // Number of inputs
        nx     = 2,       // Number of states
        nh     = nu + nx, // Number of stage outputs
        nh_N   = nx,      // Number of terminal outputs
        nc     = 0,       // Number of stage constraints
        nc_N   = 0;       // Number of terminal constraints

    mat A, B;

    Problem() : A(nx, nx), B(nx, nu) {
        A.setIdentity();
        B.setIdentity();
    }

    // Horizon length
    [[nodiscard]] length_t get_N() const { return N; }
    // Number of inputs
    [[nodiscard]] length_t get_nu() const { return nu; }
    // Number of states
    [[nodiscard]] length_t get_nx() const { return nx; }
    // Number of stage outputs
    [[nodiscard]] length_t get_nh() const { return nh; }
    // Number of terminal outputs
    [[nodiscard]] length_t get_nh_N() const { return nh_N; }
    // Number of stage constraints
    [[nodiscard]] length_t get_nc() const { return nc; }
    // Number of terminal constraints
    [[nodiscard]] length_t get_nc_N() const { return nc_N; }

    // Bound constraints on the inputs u.
    void get_U(Box &U) const {
        U.lowerbound.setConstant(-1);
        U.upperbound.setConstant(+1);
    }
    // Bound constraints on function of the states
    void get_D([[maybe_unused]] Box &D) const {}
    // Bound constraints on function of the terminal state
    void get_D_N([[maybe_unused]] Box &D) const {}

    // Get the initial state
    void get_x_init(rvec x_init) const { x_init.setConstant(10.); }

    // Discrete-time dynamics x⁺ = f(x, u)
    void eval_f([[maybe_unused]] index_t timestep, crvec x, crvec u,
                rvec fxu) const {
        fxu.noalias() = A * x + B * u;
    }
    // Jacobian of dynamics
    void eval_jac_f([[maybe_unused]] index_t timestep, [[maybe_unused]] crvec x,
                    [[maybe_unused]] crvec u, rmat J_fxu) const {
        J_fxu.leftCols(nx).noalias()  = A;
        J_fxu.rightCols(nu).noalias() = B;
    }
    // Gradient-vector product of dynamics
    void eval_grad_f_prod([[maybe_unused]] index_t timestep,
                          [[maybe_unused]] crvec x, [[maybe_unused]] crvec u,
                          crvec p, rvec grad_fxu_p) const {
        grad_fxu_p.topRows(nx).noalias()    = A.transpose() * p;
        grad_fxu_p.bottomRows(nu).noalias() = B.transpose() * p;
    }
    // Output mapping
    void eval_h([[maybe_unused]] index_t timestep, crvec x, crvec u,
                rvec h) const {
        h.topRows(nx)    = x;
        h.bottomRows(nu) = u;
    }
    // Terminal output mapping
    void eval_h_N(crvec x, rvec h) const { h = x; }
    // Stage cost
    [[nodiscard]] real_t eval_l([[maybe_unused]] index_t timestep,
                                crvec h) const {
        return 0.5 * h.squaredNorm();
    }
    // Terminal cost
    [[nodiscard]] real_t eval_l_N(crvec h) const {
        return 5.0 * h.squaredNorm();
    }
    // Gradient of the cost l(h(x, u))
    void eval_qr([[maybe_unused]] index_t timestep, [[maybe_unused]] crvec xu,
                 crvec h, rvec qr) const {
        auto Jh_xu    = mat::Identity(nx + nu, nx + nu);
        auto &&grad_l = h;
        qr            = Jh_xu.transpose() * grad_l;
    }
    // Gradient of the terminal cost l_N(h_N(x))
    void eval_q_N([[maybe_unused]] crvec x, crvec h, rvec q) const {
        auto Jh_x     = mat::Identity(nx, nx);
        auto &&grad_l = 10 * h;
        q             = Jh_x.transpose() * grad_l;
    }
    // Hessian of stage cost w.r.t. x
    void eval_add_Q([[maybe_unused]] index_t timestep,
                    [[maybe_unused]] crvec xu, [[maybe_unused]] crvec h,
                    rmat Q) const {
        Q += mat::Identity(nx, nx);
    }
    // Hessian of terminal cost
    void eval_add_Q_N([[maybe_unused]] crvec x, [[maybe_unused]] crvec h,
                      rmat Q) const {
        Q += 10 * mat::Identity(nx, nx);
    }
    // Hessian of stage cost w.r.t. u
    void eval_add_R_masked([[maybe_unused]] index_t timestep,
                           [[maybe_unused]] crvec xu, [[maybe_unused]] crvec h,
                           crindexvec mask, rmat R,
                           [[maybe_unused]] rvec work) const {
        auto R_full = mat::Identity(nu, nu);
        R.noalias() += R_full(mask, mask);
    }
    // Hessian of stage cost w.r.t. x and u
    void eval_add_S_masked([[maybe_unused]] index_t timestep,
                           [[maybe_unused]] crvec xu, [[maybe_unused]] crvec h,
                           crindexvec mask, rmat S,
                           [[maybe_unused]] rvec work) const {
        // Mixed derivatives are zero, so the following has no effect
        using Eigen::indexing::all;
        auto S_full = mat::Zero(nu, nx);
        S += S_full(mask, all);
    }
    // Hessian-vector product of stage cost w.r.t. u
    void eval_add_R_prod_masked([[maybe_unused]] index_t timestep,
                                [[maybe_unused]] crvec xu,
                                [[maybe_unused]] crvec h, crindexvec mask_J,
                                crindexvec mask_K, crvec v, rvec out,
                                [[maybe_unused]] rvec work) const {
        // R is diagonal, and J ∩ K = ∅, so the following has no effect
        auto R_full = mat::Identity(nu, nu);
        out.noalias() += R_full(mask_J, mask_K) * v(mask_K);
    }
    // Hessian-vector product of stage cost w.r.t. u and x
    void eval_add_S_prod_masked([[maybe_unused]] index_t timestep,
                                [[maybe_unused]] crvec xu,
                                [[maybe_unused]] crvec h,
                                [[maybe_unused]] crindexvec mask_K,
                                [[maybe_unused]] crvec v,
                                [[maybe_unused]] rvec out,
                                [[maybe_unused]] rvec work) const {
        // Mixed derivatives are zero, so the following has no effect
        using Eigen::indexing::all;
        auto Sᵀ = mat::Zero(nx, nu);
        out.noalias() += Sᵀ(all, mask_K) * v(mask_K);
    }
    [[nodiscard]] length_t get_R_work_size() const {
        // No workspace needed to share data between eval_add_R_masked and
        // eval_add_R_prod_masked
        return 0;
    }
    [[nodiscard]] length_t get_S_work_size() const {
        // No workspace needed to share data between eval_add_S_masked and
        // eval_add_S_prod_masked
        return 0;
    }

    // Project the Lagrange multipliers for ALM
    void eval_proj_multipliers([[maybe_unused]] rvec y,
                               [[maybe_unused]] real_t M) const {
        // Empty for unconstrained problems
    }

    // Get the distance to the feasible set
    void eval_proj_diff_g([[maybe_unused]] crvec z,
                          [[maybe_unused]] rvec e) const {
        // Empty for unconstrained problems
    }

    void check() const {
        // You could do some sanity checks here
    }
};

int main() {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    // Problem
    alpaqa::ControlProblemWithCounters<Problem> problem;

    // Problem dimensions
    const auto n = problem.get_N() * problem.get_nu(),
               m = problem.get_N() * problem.get_nc() + problem.get_nc_N();

    // Initial guess and other solver inputs
    vec u = vec::Zero(n); // Inputs (single shooting)
    vec y = vec::Zero(m); // Lagrange multipliers
    vec μ = vec::Ones(m); // Penalty factors
    vec e(m);             // Constraint violation

    // Solver
    alpaqa::PANOCOCPParams<config_t> params;
    params.stop_crit      = alpaqa::PANOCStopCrit::ProjGradUnitNorm;
    params.gn_interval    = 1;
    params.print_interval = 1;
    alpaqa::PANOCOCPSolver<config_t> solver{params};

    // Solve
    auto stats = solver(problem, {.tolerance = 1e-8}, u, y, μ, e);

    // Print statistics
    auto δ      = e.lpNorm<Eigen::Infinity>();
    auto time_s = std::chrono::duration<double>(stats.elapsed_time).count();
    std::cout << '\n'
              << *problem.evaluations << '\n'
              << "solver:  " << solver.get_name() << '\n'
              << "status:  " << stats.status << '\n'
              << "ψ = " << alpaqa::float_to_str(stats.final_ψ) << '\n'
              << "ε = " << alpaqa::float_to_str(stats.ε) << '\n'
              << "δ = " << alpaqa::float_to_str(δ) << '\n'
              << "time: " << alpaqa::float_to_str(time_s, 3) << " s\n"
              << "iter:      " << std::setw(6) << stats.iterations << '\n'
              << "line search backtrack: " << std::setw(6)
              << stats.linesearch_backtracks << '\n'
              << "step size backtrack:   " << std::setw(6)
              << stats.stepsize_backtracks << '\n'
              << "solution: ";
    alpaqa::print_python(std::cout, u) << std::endl;
}