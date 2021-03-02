#include <cfenv>
#include <iomanip>
#include <iostream>
#include <panoc-alm/alm.hpp>
#include <panoc-alm/problem.hpp>
#include <panoc-alm/solverstatus.hpp>
#include <panoc-alm/vec.hpp>

namespace pa {

using std::chrono::duration_cast;
using std::chrono::microseconds;

namespace detail {

void project_y(vec &y,          // inout
               const vec &z_lb, // in
               const vec &z_ub, // in
               real_t M         // in
) {
    constexpr real_t inf = std::numeric_limits<real_t>::infinity();
    // TODO: Handle NaN correctly
    auto max_lb = [M](real_t y, real_t z_lb) {
        real_t y_lb = z_lb == -inf ? 0 : -M;
        return std::max(y, y_lb);
    };
    y = y.binaryExpr(z_lb, max_lb);

    auto min_ub = [M](real_t y, real_t z_ub) {
        real_t y_ub = z_ub == inf ? 0 : M;
        return std::min(y, y_ub);
    };
    y = y.binaryExpr(z_ub, min_ub);
}

void update_penalty_weights(const ALMParams &params, unsigned iteration, vec &e,
                            vec &old_e, real_t norm_e, vec &Σ) {
    if (norm_e <= params.δ)
        return;
#if 1
    for (Eigen::Index i = 0; i < e.rows(); ++i) {
        if (iteration == 0 || std::abs(e(i)) > params.θ * std::abs(old_e(i))) {
            Σ(i) = std::fmin(params.σₘₐₓ,
                             std::fmax(params.Δ * std::abs(e(i)) / norm_e, 1) *
                                 Σ(i));
        }
    }
    // std::cout << "[ALM]   Σ = " << Σ.transpose() << std::endl;
#else
    Σ *= params.Δ; // OpEn-style penalty update
    std::cout << "[ALM]   Σ = " << Σ(0) << std::endl;
#endif
}

void initialize_penalty(const Problem &p, const ALMParams &params,
                        const vec &x0, vec &Σ) {
    real_t f0 = p.f(x0);
    vec g0(p.m);
    p.g(x0, g0);
    // TODO: reuse evaluations of f ang g in PANOC
    real_t σ = params.σ₀ * std::abs(f0) / g0.squaredNorm();
    σ        = std::max(σ, params.σ₀);
    Σ.fill(σ);
}

} // namespace detail

ALMSolver::Stats ALMSolver::operator()(const Problem &problem, vec &y, vec &x) {
    auto start_time = std::chrono::steady_clock::now();

    auto sigNaN   = std::numeric_limits<real_t>::signaling_NaN();
    vec Σ         = vec::Constant(problem.m, sigNaN);
    vec z         = vec::Constant(problem.m, sigNaN);
    vec error     = vec::Constant(problem.m, sigNaN);
    vec error_old = vec::Constant(problem.m, sigNaN);

    Stats s;

    Problem prec_problem = problem;
    real_t prec_f;
    vec prec_g;

    if (params.preconditioning) {
        vec grad_f(problem.n);
        vec v = vec::Zero(problem.m);
        vec grad_g(problem.n);
        problem.grad_f(x, grad_f);
        prec_f = 1. / std::max(grad_f.lpNorm<Eigen::Infinity>(), 1.);
        prec_g.resize(problem.m);

        for (Eigen::Index i = 0; i < problem.m; ++i) {
            v(i) = 1;
            problem.grad_g(x, v, grad_g);
            v(i)      = 0;
            prec_g(i) = 1. / std::max(grad_g.lpNorm<Eigen::Infinity>(), 1.);
        }

        prec_problem.f = [&](const vec &x) { return problem.f(x) * prec_f; };
        prec_problem.grad_f = [&](const vec &x, vec &grad_f) {
            problem.grad_f(x, grad_f);
            grad_f *= prec_f;
        };
        prec_problem.g = [&](const vec &x, vec &g) {
            problem.g(x, g);
            g = prec_g.asDiagonal() * g;
        };
        prec_problem.grad_g = [&](const vec &x, const vec &v, vec &grad_g) {
            vec prec_v = prec_g.asDiagonal() * v;
            problem.grad_g(x, prec_v, grad_g);
        };
        prec_problem.D.lowerbound = prec_g.asDiagonal() * problem.D.lowerbound;
        prec_problem.D.upperbound = prec_g.asDiagonal() * problem.D.upperbound;

        std::cout << "prec_g: " << prec_g.transpose() << std::endl;
        std::cout << "prec_f: " << prec_f << std::endl;
    }
    const auto &p = params.preconditioning ? prec_problem : problem;

    // Initialize the penalty weights
    if (params.Σ₀ > 0) {
        Σ.fill(params.Σ₀);
    }
    // Initial penalty weights from problem
    else {
        detail::initialize_penalty(p, params, x, Σ);
    }

    real_t ε = params.ε₀;

    for (unsigned int i = 0; i < params.max_iter; ++i) {
#ifdef PRINT_DEBUG_COUT
        std::cout << std::endl;
        std::cout << "[ALM]   "
                  << "Iteration #" << i << std::endl;
#endif
        detail::project_y(y, p.D.lowerbound, p.D.upperbound, params.M);
        auto ps = panoc(p, x, z, y, error, Σ, ε);
        s.inner_iterations += ps.iterations;
        s.inner_linesearch_failures += ps.linesearch_failures;
        s.inner_lbfgs_failures += ps.lbfgs_failures;
        s.inner_lbfgs_rejected += ps.lbfgs_rejected;
        s.inner_convergence_failures += ps.status != SolverStatus::Converged;
        real_t norm_e = vec_util::norm_inf(error);

        if (params.print_interval != 0 && i % params.print_interval == 0) {
            std::cout << "[\x1b[0;34mALM\x1b[0m]   " << std::setw(5) << i
                      << ": ‖Σ‖ = " << std::setw(13) << Σ.norm()
                      << ", δ = " << std::setw(13) << norm_e
                      << ", ε = " << std::setw(13) << ps.ε << "\r\n";
        }

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        // TODO: check penalty size?
        if (ps.status == SolverStatus::NotFinite ||
            ps.status == SolverStatus::Interrupted) {
            s.ε                = ps.ε;
            s.δ                = norm_e;
            s.norm_penalty     = Σ.norm();
            s.outer_iterations = i;
            s.elapsed_time     = duration_cast<microseconds>(time_elapsed);
            s.status           = ps.status;
            if (params.preconditioning)
                y = prec_g.asDiagonal() * y / prec_f;
            return s;
        }

        bool converged = ps.ε <= params.ε &&
                         ps.status == SolverStatus::Converged &&
                         norm_e <= params.δ;
        bool out_of_time = time_elapsed > params.max_time;
        if (converged || i + 1 == params.max_iter || out_of_time) {
            s.ε                = ps.ε;
            s.δ                = norm_e;
            s.norm_penalty     = Σ.norm();
            s.outer_iterations = i + 1;
            s.elapsed_time     = duration_cast<microseconds>(time_elapsed);
            s.status           = converged     ? SolverStatus::Converged
                                 : out_of_time ? SolverStatus::MaxTime
                                               : SolverStatus::MaxIter;
            if (params.preconditioning)
                y = prec_g.asDiagonal() * y / prec_f;
            return s;
        }
        detail::update_penalty_weights(params, i, error, error_old, norm_e, Σ);
        ε = std::fmax(params.ρ * ε, params.ε);
        std::swap(error_old, error);
    }
    throw std::logic_error("[ALM]   loop error");
}

} // namespace pa