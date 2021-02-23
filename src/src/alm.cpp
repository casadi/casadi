#include <cfenv>
#include <panoc-alm/alm.hpp>

namespace pa {
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

} // namespace detail

ALMSolver::Stats ALMSolver::operator()(const Problem &problem, vec &y, vec &x) {
    // feenableexcept(FE_INVALID); // TODO
    auto sigNaN   = std::numeric_limits<real_t>::signaling_NaN();
    vec Σ         = vec::Constant(problem.m, sigNaN);
    vec z         = vec::Constant(problem.m, sigNaN);
    vec error     = vec::Constant(problem.m, sigNaN);
    vec error_old = vec::Constant(problem.m, sigNaN);

    Stats s;

    // Initialize the penalty weights
    Σ.fill(params.Σ₀);
    real_t ε = params.ε₀;

    for (unsigned int i = 0; i < params.max_iter; ++i) {
#ifdef PRINT_DEBUG_COUT
        std::cout << std::endl;
        std::cout << "[ALM]   "
                  << "Iteration #" << i << std::endl;
#endif
        detail::project_y(y, problem.D.lowerbound, problem.D.upperbound,
                          params.M);
        auto ps = panoc(problem, x, z, y, error, Σ, ε);
        s.inner_iterations += ps.iterations;
        real_t norm_e = vec_util::norm_inf(error);

        if (ε <= params.ε && norm_e <= params.δ) {
            s.outer_iterations = i;
            return s;
        }
        detail::update_penalty_weights(params, i, error, error_old, norm_e, Σ);
        ε = std::fmax(params.ρ * ε, params.ε);
        std::swap(error_old, error);
    }
    throw std::runtime_error("[ALM]   Maximum iterations exceeded");
}

} // namespace pa