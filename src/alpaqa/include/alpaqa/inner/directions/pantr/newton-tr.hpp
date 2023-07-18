#pragma once

#include <alpaqa/accelerators/steihaugcg.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <alpaqa/util/alloc-check.hpp>
#include <alpaqa/util/index-set.hpp>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace alpaqa {

/// Parameters for the @ref NewtonTRDirection class.
template <Config Conf>
struct NewtonTRDirectionParams {
    USING_ALPAQA_CONFIG(Conf);
    bool rescale_on_step_size_changes = false;
    real_t hessian_vec_factor         = real_t(0.5);
    bool finite_diff                  = false;
    real_t finite_diff_stepsize =
        std::sqrt(std::numeric_limits<real_t>::epsilon());
};

/// @ingroup grp_DirectionProviders
template <Config Conf>
struct NewtonTRDirection {
    USING_ALPAQA_CONFIG(Conf);

    using Problem           = TypeErasedProblem<config_t>;
    using AcceleratorParams = SteihaugCGParams<config_t>;
    using DirectionParams   = NewtonTRDirectionParams<config_t>;

    struct Params {
        AcceleratorParams accelerator = {};
        DirectionParams direction     = {};
    };

    NewtonTRDirection() = default;
    NewtonTRDirection(const Params &params)
        : steihaug(params.accelerator), direction_params(params.direction) {}
    NewtonTRDirection(const AcceleratorParams &params,
                      const DirectionParams &directionparams = {})
        : steihaug(params), direction_params(directionparams) {}

    /// @see @ref PANTRDirection::initialize
    void initialize(const Problem &problem, [[maybe_unused]] crvec y,
                    [[maybe_unused]] crvec Σ, [[maybe_unused]] real_t γ_0,
                    [[maybe_unused]] crvec x_0, [[maybe_unused]] crvec x̂_0,
                    [[maybe_unused]] crvec p_0,
                    [[maybe_unused]] crvec grad_ψx_0) {
        if (!direction_params.finite_diff &&
            !problem.provides_eval_hess_ψ_prod())
            throw std::invalid_argument("NewtonTR without finite differences "
                                        "requires Problem::eval_hess_ψ_prod()");
        if (!problem.provides_eval_inactive_indices_res_lna())
            throw std::invalid_argument(
                "NewtonTR requires "
                "Problem::eval_inactive_indices_res_lna()");
        // Store references to problem and ALM variables
        this->problem = &problem;
        this->y.emplace(y);
        this->Σ.emplace(Σ);
        // Resize workspaces
        const auto n = problem.get_n(), m = problem.get_m();
        JK_sto.resize(n);
        rJ_sto.resize(n);
        qJ_sto.resize(n);
        work.resize(n);
        work_2.resize(n);
        steihaug.resize(n);
        if (direction_params.finite_diff) {
            work_n_fd.resize(n);
            work_m_fd.resize(m);
        }
    }

    /// @see @ref PANTRDirection::has_initial_direction
    bool has_initial_direction() const { return true; }

    /// @see @ref PANTRDirection::update
    bool update([[maybe_unused]] real_t γₖ, [[maybe_unused]] real_t γₙₑₓₜ,
                [[maybe_unused]] crvec xₖ, [[maybe_unused]] crvec xₙₑₓₜ,
                [[maybe_unused]] crvec pₖ, [[maybe_unused]] crvec pₙₑₓₜ,
                [[maybe_unused]] crvec grad_ψxₖ,
                [[maybe_unused]] crvec grad_ψxₙₑₓₜ) {
        return true;
    }

    /// @see @ref PANTRDirection::apply
    real_t apply([[maybe_unused]] real_t γₖ, [[maybe_unused]] crvec xₖ,
                 [[maybe_unused]] crvec x̂ₖ, crvec pₖ,
                 [[maybe_unused]] crvec grad_ψxₖ, real_t radius,
                 rvec qₖ) const {

        if (!std::isfinite(radius))
            throw std::logic_error("Invalid trust radius");
        if (radius < std::numeric_limits<real_t>::epsilon())
            throw std::logic_error("Trust radius too small");

        // Newton with exact Hessian

        // Find inactive and active constraints
        const auto n = problem->get_n();
        index_t nJ =
            problem->eval_inactive_indices_res_lna(γₖ, xₖ, grad_ψxₖ, JK_sto);
        crindexvec J = JK_sto.topRows(nJ);
        rindexvec K  = JK_sto.bottomRows(n - nJ);
        detail::IndexSet<config_t>::compute_complement(J, K, n);
        auto rJ = rJ_sto.topRows(nJ);
        auto qJ = qJ_sto.topRows(nJ);
        rJ      = (-real_t(1) / γₖ) * pₖ(J);
        qₖ(K)   = pₖ(K);
        qₖ(J).setZero();
        real_t norm_qK_sq = pₖ(K).squaredNorm();

        // Hessian-vector term
        if (direction_params.hessian_vec_factor != 0) {
            if (direction_params.finite_diff) {
                real_t ε = (1 + grad_ψxₖ.norm()) *
                           direction_params.finite_diff_stepsize;
                work = xₖ + ε * qₖ;
                problem->eval_grad_ψ(work, *y, *Σ, work_2, work_n_fd,
                                     work_m_fd);
                rJ.noalias() += (work_2 - grad_ψxₖ)(J) *
                                (direction_params.hessian_vec_factor / ε);
            } else {
                problem->eval_hess_ψ_prod(xₖ, *y, *Σ, 1, qₖ, work);
                rJ.noalias() += work(J) * direction_params.hessian_vec_factor;
            }
        }

        // Hessian-vector product on subset J
        auto hess_vec_mult = [&](crvec p, rvec Bp) {
            if (direction_params.finite_diff) {
                real_t ε = (1 + grad_ψxₖ.norm()) *
                           direction_params.finite_diff_stepsize;
                work = xₖ;
                work(J) += ε * p;
                problem->eval_grad_ψ(work, *y, *Σ, work_2, work_n_fd,
                                     work_m_fd);
                Bp.topRows(nJ) = (work_2 - grad_ψxₖ)(J) / ε;
            } else {
                work.setZero();
                work(J) = p;
                problem->eval_hess_ψ_prod(xₖ, *y, *Σ, 1, work, work_2);
                Bp.topRows(nJ) = work_2(J);
            }
        };

        // Steihaug conjugate gradients
        real_t qJ_model = steihaug.solve(rJ, hess_vec_mult, radius, qJ);
        qₖ(J)           = qJ;
        return qJ_model - norm_qK_sq / (2 * γₖ);
    }

    /// @see @ref PANTRDirection::changed_γ
    void changed_γ([[maybe_unused]] real_t γₖ, [[maybe_unused]] real_t old_γₖ) {
        if (direction_params.rescale_on_step_size_changes)
            throw std::invalid_argument("NewtonTRDirection does not support "
                                        "rescale_on_step_size_changes");
    }

    /// @see @ref PANTRDirection::reset
    void reset() {}

    /// @see @ref PANTRDirection::get_name
    std::string get_name() const {
        return "NewtonTRDirection<" + std::string(config_t::get_name()) + '>';
    }

    auto get_params() const {
        return std::tie(steihaug.params, direction_params);
    }

    SteihaugCG<config_t> steihaug;
    DirectionParams direction_params;
    const Problem *problem = nullptr;
#ifndef _WIN32
    std::optional<crvec> y = std::nullopt;
    std::optional<crvec> Σ = std::nullopt;
#else
    std::optional<vec> y = std::nullopt;
    std::optional<vec> Σ = std::nullopt;
#endif
    mutable indexvec JK_sto;
    mutable vec rJ_sto;
    mutable vec qJ_sto;
    mutable vec work, work_2, work_n_fd, work_m_fd;
};

} // namespace alpaqa