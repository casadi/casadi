#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/alloc-check.hpp>

namespace alpaqa {

template <Config Conf>
struct SteihaugCGParams {
    USING_ALPAQA_CONFIG(Conf);
    real_t tol_scale       = 1;
    real_t tol_scale_root  = real_t(0.5);
    real_t tol_max         = inf<config_t>;
    real_t max_iter_factor = 1;
};

/// Steihaug conjugate gradients procedure based on
/// https://github.com/scipy/scipy/blob/583e70a50573169fc352b5dc6d94588a97c7389a/scipy/optimize/_trustregion_ncg.py#L44
template <Config Conf>
struct SteihaugCG {
    USING_ALPAQA_CONFIG(Conf);

    using Params = SteihaugCGParams<config_t>;
    Params params;

    SteihaugCG() = default;
    SteihaugCG(const Params &params) : params{params} {}

    mutable vec z, r, d, Bd, work_eval;

    void resize(length_t n) {
        z.resize(n);
        r.resize(n);
        d.resize(n);
        Bd.resize(n);
        work_eval.resize(n);
    }

    template <class HessFun>
    real_t solve(const auto &grad, const HessFun &hess_prod,
                 real_t trust_radius, rvec step) const {
        length_t n = grad.size();
        // get the norm of jacobian and define the origin
        auto v = [n](auto &v) { return v.topRows(n); };
        auto z = v(this->z), r = v(this->r), d = v(this->d), Bd = v(this->Bd);
        auto g = v(grad);
        auto s = v(step);
        // init the state for the first iteration
        z.setZero();
        r               = g;
        d               = -r;
        real_t r_sq     = r.squaredNorm();
        real_t grad_mag = g.norm();

        // define a default tolerance
        real_t tolerance =
            std::fmin(params.tol_max, params.tol_scale * grad_mag *
                                          std::fmin(params.tol_scale_root,
                                                    std::sqrt(grad_mag)));

        // Workspaces and function evaluation
        auto eval = [&](crvec p) {
            hess_prod(p, work_eval);
            return p.dot(g) + real_t(0.5) * p.dot(v(work_eval));
        };

        // Search for the min of the approximation of the objective function.
        index_t i           = 0;
        const auto max_iter = static_cast<index_t>(
            std::round(static_cast<real_t>(n) * params.max_iter_factor));
        while (true) {
            // do an iteration
            hess_prod(d, Bd);
            real_t dBd = d.dot(Bd);
            if (dBd <= 0) {
                // Look at the two boundary points.
                // Find both values of t to get the boundary points such that
                // ||z + t d|| == trust_radius
                // and then choose the one with the predicted min value.
                auto [ta, tb] =
                    get_boundaries_intersections(z, d, trust_radius);
                auto &pa   = r; // Reuse storage
                auto &pb   = d; // Reuse storage
                pa         = z + ta * d;
                pb         = z + tb * d;
                real_t q_a = eval(pa), q_b = eval(pb);
                real_t q_min = std::fmin(q_a, q_b);
                if (q_a == q_min) {
                    s = pa;
                    return q_a;
                } else {
                    s = pb;
                    return q_b;
                }
            }

            real_t alpha = r_sq / dBd;
            s            = z + alpha * d;
            if (s.norm() >= trust_radius) {
                // Find t >= 0 to get the boundary point such that
                // ||z + t d|| == trust_radius
                auto [ta, tb] =
                    get_boundaries_intersections(z, d, trust_radius);
                s = z + tb * d;
                return eval(s);
            }
            r += alpha * Bd;
            real_t r_next_sq = r.squaredNorm();
            if (std::sqrt(r_next_sq) < tolerance || i > max_iter)
                return eval(s);
            real_t beta_next = r_next_sq / r_sq;
            r_sq             = r_next_sq;
            d                = beta_next * d - r;
            z                = s;
            ++i;
        }
    }

    /// Solve the scalar quadratic equation ||z + t d|| == trust_radius.
    /// This is like a line-sphere intersection.
    /// Return the two values of t, sorted from low to high.
    static auto get_boundaries_intersections(crvec z, crvec d,
                                             real_t trust_radius) {
        real_t a = d.squaredNorm();
        real_t b = 2 * z.dot(d);
        real_t c = z.squaredNorm() - trust_radius * trust_radius;
        real_t sqrt_discriminant = std::sqrt(b * b - 4 * a * c);

        // The following calculation is mathematically
        // equivalent to:
        // ta = (-b - sqrt_discriminant) / (2*a)
        // tb = (-b + sqrt_discriminant) / (2*a)
        // but produce smaller round off errors.
        // Look at Matrix Computation p.97
        // for a better justification.
        real_t aux = b + std::copysign(sqrt_discriminant, b);
        real_t ta  = -aux / (2 * a);
        real_t tb  = -2 * c / aux;
        return std::make_tuple(std::fmin(ta, tb), std::fmax(ta, tb));
    }
};

} // namespace alpaqa