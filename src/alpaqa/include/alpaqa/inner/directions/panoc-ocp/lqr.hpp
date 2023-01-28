#pragma once

#include <alpaqa/config/config.hpp>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <cassert>

namespace alpaqa {

template <Config Conf>
struct Dim {
    USING_ALPAQA_CONFIG(Conf);
    length_t N, nx, nu;
    struct Horizon {
        length_t N;
        struct Iter {
            index_t i;
            Iter &operator++() {
                ++i;
                return *this;
            }
            Iter operator++(int) const {
                Iter t = *this;
                ++i;
                return t;
            }
            index_t &operator*() { return i; }
            const index_t &operator*() const { return i; }
            friend auto operator<=>(const Iter &, const Iter &) = default;
        };
        static Iter begin() { return {0}; }
        Iter end() const { return {N}; }
    };
    Horizon horizon() const { return {N}; }
};

template <Config Conf>
struct StatefulLQRFactor {
    USING_ALPAQA_CONFIG(Conf);

    using Dim = alpaqa::Dim<config_t>;

    StatefulLQRFactor(Dim d) : dim{d} {}
    Dim dim;
    mat P{dim.nx, dim.nx};
    mat gain_K{dim.nu * dim.nx, dim.N};
    mat e{dim.nu, dim.N};
    vec s{dim.nx};
    vec c{dim.nx};
    vec y{dim.nx};
    vec t{dim.nu};
    vec R̅_sto{dim.nu * dim.nu};
    vec S̅_sto{dim.nu * dim.nx};
    vec BiJ_sto{dim.nx * dim.nu};
    vec PBiJ_sto{dim.nx * dim.nu};
    mat PA{dim.nx, dim.nx};
    real_t min_rcond = 1;

    void factor_masked(auto &&AB,        ///< System matrix A & input matrix B
                       auto &&Q,         ///< State cost matrix Q
                       auto &&R,         ///< Input cost matrix R
                       auto &&S,         ///< Cross cost matrix S
                       auto &&R_prod,    ///< Product with input cost matrix R
                       auto &&S_prod,    ///< Product with cross cost matrix S
                       auto &&q,         ///< Linear state factor q
                       auto &&r,         ///< Linear input factor r
                       auto &&u,         ///< Fixed inputs u
                       auto &&J,         ///< Index set of inactive constraints
                       auto &&K,         ///< Index set of active constraints
                       bool use_cholesky ///< Use Cholesky instead of LU solver
    ) {
        using mmat = Eigen::Map<mat>;
        using Eigen::indexing::all;
        auto [N, nx, nu] = dim;

        min_rcond = 1;
        P.setZero();
        Q(N)(P);
        s = q(N);
        for (index_t i = N; i-- > 0;) {
            auto &&ABi  = AB(i);
            auto &&Ai   = ABi.leftCols(nx);
            auto &&Bi   = ABi.rightCols(nu);
            auto &&ui   = u(i);
            auto &&Ji   = J(i);
            auto &&Ki   = K(i);
            length_t nJ = Ji.size(); // number of inactive constraints
            mmat R̅{R̅_sto.data(), nJ, nJ};
            mmat S̅{S̅_sto.data(), nJ, nx};
            mmat BiJ{BiJ_sto.data(), nx, nJ};
            mmat PBiJ{PBiJ_sto.data(), nx, nJ};
            auto &&ti = t.topRows(nJ);
            mmat gain_Ki{gain_K.col(i).data(), nJ, nx};
            auto &&ei = e.col(i).topRows(nJ);
            // R̅ ← R + Bᵀ P B
            BiJ.noalias()  = Bi(all, Ji);
            PBiJ.noalias() = P * BiJ;
            R̅.noalias()    = BiJ.transpose() * PBiJ;
            R(i)(Ji, R̅);
            // S̅ ← S + Bᵀ P A
            PA.noalias() = P * Ai;
            S̅.noalias()  = BiJ.transpose() * PA;
            S(i)(Ji, S̅);
            // c = B(·,K) u(K), y ← P c + s
            c.noalias() = Bi(all, Ki) * ui(Ki);
            y.noalias() = P * c;
            y += s;
            // t ← Bᵀy + r + R(J,K) u(K)
            ti.noalias() = BiJ.transpose() * y;
            ti += r(i)(Ji);
            R_prod(i)(Ji, Ki, ui, ti);
            // Factor R̅
            if (use_cholesky) {
#ifdef EIGEN_RUNTIME_NO_MALLOC
                bool prev = Eigen::internal::is_malloc_allowed();
                Eigen::internal::set_is_malloc_allowed(true); // TODO
#endif
                Eigen::LDLT<rmat> R̅LU{R̅};
                min_rcond = std::min(R̅LU.rcond(), min_rcond);
#ifdef EIGEN_RUNTIME_NO_MALLOC
                Eigen::internal::set_is_malloc_allowed(prev);
#endif
                // K ← -R̅⁻¹S̅
                gain_Ki.noalias() = R̅LU.solve(S̅);
                // e ← -R̅⁻¹(Bᵀy + r)
                ei.noalias() = R̅LU.solve(ti);
            } else {
#ifdef EIGEN_RUNTIME_NO_MALLOC
                bool prev = Eigen::internal::is_malloc_allowed();
                Eigen::internal::set_is_malloc_allowed(true); // TODO
#endif
                Eigen::PartialPivLU<rmat> R̅LU{R̅};
                min_rcond = std::min(R̅LU.rcond(), min_rcond);
#ifdef EIGEN_RUNTIME_NO_MALLOC
                Eigen::internal::set_is_malloc_allowed(prev);
#endif
                // K ← -R̅⁻¹S̅
                gain_Ki.noalias() = R̅LU.solve(S̅);
                // e ← -R̅⁻¹(Bᵀy + r)
                ei.noalias() = R̅LU.solve(ti);
            }
            gain_Ki = -gain_Ki;
            ei      = -ei;
            if (i > 0) {
                // P ← Q + Aᵀ P A + S̅ᵀ K
                P.noalias() = Ai.transpose() * PA;
                P.noalias() += S̅.transpose() * gain_Ki;
                // s ← S̅ᵀ e + Aᵀ y + q + Sᵀ(·,K) u(K)
                s.noalias() = S̅.transpose() * ei;
                s.noalias() += Ai.transpose() * y;
                s += q(i);
                S_prod(i)(Ki, ui, s);
                Q(i)(P);
            }
        }
    }

    void solve_masked(auto &&AB, auto &&J, rvec Δu_eq, rvec Δx) {
        auto [N, nx, nu] = dim;
        assert(Δx.size() == 2 * nx);
        Δx.topRows(nx).setZero();
        for (index_t i = 0; i < N; ++i) {
            auto &&ABi     = AB(i);
            auto &&Ai      = ABi.leftCols(nx);
            auto &&Bi      = ABi.rightCols(nu);
            auto &&Ji      = J(i);
            auto &&Δxi     = Δx.segment((i % 2) * nx, nx);
            auto &&Δx_next = Δx.segment(((i + 1) % 2) * nx, nx);
            length_t nJ    = Ji.size();
            mmat Ki{gain_K.col(i).data(), nJ, nx};
            auto &&ei  = e.col(i).topRows(nJ);
            auto &&Δui = Δu_eq.segment(i * nu, nu);
            ei.noalias() += Ki * Δxi;
            Δui(Ji).noalias() = ei;
            Δx_next.noalias() = Ai * Δxi;
            Δx_next.noalias() += Bi * Δui;
        }
    }
};

} // namespace alpaqa
