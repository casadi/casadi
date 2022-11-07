#pragma once

#include <alpaqa/config/config.hpp>
#include <Eigen/LU>

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
struct LQRFactor {
    USING_ALPAQA_CONFIG(Conf);

    using Dim = alpaqa::Dim<config_t>;

    static auto
    factor_masked_impl(Dim d,    ///< Problem dimensions
                       auto &&A, ///< System matrix A
                       auto &&B, ///< Input matrix B
                       auto &&Q, ///< State cost matrix Q
                       auto &&R, ///< Input cost matrix R
                       auto &&q, ///< Linear state factor q
                       auto &&r, ///< Linear input factor r
                       auto &&u, ///< Fixed inputs u
                       auto &&J, ///< Index set of inactive constraints
                       auto &&K  ///< Index set of active constraints
    ) {
        using mmat       = Eigen::Map<mat>;
        auto [N, nx, nu] = d;
        mat P{nx, nx};
        mat gain_K{nu, N * nx};
        mat e{nu, N};
        vec s{nx};
        vec c{nx};
        vec y{nx};
        vec t{nu};
        vec R̅_sto{nu * nu};
        vec S̅_sto{nu * nx};
        vec BiJ_sto{nx * nu};
        vec PBiJ_sto{nx * nu};
        mat PA{nx, nx};

        assign_possibly_diagonal(P, Q(N));
        s = q(N);
        for (index_t i = N; i-- > 0;) {
            auto &&Ai   = A(i);
            auto &&Bi   = B(i);
            auto &&Ri   = R(i);
            auto &&Qi   = Q(i);
            auto &&ui   = u(i);
            auto &&Ji   = J(i);
            auto &&Ki   = K(i);
            length_t nJ = Ji.size(); // number of inactive constraints
            mmat R̅{R̅_sto.data(), nJ, nJ};
            mmat S̅{S̅_sto.data(), nJ, nx};
            mmat BiJ{BiJ_sto.data(), nx, nJ};
            mmat PBiJ{PBiJ_sto.data(), nx, nJ};
            auto &&ti = t.topRows(nJ);
            // todo: make K contiguous
            auto &&gain_Ki = gain_K.middleCols(i * nx, nx).topRows(nJ);
            auto &&ei      = e.col(i).topRows(nJ);
            // R̅ ← R + Bᵀ P B
            BiJ.noalias()  = Bi(Eigen::all, Ji);
            PBiJ.noalias() = P * BiJ;
            R̅.noalias()    = BiJ.transpose() * PBiJ;
            add_possibly_diagonal_masked(R̅, Ri, Ji);
            // S̅ ← S + Bᵀ P A
            PA.noalias() = P * Ai;
            S̅.noalias()  = BiJ.transpose() * PA;
            // y ← Pc + s
            c = Bi(Eigen::all, Ki) * ui(Ki);
            y = P * c + s;
            // K ← -R̅⁻¹S̅
            Eigen::PartialPivLU<rmat> R̅LU{R̅};
            gain_Ki.noalias() = -R̅LU.solve(S̅);
            // e ← -R̅⁻¹(Bᵀy + r)
            ti = BiJ.transpose() * y + r(i)(Ji);
            if (Ri.cols() > 1 && Ri.rows() > 1)
                ti += Ri(Ji, Ki) * ui(Ki);
            ei = -R̅LU.solve(ti);
            if (i > 0) {
                // P ← Q + Aᵀ P A + S̅ᵀ K
                P.noalias() = Ai.transpose() * PA + S̅.transpose() * gain_Ki;
                // s ← S̅ᵀ e + Aᵀ y + q
                s = S̅.transpose() * ei + Ai.transpose() * y + q(i);
                add_possibly_diagonal(P, Qi);
            }
        }

        return std::make_tuple(std::move(gain_K), std::move(e));
    }

    static auto solve_masked_impl(Dim d, auto &&A, auto &&B, auto &&J,
                                  crmat Δu_eq, crmat K, crmat e) {
        auto [N, nx, nu] = d;
        mat Δx{nx, N + 1};
        mat Δu = Δu_eq.reshaped(nu, N);
        Δx.topRows(nx).setZero();
        for (index_t i = 0; i < N; ++i) {
            auto &&Ji     = J(i);
            auto &&Ai     = A(i);
            auto &&Bi     = B(i);
            length_t nJ   = Ji.size();
            auto &&Ki     = K.middleCols(i * nx, nx).topRows(nJ);
            auto &&ei     = e.col(i).topRows(nJ);
            auto &&Δxi    = Δx.col(i);
            auto &&Δui    = Δu.col(i);
            Δui(Ji)       = Ki * Δxi + ei;
            Δx.col(i + 1) = Ai * Δxi + Bi * Δui;
        }
        return std::make_tuple(std::move(Δx), std::move(Δu));
    }

    /// If `src` is a matrix, performs `dest = src`, if `src` is a vector,
    /// performs `dest = src.asDiagonal()`.
    static void assign_possibly_diagonal(rmat dest, crmat src);
    /// If `src` is a matrix, performs `dest += src`, if `src` is a vector,
    /// performs `dest += src.asDiagonal()`.
    static void add_possibly_diagonal(rmat dest, crmat src);
    /// If `src` is a matrix, performs `dest += src(V,V)`, if `src` is a vector,
    /// performs `dest += src(V).asDiagonal()`.
    static void add_possibly_diagonal_masked(rmat dest, crmat src,
                                             const auto &V);
    /// If `M` is a matrix, performs `dest = M * v`, if `M` is a vector,
    /// performs `dest = M.asDiagonal() * v`.
    static void mat_vec_possibly_diagonal(rvec dest, crmat M, crvec v);
};

template <Config Conf>
struct StatefulLQRFactor {
    USING_ALPAQA_CONFIG(Conf);

    using Dim = alpaqa::Dim<config_t>;

    StatefulLQRFactor(Dim d) : dim{d} {}
    Dim dim;
    mat P{dim.nx, dim.nx};
    mat gain_K{dim.nu, dim.N *dim.nx};
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

    void factor_masked(auto &&A, ///< System matrix A
                       auto &&B, ///< Input matrix B
                       auto &&Q, ///< State cost matrix Q
                       auto &&R, ///< Input cost matrix R
                       auto &&q, ///< Linear state factor q
                       auto &&r, ///< Linear input factor r
                       auto &&u, ///< Fixed inputs u
                       auto &&J, ///< Index set of inactive constraints
                       auto &&K  ///< Index set of active constraints
    ) {
        using mmat       = Eigen::Map<mat>;
        auto [N, nx, nu] = dim;

        assign_possibly_diagonal(P, Q(N));
        s = q(N);
        for (index_t i = N; i-- > 0;) {
            auto &&Ai   = A(i);
            auto &&Bi   = B(i);
            auto &&Ri   = R(i);
            auto &&Qi   = Q(i);
            auto &&ui   = u(i);
            auto &&Ji   = J(i);
            auto &&Ki   = K(i);
            length_t nJ = Ji.size(); // number of inactive constraints
            mmat R̅{R̅_sto.data(), nJ, nJ};
            mmat S̅{S̅_sto.data(), nJ, nx};
            mmat BiJ{BiJ_sto.data(), nx, nJ};
            mmat PBiJ{PBiJ_sto.data(), nx, nJ};
            auto &&ti = t.topRows(nJ);
            // todo: make K contiguous
            auto &&gain_Ki = gain_K.middleCols(i * nx, nx).topRows(nJ);
            auto &&ei      = e.col(i).topRows(nJ);
            // R̅ ← R + Bᵀ P B
            BiJ.noalias()  = Bi(Eigen::all, Ji);
            PBiJ.noalias() = P * BiJ;
            R̅.noalias()    = BiJ.transpose() * PBiJ;
            add_possibly_diagonal_masked(R̅, Ri, Ji);
            // S̅ ← S + Bᵀ P A
            PA.noalias() = P * Ai;
            S̅.noalias()  = BiJ.transpose() * PA;
            // y ← Pc + s
            c = Bi(Eigen::all, Ki) * ui(Ki);
            y = P * c + s;
            // K ← -R̅⁻¹S̅
            Eigen::PartialPivLU<rmat> R̅LU{R̅};
            gain_Ki.noalias() = -R̅LU.solve(S̅);
            // e ← -R̅⁻¹(Bᵀy + r)
            ti = BiJ.transpose() * y + r(i)(Ji);
            if (Ri.cols() > 1 && Ri.rows() > 1)
                ti += Ri(Ji, Ki) * ui(Ki);
            ei = -R̅LU.solve(ti);
            if (i > 0) {
                // P ← Q + Aᵀ P A + S̅ᵀ K
                P.noalias() = Ai.transpose() * PA + S̅.transpose() * gain_Ki;
                // s ← S̅ᵀ e + Aᵀ y + q
                s = S̅.transpose() * ei + Ai.transpose() * y + q(i);
                add_possibly_diagonal(P, Qi);
            }
        }
    }

    void solve_masked(auto &&A, auto &&B, auto &&J, rvec Δxu_eq) {
        auto [N, nx, nu] = dim;
        for (index_t i = 0; i < N; ++i) {
            auto &&Ji      = J(i);
            length_t nJ    = Ji.size();
            auto &&Ki      = gain_K.middleCols(i * nx, nx).topRows(nJ);
            auto &&ei      = e.col(i).topRows(nJ);
            auto &&Δui     = Δxu_eq.segment(i * (nx + nu) + nx, nu);
            auto &&Δxi     = Δxu_eq.segment(i * (nx + nu), nx);
            auto &&Δx_next = Δxu_eq.segment((i + 1) * (nx + nu), nx);
            Δui(Ji)        = Ki * Δxi + ei;
            Δx_next        = A(i) * Δxi + B(i) * Δui;
        }
    }

    /// If `src` is a matrix, performs `dest = src`, if `src` is a vector,
    /// performs `dest = src.asDiagonal()`.
    static void assign_possibly_diagonal(rmat dest, crmat src);
    /// If `src` is a matrix, performs `dest += src`, if `src` is a vector,
    /// performs `dest += src.asDiagonal()`.
    static void add_possibly_diagonal(rmat dest, crmat src);
    /// If `src` is a matrix, performs `dest += src(V,V)`, if `src` is a vector,
    /// performs `dest += src(V).asDiagonal()`.
    static void add_possibly_diagonal_masked(rmat dest, crmat src,
                                             const auto &V);
    /// If `M` is a matrix, performs `dest = M * v`, if `M` is a vector,
    /// performs `dest = M.asDiagonal() * v`.
    static void mat_vec_possibly_diagonal(rvec dest, crmat M, crvec v);
};

} // namespace alpaqa

#include <alpaqa/inner/directions/panoc-ocp/src/lqr.tpp>