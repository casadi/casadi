#include <panoc-alm/inner/detail/panoc-helpers.hpp>
#include <panoc-alm/inner/directions/lbfgs.hpp>
#include <panoc-alm/interop/casadi/CasADiLoader.hpp>
#include <panoc-alm/interop/cutest/CUTEstLoader.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>

#include <Eigen/Cholesky>

int main() {
    using namespace pa;

    const auto so_name = "examples/experiments/libhimmelblau_functions.so";
    Problem p;
    p.n = 2;
    p.m = 2;
    p.C = Box{
        vec::Constant(p.n, inf),
        vec::Constant(p.n, 2.1),
    };
    p.D = Box{
        vec::Constant(p.m, 2.),
        vec::Constant(p.m, 0.),
    };
    p.f           = load_CasADi_objective(so_name);
    p.grad_f      = load_CasADi_gradient_objective(so_name);
    p.g           = load_CasADi_constraints(so_name);
    p.grad_g_prod = load_CasADi_gradient_constraints_prod(so_name);
    p.hess_L      = load_CasADi_hessian_lagrangian(so_name);

    LBFGSParams lbfgsparam;
    lbfgsparam.memory = 5;
    LBFGS lbfgs(lbfgsparam, p.n);

    mat H(p.n, p.n);
    vec x(p.n), y(p.m), g(p.m), Σ(p.m), ŷ(p.m), grad_f(p.n), grad_g(p.n);
    x(0) = 3.3;
    x(1) = 3;
    y.fill(0);
    Σ.fill(1e-16);
    real_t γ = 5e-3;

    std::cout << std::setprecision(17);

    bool gradient_descent = true;
    bool use_lbfgs        = true;

    for (size_t loop_idx = 0; loop_idx < 20; ++loop_idx) {

        std::cout << x(0) << ", " << x(1) << "," << std::endl;
        [[maybe_unused]] real_t ψ = detail::calc_ψ_ŷ(p, x, y, Σ, ŷ);
        // std::cout << ψ << std::endl;

        p.hess_L(x, ŷ, H);
        p.grad_f(x, grad_f);

        p.g(x, g);
        for (vec::Index i = 0; i < p.m; ++i) {
            real_t ζ      = g(i) - ŷ(i) / Σ(i);
            bool inactive = p.D.lowerbound(i) < ζ && ζ < p.D.upperbound(i);
            vec ei        = vec::Zero(p.m);
            ei(i)         = 1;
            p.grad_g_prod(x, ei, grad_g);
            if (inactive)
                H += Σ(i) * grad_g * grad_g.transpose();
        }

        vec grad_L(p.n), work_n(p.n);
        detail::calc_grad_ψ_from_ŷ(p, x, ŷ, grad_L, work_n);
        // std::cout << "∇L = " << grad_L.transpose() << std::endl;
        // p.grad_f(x, grad_f);
        // std::cout << "∇f = " << grad_f.transpose() << std::endl;

        using indexvec  = std::vector<vec::Index>;
        using bvec      = Eigen::Matrix<bool, Eigen::Dynamic, 1>;
        bvec inactive_C = bvec::Zero(p.n);
        vec R(p.n);

        indexvec J, K;
        J.reserve(p.n);
        K.reserve(p.n);
        for (vec::Index i = 0; i < p.n; ++i) {
            real_t gd_step = x(i) - γ * grad_L(i);
            if (gd_step < p.C.lowerbound(i)) {
                R(i) = x(i) - p.C.lowerbound(i);
                K.push_back(i);
            } else if (p.C.upperbound(i) < gd_step) {
                R(i) = x(i) - p.C.upperbound(i);
                K.push_back(i);
            } else {
                R(i)          = γ * grad_L(i);
                inactive_C(i) = true;
                J.push_back(i);
            }
        }
        vec d = -R;
        if (not J.empty() && not gradient_descent) {
            vec rhs(J.size());
            vec::Index ri = 0;
            for (auto j : J) {
                real_t hess_d = 0;
                for (auto k : K)
                    hess_d += H(j, k) * d(k);
                rhs(ri++) = -grad_L(j) - hess_d;
            }
            mat hess_Ljj(J.size(), J.size());

            vec::Index r = 0;
            for (auto rj : J) {
                vec::Index c = 0;
                for (auto cj : J) {
                    hess_Ljj(r, c) = H(rj, cj);
                    ++c;
                }
                ++r;
            }

            auto ldl = hess_Ljj.ldlt();
            vec dj   = ldl.solve(rhs);

            r = 0;
            for (auto j : J)
                d(j) = dj(r++);
        }

        if (use_lbfgs) {
            lbfgs.apply(d, 1);
            vec pₖ   = detail::projected_gradient_step(p.C, γ, x, grad_L);
            vec xₖ₊₁ = x + d;
            vec grad_Lₖ₊₁(p.n);
            detail::calc_grad_ψ(p, xₖ₊₁, y, Σ, grad_Lₖ₊₁, work_n, ŷ);
            vec pₖ₊₁ = detail::projected_gradient_step(p.C, γ, xₖ₊₁, grad_Lₖ₊₁);
            lbfgs.update(x, xₖ₊₁, pₖ, pₖ₊₁, LBFGS::Sign::Negative);
        }

        x += d;
    }
    std::cout << x(0) << ", " << x(1) << "," << std::endl;
    vec grad_L(p.n), work_n(p.n), work_m(p.n);
    detail::calc_grad_ψ(p, x, y, Σ, grad_L, work_n, work_m);
    std::cout << grad_L.transpose() << std::endl;

    {
        mat H(p.n, p.n);
        mat Hh(p.n + 5, p.n + 3);
        Hh << 11, 12, 13, 14, 15, //
            21, 22, 23, 24, 25,   //
            31, 32, 33, 34, 35,   //
            41, 42, 43, 44, 45,   //
            51, 52, 53, 54, 55,   //
            61, 62, 63, 64, 65,   //
            71, 72, 73, 74, 75;
        p.hess_L(x, y, H);
        p.hess_L(x, y, Hh);
        std::cout << H << std::endl << std::endl;
        std::cout << Hh << std::endl << std::endl;
    }
}