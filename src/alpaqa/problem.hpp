#pragma once

#include <alpaqa/util/problem.hpp>

inline auto prob_getter_f() {
    return [](const alpaqa::Problem &p) -> std::function<alpaqa::real_t(alpaqa::crvec)> {
        return [n{p.n}, f{p.f}](alpaqa::crvec x) {
            if (x.size() != n)
                throw std::out_of_range("Dimension of x not consistent "
                                        "with problem dimension n");
            return f(x);
        };
    };
}
inline auto prob_setter_f() {
    return [](alpaqa::Problem &p,
              std::function<alpaqa::real_t(alpaqa::crvec)> fun) -> void { p.f = fun; };
}
inline auto prob_getter_grad_f() {
    return [](const alpaqa::Problem &p) -> std::function<alpaqa::vec(alpaqa::crvec)> {
        return [n{p.n}, grad_f{p.grad_f}](alpaqa::crvec x) {
            if (x.size() != n)
                throw std::out_of_range("Dimension of x not consistent "
                                        "with problem dimension n");
            alpaqa::vec gr(n);
            grad_f(x, gr);
            return gr;
        };
    };
}
inline auto prob_setter_grad_f() {
    return [](alpaqa::Problem &p, std::function<alpaqa::vec(alpaqa::crvec)> fun) -> void {
        p.grad_f = [n{p.n}, fun{std::move(fun)}](alpaqa::crvec x, alpaqa::rvec gr) {
            auto &&res = fun(x);
            if (res.size() != n)
                throw std::out_of_range(
                    "Dimension of result of grad_f not consistent "
                    "with problem dimension n");
            gr = std::move(res);
        };
    };
}
inline auto prob_getter_g() {
    return [](const alpaqa::Problem &p) -> std::function<alpaqa::vec(alpaqa::crvec)> {
        return [n{p.n}, m{p.m}, g{p.g}](alpaqa::crvec x) {
            if (x.size() != n)
                throw std::out_of_range("Dimension of x not consistent "
                                        "with problem dimension n");
            alpaqa::vec gg(m);
            g(x, gg);
            return gg;
        };
    };
}
inline auto prob_setter_g() {
    return [](alpaqa::Problem &p, std::function<alpaqa::vec(alpaqa::crvec)> fun) -> void {
        p.g = [m{p.m}, fun{std::move(fun)}](alpaqa::crvec x, alpaqa::rvec gg) {
            auto &&res = fun(x);
            if (res.size() != m)
                throw std::out_of_range(
                    "Dimension of result of g not consistent "
                    "with problem dimension m");
            gg = std::move(res);
        };
    };
}
inline auto prob_getter_grad_g_prod() {
    return [](const alpaqa::Problem &p)
               -> std::function<alpaqa::vec(alpaqa::crvec, alpaqa::crvec)> {
        return [n{p.n}, m{p.m}, grad_g_prod{p.grad_g_prod}](alpaqa::crvec x,
                                                            alpaqa::crvec y) {
            if (x.size() != n)
                throw std::out_of_range("Dimension of x not consistent "
                                        "with problem dimension n");
            if (y.size() != m)
                throw std::out_of_range("Dimension of y not consistent "
                                        "with problem dimension m");
            alpaqa::vec gy(n);
            grad_g_prod(x, y, gy);
            return gy;
        };
    };
}
inline auto prob_setter_grad_g_prod() {
    return [](alpaqa::Problem &p,
              std::function<alpaqa::vec(alpaqa::crvec, alpaqa::crvec)> fun) -> void {
        p.grad_g_prod = [n{p.n}, fun{std::move(fun)}](alpaqa::crvec x, alpaqa::crvec y,
                                                      alpaqa::rvec gy) {
            auto &&res = fun(x, y);
            if (res.size() != n)
                throw std::out_of_range(
                    "Dimension of result of grad_g_prod not consistent "
                    "with problem dimension n");
            gy = std::move(res);
        };
    };
}
inline auto prob_getter_grad_gi() {
    return [](const alpaqa::Problem &p)
               -> std::function<alpaqa::vec(alpaqa::crvec, unsigned)> {
        return [n{p.n}, m{p.m}, grad_gi{p.grad_gi}](alpaqa::crvec x, unsigned i) {
            if (x.size() != n)
                throw std::out_of_range("Dimension of x not consistent "
                                        "with problem dimension n");
            if (i < m)
                throw std::out_of_range("Constraint index greater or "
                                        "equal to problem dimension m");
            alpaqa::vec gg(n);
            grad_gi(x, i, gg);
            return gg;
        };
    };
}
inline auto prob_setter_grad_gi() {
    return [](alpaqa::Problem &p,
              std::function<alpaqa::vec(alpaqa::crvec, unsigned)> fun) -> void {
        p.grad_gi = [n{p.n}, fun{std::move(fun)}](alpaqa::crvec x, unsigned i,
                                                  alpaqa::rvec gg) {
            auto &&res = fun(x, i);
            if (res.size() != n)
                throw std::out_of_range(
                    "Dimension of result of grad_gi not consistent "
                    "with problem dimension n");
            gg = std::move(res);
        };
    };
}
inline auto prob_getter_hess_L() {
    return [](const alpaqa::Problem &p)
               -> std::function<alpaqa::mat(alpaqa::crvec, alpaqa::crvec)> {
        return [n{p.n}, m{p.m}, hess_L{p.hess_L}](alpaqa::crvec x, alpaqa::crvec y) {
            if (x.size() != n)
                throw std::out_of_range("Dimension of x not consistent "
                                        "with problem dimension n");
            if (y.size() != m)
                throw std::out_of_range("Dimension of y not consistent "
                                        "with problem dimension m");
            alpaqa::mat H(n, n);
            hess_L(x, y, H);
            return H;
        };
    };
}
inline auto prob_setter_hess_L() {
    return [](alpaqa::Problem &p,
              std::function<alpaqa::mat(alpaqa::crvec, alpaqa::crvec)> fun) -> void {
        p.hess_L = [n{p.n}, fun{std::move(fun)}](alpaqa::crvec x, alpaqa::crvec y,
                                                 alpaqa::rmat H) {
            auto &&res = fun(x, y);
            if (res.rows() != n)
                throw std::out_of_range(
                    "Number of rows of result of hess_L not consistent "
                    "with problem dimension n");
            if (res.cols() != n)
                throw std::out_of_range("Number of columns of result "
                                        "of hess_L not consistent "
                                        "with problem dimension n");
            H = std::move(res);
        };
    };
}
inline auto prob_getter_hess_L_prod() {
    return [](const alpaqa::Problem &p)
               -> std::function<alpaqa::vec(alpaqa::crvec, alpaqa::crvec, alpaqa::crvec)> {
        return [n{p.n}, m{p.m}, hess_L_prod{p.hess_L_prod}](
                   alpaqa::crvec x, alpaqa::crvec y, alpaqa::crvec v) {
            if (x.size() != n)
                throw std::out_of_range("Dimension of x not consistent "
                                        "with problem dimension n");
            if (y.size() != m)
                throw std::out_of_range("Dimension of y not consistent "
                                        "with problem dimension m");
            if (v.size() != n)
                throw std::out_of_range("Dimension of v not consistent "
                                        "with problem dimension n");
            alpaqa::vec Hv(n);
            hess_L_prod(x, y, v, Hv);
            return Hv;
        };
    };
}
inline auto prob_setter_hess_L_prod() {
    return [](alpaqa::Problem &p,
              std::function<alpaqa::vec(alpaqa::crvec, alpaqa::crvec, alpaqa::crvec)> fun)
               -> void {
        p.hess_L_prod = [n{p.n}, fun{std::move(fun)}](alpaqa::crvec x, alpaqa::crvec y,
                                                      alpaqa::crvec v,
                                                      alpaqa::rvec Hv) {
            auto &&res = fun(x, y, v);
            if (res.rows() != n)
                throw std::out_of_range(
                    "Dimension of result of hess_L_prod not consistent "
                    "with problem dimension n");
            Hv = std::move(res);
        };
    };
}