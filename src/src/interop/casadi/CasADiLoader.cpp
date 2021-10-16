#include <casadi/core/external.hpp>
#include <optional>
#include <panoc-alm/interop/casadi/CasADiFunctionWrapper.hpp>
#include <panoc-alm/interop/casadi/CasADiLoader.hpp>
#include <stdexcept>

namespace pa {

#if 0
std::function<pa::Problem::f_sig>
load_CasADi_objective(const std::string &so_name, const std::string &fun_name) {
    return CasADiFun_1Vi1So(casadi::external(fun_name, so_name));
}
std::function<pa::Problem::grad_f_sig>
load_CasADi_gradient_objective(const std::string &so_name,
                               const std::string &fun_name) {
    return CasADiFun_1Vi1Vo(casadi::external(fun_name, so_name));
}
std::function<pa::Problem::g_sig>
load_CasADi_constraints(const std::string &so_name,
                        const std::string &fun_name) {
    return CasADiFun_1Vi1Vo(casadi::external(fun_name, so_name));
}
std::function<pa::Problem::grad_g_prod_sig>
load_CasADi_gradient_constraints_prod(const std::string &so_name,
                                      const std::string &fun_name) {
    return [csf{CasADiFun_2Vi1Vo(casadi::external(fun_name, so_name))}] //
        (pa::crvec x, pa::crvec y, pa::rvec gradprod) {                 //
            if (y.size() == 0)
                gradprod.setZero();
            else
                csf(x, y, gradprod);
        };
}
std::function<pa::Problem::hess_L_sig>
load_CasADi_hessian_lagrangian(const std::string &so_name,
                               const std::string &fun_name) {
    return [csf{CasADiFun_2Vi1Mo(casadi::external(fun_name, so_name))}] //
        (pa::crvec x, pa::crvec y, pa::rmat H) {                        //
            // Fix the stride if the matrix is larger than n
            if (x.rows() != H.rows()) { // TODO: this is probably unnecessary
                for (auto c = x.rows(); c-- > 1;)
                    for (auto r = x.rows(); r-- > 0;)
                        std::swap(H(r, c), H.data()[r + x.rows() * c]);
            }
            csf(x, y, H);
            // Fix the stride if the matrix is larger than n
            if (x.rows() != H.rows()) {
                for (auto c = x.rows(); c-- > 1;)
                    for (auto r = x.rows(); r-- > 0;)
                        std::swap(H(r, c), H.data()[r + x.rows() * c]);
            }
        };
}
std::function<pa::Problem::hess_L_prod_sig>
load_CasADi_hessian_lagrangian_prod(const std::string &so_name,
                                    const std::string &fun_name) {
    return CasADiFun_3Vi1Vo(casadi::external(fun_name, so_name));
}
#endif

pa::Problem load_CasADi_problem(const std::string &so_name, unsigned n,
                                unsigned m, bool second_order) {

    auto load = [&](const std::string &name) {
        return casadi::external(name, so_name);
    };
    auto wrap_load = [&so_name](const char *name, auto f) {
        try {
            f();
        } catch (const std::invalid_argument &e) {
            throw std::invalid_argument("Unable to load function '" + so_name +
                                        ":" + name + "': " + e.what());
        }
    };

    std::optional<CasADiFunctionEvaluator<1, 1>> g;
    // If not all dimensions are specified, load the function "g" to determine
    // the missing dimensions.
    if (n == 0 || m == 0) {
        wrap_load("g", [&] {
            g = load("g");
            if (g->fun.size2_in(0) != 1)
                throw std::invalid_argument(
                    "First input argument should be a column vector.");
            if (g->fun.size2_out(0) != 1)
                throw std::invalid_argument(
                    "First output argument should be a column vector.");
            if (n == 0)
                n = g->fun.size1_in(0);
            if (m == 0)
                m = g->fun.size1_out(0);
            g->validate_dimensions({{n, 1}}, {{m, 1}});
        });
    }
    // Otherwise, load the function "g" and compare its dimensions to the
    // dimensions specified by the user.
    else {
        wrap_load("g", [&] { g = {load("g"), {{n, 1}}, {{m, 1}}}; });
    }

    auto prob = pa::Problem(n, m);
    pa::vec w = pa::vec::Zero(m);

    wrap_load("f", [&] {
        prob.f = [csf{CasADiFun_1Vi1So(load("f"), n)}] //
            (pa::crvec x) {                            //
                return csf(x);
            };
    });
    wrap_load("grad_f", [&] {
        prob.grad_f = [csf{CasADiFun_1Vi1Vo(load("grad_f"), n, n)}] //
            (pa::crvec x, pa::rvec gr) {                            //
                csf(x, gr);
            };
    });
    wrap_load("g", [&] {
        prob.g = [csf{CasADiFun_1Vi1Vo(std::move(*g))}] //
            (pa::crvec x, pa::rvec g) {                 //
                csf(x, g);
            };
    });
    wrap_load("grad_g", [&] {
        prob.grad_g_prod = [csf{CasADiFun_2Vi1Vo(load("grad_g"), {n, m}, n)}] //
            (pa::crvec x, pa::crvec y, pa::rvec g) {                          //
                if (y.size() == 0)
                    g.setZero();
                else
                    csf(x, y, g);
            };
    });
    if (second_order) {
        wrap_load("grad_gi", [&] {
            prob.grad_gi =
                [csf{CasADiFun_2Vi1Vo(load("grad_g"), {n, m}, n)}, w] //
                (pa::crvec x, unsigned i, pa::rvec g) mutable {       //
                    if (w.size() == 0) {
                        g.setZero();
                    } else {
                        w(i) = 1;
                        csf(x, w, g);
                        w(i) = 0;
                    }
                };
        });
        wrap_load("hess_L", [&] {
            prob.hess_L =
                [csf{CasADiFun_2Vi1Mo(load("hess_L"), {n, m}, {n, n})}] //
                (pa::crvec x, pa::crvec y, pa::rvec g) {                //
                    csf(x, y, g);
                };
        });
        wrap_load("hess_L_prod", [&] {
            prob.hess_L_prod =
                [csf{CasADiFun_3Vi1Vo(load("hess_L_prod"), {n, m, n}, n)}] //
                (pa::crvec x, pa::crvec y, pa::crvec v, pa::rvec g) {      //
                    csf(x, y, v, g);
                };
        });
    }
    return prob;
}

ProblemWithParam load_CasADi_problem_with_param(const std::string &so_name,
                                                unsigned n, unsigned m,
                                                unsigned p, bool second_order) {

    auto load = [&](const std::string &name) {
        return casadi::external(name, so_name);
    };
    auto wrap_load = [&so_name](const char *name, auto f) {
        try {
            f();
        } catch (const std::invalid_argument &e) {
            throw std::invalid_argument("Unable to load function '" + so_name +
                                        ":" + name + "': " + e.what());
        }
    };

    std::optional<CasADiFunctionEvaluator<2, 1>> g;
    // If not all dimensions are specified, load the function "g" to determine
    // the missing dimensions.
    if (n == 0 || m == 0 || p == 0) {
        wrap_load("g", [&] {
            g = load("g");
            if (g->fun.size2_in(0) != 1)
                throw std::invalid_argument(
                    "First input argument should be a column vector.");
            if (g->fun.size2_in(1) != 1)
                throw std::invalid_argument(
                    "Second input argument should be a column vector.");
            if (g->fun.size2_out(0) != 1)
                throw std::invalid_argument(
                    "First output argument should be a column vector.");
            if (n == 0)
                n = g->fun.size1_in(0);
            if (m == 0)
                m = g->fun.size1_out(0);
            if (p == 0)
                p = g->fun.size1_in(1);
            g->validate_dimensions({{n, 1}, {p, 1}}, {{m, 1}});
        });
    }
    // Otherwise, load the function "g" and compare its dimensions to the
    // dimensions specified by the user.
    else {
        wrap_load("g", [&] { g = {load("g"), {{n, 1}, {p, 1}}, {{m, 1}}}; });
    }

    auto prob        = ProblemWithParam(n, m, p);
    pa::vec w        = pa::vec::Zero(m);
    const auto param = prob.get_param_ptr();

    wrap_load("f", [&] {
        prob.f = [csf{CasADiFun_2Vi1So(load("f"), {n, p})}, p{param}] //
            (pa::crvec x) {                                           //
                return csf(x, *p);
            };
    });
    wrap_load("grad_f", [&] {
        prob.grad_f =
            [csf{CasADiFun_2Vi1Vo(load("grad_f"), {n, p}, n)}, p{param}] //
            (pa::crvec x, pa::rvec gr) {                                 //
                csf(x, *p, gr);
            };
    });
    wrap_load("g", [&] {
        prob.g = [csf{CasADiFun_2Vi1Vo(std::move(*g))}, p{param}] //
            (pa::crvec x, pa::rvec g) {                           //
                csf(x, *p, g);
            };
    });
    wrap_load("grad_g", [&] {
        prob.grad_g_prod =
            [csf{CasADiFun_3Vi1Vo(load("grad_g"), {n, p, m}, n)}, p{param}] //
            (pa::crvec x, pa::crvec y, pa::rvec g) {                        //
                if (y.size() == 0)
                    g.setZero();
                else
                    csf(x, *p, y, g);
            };
    });
    if (second_order) {
        wrap_load("grad_g", [&] {
            prob.grad_gi = [csf{CasADiFun_3Vi1Vo(load("grad_g"), {n, p, m}, n)},
                            p{param}, w]                        //
                (pa::crvec x, unsigned i, pa::rvec g) mutable { //
                    if (w.size() == 0) {
                        g.setZero();
                    } else {
                        w(i) = 1;
                        csf(x, *p, w, g);
                        w(i) = 0;
                    }
                };
        });
        wrap_load("hess_L", [&] {
            prob.hess_L =
                [csf{CasADiFun_3Vi1Mo(load("hess_L"), {n, p, m}, {n, n})},
                 p{param}]                               //
                (pa::crvec x, pa::crvec y, pa::rvec g) { //
                    csf(x, *p, y, g);
                };
        });
        wrap_load("hess_L_prod", [&] {
            prob.hess_L_prod =
                [csf{CasADiFun_4Vi1Vo(load("hess_L_prod"), {n, p, m, n}, n)},
                 p{param}]                                            //
                (pa::crvec x, pa::crvec y, pa::crvec v, pa::rvec g) { //
                    csf(x, *p, y, v, g);
                };
        });
    }
    return prob;
}

} // namespace pa