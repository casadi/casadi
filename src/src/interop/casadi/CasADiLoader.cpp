#include <casadi/core/external.hpp>
#include <panoc-alm/interop/casadi/CasADiFunctionWrapper.hpp>
#include <panoc-alm/interop/casadi/CasADiLoader.hpp>

std::function<pa::Problem::f_sig> load_CasADi_objective(const char *so_name,
                                                        const char *fun_name) {
    return CasADiFun_1Vi1So(casadi::external(fun_name, so_name));
}
std::function<pa::Problem::grad_f_sig>
load_CasADi_gradient_objective(const char *so_name, const char *fun_name) {
    return CasADiFun_1Vi1Vo(casadi::external(fun_name, so_name));
}
std::function<pa::Problem::g_sig>
load_CasADi_constraints(const char *so_name, const char *fun_name) {
    return CasADiFun_1Vi1Vo(casadi::external(fun_name, so_name));
}
std::function<pa::Problem::grad_g_prod_sig>
load_CasADi_gradient_constraints_prod(const char *so_name,
                                      const char *fun_name) {
    return CasADiFun_2Vi1Vo(casadi::external(fun_name, so_name));
}
std::function<pa::Problem::hess_L_sig>
load_CasADi_hessian_lagrangian(const char *so_name, const char *fun_name) {
    return [csf{CasADiFun_2Vi1Mo(casadi::external(fun_name, so_name))}](
               pa::crvec x, pa::crvec y, pa::rmat H) {
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
load_CasADi_hessian_lagrangian_prod(const char *so_name, const char *fun_name) {
    return CasADiFun_3Vi1Vo(casadi::external(fun_name, so_name));
}

pa::Problem load_CasADi_problem(const char *so_name, unsigned n, unsigned m) {
    auto prob = pa::Problem(n, m);
    pa::vec w = pa::vec::Zero(m);
    auto load = [&](const char *name) {
        return casadi::external(name, so_name);
    };
    prob.f           = [f{CasADiFun_1Vi1So(load("f"))} //
    ](pa::crvec x) { return f(x); };
    prob.grad_f      = [f{CasADiFun_1Vi1Vo(load("grad_f"))} //
    ](pa::crvec x, pa::rvec gr) { f(x, gr); };
    prob.g           = [f{CasADiFun_1Vi1Vo(load("g"))} //
    ](pa::crvec x, pa::rvec g) { f(x, g); };
    prob.grad_g_prod = [f{CasADiFun_2Vi1Vo(load("grad_g"))} //
    ](pa::crvec x, pa::crvec y, pa::rvec g) { f(x, y, g); };
    prob.grad_gi     = [f{CasADiFun_2Vi1Vo(load("grad_g"))}, w //
    ](pa::crvec x, unsigned i, pa::rvec g) mutable {
        w(i) = 1;
        f(x, w, g);
        w(i) = 0;
    };
    prob.hess_L      = [f{CasADiFun_2Vi1Mo(load("hess_L"))} //
    ](pa::crvec x, pa::crvec y, pa::rvec g) { f(x, y, g); };
    prob.hess_L_prod = [f{CasADiFun_3Vi1Vo(load("hess_L_prod"))} //
    ](pa::crvec x, pa::crvec y, pa::crvec v, pa::rvec g) { f(x, y, v, g); };
    return prob;
}

ProblemWithParam load_CasADi_problem_with_param(const char *so_name, unsigned n,
                                                unsigned m) {
    auto prob        = ProblemWithParam(n, m);
    pa::vec w        = pa::vec::Zero(m);
    const auto param = prob.get_param_ptr();
    auto load        = [&](const char *name) {
        return casadi::external(name, so_name);
    };
    prob.f           = [f{CasADiFun_2Vi1So(load("f"))}, p{param} //
    ](pa::crvec x) { return f(x, *p); };
    prob.grad_f      = [f{CasADiFun_2Vi1Vo(load("grad_f"))}, p{param} //
    ](pa::crvec x, pa::rvec gr) { f(x, *p, gr); };
    prob.g           = [f{CasADiFun_2Vi1Vo(load("g"))}, p{param} //
    ](pa::crvec x, pa::rvec g) { f(x, *p, g); };
    prob.grad_g_prod = [f{CasADiFun_3Vi1Vo(load("grad_g"))}, p{param} //
    ](pa::crvec x, pa::crvec y, pa::rvec g) { f(x, *p, y, g); };
    prob.grad_gi     = [f{CasADiFun_3Vi1Vo(load("grad_g"))}, p{param}, w //
    ](pa::crvec x, unsigned i, pa::rvec g) mutable {
        w(i) = 1;
        f(x, *p, w, g);
        w(i) = 0;
    };
    prob.hess_L      = [f{CasADiFun_3Vi1Mo(load("hess_L"))}, p{param} //
    ](pa::crvec x, pa::crvec y, pa::rvec g) { f(x, *p, y, g); };
    prob.hess_L_prod = [f{CasADiFun_4Vi1Vo(load("hess_L_prod"))}, p{param} //
    ](pa::crvec x, pa::crvec y, pa::crvec v, pa::rvec g) { f(x, *p, y, v, g); };
    return prob;
}
