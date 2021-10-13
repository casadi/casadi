#include <casadi/core/external.hpp>
#include <panoc-alm/interop/casadi/CasADiFunctionWrapper.hpp>
#include <panoc-alm/interop/casadi/CasADiLoader.hpp>

namespace pa {

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

pa::Problem load_CasADi_problem(const std::string &so_name, unsigned n,
                                unsigned m, bool second_order) {
    auto prob = pa::Problem(n, m);
    pa::vec w = pa::vec::Zero(m);
    auto load = [&](const std::string &name) {
        return casadi::external(name, so_name);
    };
    prob.f = [csf{CasADiFun_1Vi1So(load("f"))}] //
        (pa::crvec x) {                         //
            return csf(x);
        };
    prob.grad_f = [csf{CasADiFun_1Vi1Vo(load("grad_f"))}] //
        (pa::crvec x, pa::rvec gr) {                      //
            csf(x, gr);
        };
    prob.g = [csf{CasADiFun_1Vi1Vo(load("g"))}] //
        (pa::crvec x, pa::rvec g) {             //
            csf(x, g);
        };
    prob.grad_g_prod = [csf{CasADiFun_2Vi1Vo(load("grad_g"))}] //
        (pa::crvec x, pa::crvec y, pa::rvec g) {               //
            if (y.size() == 0)
                g.setZero();
            else
                csf(x, y, g);
        };
    if (second_order) {
        prob.grad_gi = [csf{CasADiFun_2Vi1Vo(load("grad_g"))}, w] //
            (pa::crvec x, unsigned i, pa::rvec g) mutable {       //
                if (w.size() == 0) {
                    g.setZero();
                } else {
                    w(i) = 1;
                    csf(x, w, g);
                    w(i) = 0;
                }
            };
        prob.hess_L = [csf{CasADiFun_2Vi1Mo(load("hess_L"))}] //
            (pa::crvec x, pa::crvec y, pa::rvec g) {          //
                csf(x, y, g);
            };
        prob.hess_L_prod = [csf{CasADiFun_3Vi1Vo(load("hess_L_prod"))}] //
            (pa::crvec x, pa::crvec y, pa::crvec v, pa::rvec g) {       //
                csf(x, y, v, g);
            };
    }
    return prob;
}

ProblemWithParam load_CasADi_problem_with_param(const std::string &so_name,
                                                unsigned n, unsigned m,
                                                bool second_order) {
    auto prob        = ProblemWithParam(n, m);
    pa::vec w        = pa::vec::Zero(m);
    const auto param = prob.get_param_ptr();
    auto load        = [&](const std::string &name) {
        return casadi::external(name, so_name);
    };
    prob.f = [csf{CasADiFun_2Vi1So(load("f"))}, p{param}] //
        (pa::crvec x) {                                   //
            return csf(x, *p);
        };
    prob.grad_f = [csf{CasADiFun_2Vi1Vo(load("grad_f"))}, p{param}] //
        (pa::crvec x, pa::rvec gr) {                                //
            csf(x, *p, gr);
        };
    prob.g = [csf{CasADiFun_2Vi1Vo(load("g"))}, p{param}] //
        (pa::crvec x, pa::rvec g) {                       //
            csf(x, *p, g);
        };
    prob.grad_g_prod = [csf{CasADiFun_3Vi1Vo(load("grad_g"))}, p{param}] //
        (pa::crvec x, pa::crvec y, pa::rvec g) {                         //
            if (y.size() == 0)
                g.setZero();
            else
                csf(x, *p, y, g);
        };
    if (second_order) {
        prob.grad_gi = [csf{CasADiFun_3Vi1Vo(load("grad_g"))}, p{param}, w] //
            (pa::crvec x, unsigned i, pa::rvec g) mutable {                 //
                if (w.size() == 0) {
                    g.setZero();
                } else {
                    w(i) = 1;
                    csf(x, *p, w, g);
                    w(i) = 0;
                }
            };
        prob.hess_L = [csf{CasADiFun_3Vi1Mo(load("hess_L"))}, p{param}] //
            (pa::crvec x, pa::crvec y, pa::rvec g) {                    //
                csf(x, *p, y, g);
            };
        prob.hess_L_prod =
            [csf{CasADiFun_4Vi1Vo(load("hess_L_prod"))}, p{param}] //
            (pa::crvec x, pa::crvec y, pa::crvec v, pa::rvec g) {  //
                csf(x, *p, y, v, g);
            };
    }
    return prob;
}

} // namespace pa