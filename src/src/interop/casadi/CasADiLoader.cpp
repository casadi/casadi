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