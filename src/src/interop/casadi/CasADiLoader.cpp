#include <casadi/core/external.hpp>
#include <panoc-alm/interop/casadi/CasADiFunctionWrapper.hpp>
#include <panoc-alm/interop/casadi/CasADiLoader.hpp>

std::function<pa::Problem::f_sig> load_CasADi_objective(const char *so_name,
                                                        const char *fun_name) {
    return CasADiFun_1iso(casadi::external(fun_name, so_name));
}
std::function<pa::Problem::grad_f_sig>
load_CasADi_gradient_objective(const char *so_name, const char *fun_name) {
    return CasADiFun_1i1o(casadi::external(fun_name, so_name));
}
std::function<pa::Problem::g_sig>
load_CasADi_constraints(const char *so_name, const char *fun_name) {
    return CasADiFun_1i1o(casadi::external(fun_name, so_name));
}
std::function<pa::Problem::grad_g_sig>
load_CasADi_gradient_constraints(const char *so_name, const char *fun_name) {
    return CasADiFun_2i1o(casadi::external(fun_name, so_name));
}