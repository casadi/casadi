#pragma once

#include <alpaqa/casadi-loader-export.hpp>
#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/box.hpp>
#include <alpaqa/problem/structure.hpp>
#include <alpaqa/util/copyable_unique_ptr.hpp>

namespace alpaqa {

namespace casadi_loader {
template <Config>
struct CasADiQuadraticControlFunctionsWithParam;
} // namespace casadi_loader

template <Config Conf>
class CasADiQuadraticControlProblem {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using Box = alpaqa::Box<config_t>;
    length_t N, nx, nu;
    vec x_init;
    Box U;
    Box D;
    Box D_N;
    vec Q;     ///< State cost (nx) (Diagonal)
    vec R;     ///< Input cost (nu) (Diagonal)
    vec Q_N;   ///< Terminal cost (nx) (Diagonal)
    mat x_ref; ///< (nx×(N+1))
    mat u_ref; ///< (nu×N)
    mat μ;     ///< Penalty factors (nx×N+1)
    mat y;     ///< Lagrange multipliers (nz×(N+1)) TODO
    vec param;

    CasADiQuadraticControlProblem(const std::string &filename, length_t N,
                                  length_t nx = 0, length_t nu = 0,
                                  length_t p = 0);
    ~CasADiQuadraticControlProblem();

    CasADiQuadraticControlProblem(const CasADiQuadraticControlProblem &);
    CasADiQuadraticControlProblem &
    operator=(const CasADiQuadraticControlProblem &);
    CasADiQuadraticControlProblem(CasADiQuadraticControlProblem &&) noexcept;
    CasADiQuadraticControlProblem &
    operator=(CasADiQuadraticControlProblem &&) noexcept;

    void get_U(Box &U) const { U = this->U; }
    void get_x_init(rvec x_init) const { x_init = this->x_init; }
    void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const;
    void eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const;
    void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec v,
                          rvec grad_fxu_p) const;
    [[nodiscard]] real_t eval_l(index_t timestep, crvec h) const;
    [[nodiscard]] real_t eval_l_N(crvec h) const;
    void eval_grad_l(index_t timestep, crvec h, rvec grad_lh) const;
    void eval_grad_l_N(crvec h, rvec grad_lh) const;
    void eval_hess_l(index_t timestep, crvec h, rmat hess_lh) const;
    void eval_hess_l_N(crvec h, rmat hess_lh) const;
    [[nodiscard]] CostStructure get_l_structure() const {
        return CostStructure::DiagonalHessian;
    }
    void check() const {} // TODO

    [[nodiscard]] length_t get_N() const { return N; }
    [[nodiscard]] length_t get_nx() const { return nx; }
    [[nodiscard]] length_t get_nu() const { return nu; }
    [[nodiscard]] length_t get_nh() const { return nx + nu; }

  private:
    using Functions =
        casadi_loader::CasADiQuadraticControlFunctionsWithParam<Conf>;
    util::copyable_unique_ptr<Functions> impl;
};

CASADI_LOADER_EXPORT_EXTERN_TEMPLATE(class, CasADiQuadraticControlProblem,
                                     EigenConfigd);
CASADI_LOADER_EXPORT_EXTERN_TEMPLATE(class, CasADiQuadraticControlProblem,
                                     DefaultConfig);

} // namespace alpaqa