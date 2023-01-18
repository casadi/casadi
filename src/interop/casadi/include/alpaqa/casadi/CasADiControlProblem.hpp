#pragma once

#include <alpaqa/casadi-ocp-loader-export.hpp>
#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/box.hpp>
#include <alpaqa/util/copyable_unique_ptr.hpp>

namespace alpaqa {

namespace casadi_loader {
template <Config>
struct CasADiControlFunctionsWithParam;
} // namespace casadi_loader

template <Config Conf>
class CasADiControlProblem {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using Box = alpaqa::Box<config_t>;
    length_t N, nx, nu, nh, nh_N, nc;
    vec x_init;
    vec param;
    Box U, D, D_N;
    mutable vec work;

    CasADiControlProblem(const std::string &so_name, length_t N);
    ~CasADiControlProblem();

    CasADiControlProblem(const CasADiControlProblem &);
    CasADiControlProblem &operator=(const CasADiControlProblem &);
    CasADiControlProblem(CasADiControlProblem &&) noexcept;
    CasADiControlProblem &operator=(CasADiControlProblem &&) noexcept;

    void get_U(Box &U) const { U = this->U; }
    void get_D(Box &D) const { D = this->D; }
    void get_D_N(Box &D_N) const { D_N = this->D_N; }
    void get_x_init(rvec x_init) const { x_init = this->x_init; }

    void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const;
    void eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const;
    void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p,
                          rvec grad_fxu_p) const;
    void eval_h(index_t timestep, crvec x, crvec u, rvec h) const;
    void eval_h_N(crvec x, rvec h) const;
    [[nodiscard]] real_t eval_l(index_t timestep, crvec h) const;
    [[nodiscard]] real_t eval_l_N(crvec h) const;
    void eval_qr(index_t timestep, crvec xu, crvec h, rvec qr) const;
    void eval_q_N(crvec x, crvec h, rvec q) const;
    void eval_add_Q(index_t timestep, crvec xu, crvec h, rmat Q) const;
    void eval_add_Q_N(crvec x, crvec h, rmat Q) const;
    void eval_add_R_masked(index_t timestep, crvec xu, crvec h, crindexvec mask,
                           rmat R, rvec work) const;
    void eval_add_S_masked(index_t timestep, crvec xu, crvec h, crindexvec mask,
                           rmat S, rvec work) const;
    void eval_add_R_prod_masked(index_t timestep, crvec xu, crvec h,
                                crindexvec mask_J, crindexvec mask_K, crvec v,
                                rvec out, rvec work) const;
    void eval_add_S_prod_masked(index_t timestep, crvec xu, crvec h,
                                crindexvec mask_K, crvec v, rvec out,
                                rvec work) const;
    [[nodiscard]] length_t get_R_work_size() const;
    [[nodiscard]] length_t get_S_work_size() const;
    void eval_constr(index_t timestep, crvec x, rvec c) const;
    void eval_jac_constr(index_t timestep, crvec x, rmat J_c) const;
    void eval_grad_constr_prod(index_t timestep, crvec x, crvec p,
                               rvec grad_cx_p) const;
    void eval_add_gn_hess_constr(index_t timestep, crvec x, crvec M,
                                 rmat out) const;

    void check() const {} // TODO

    [[nodiscard]] length_t get_N() const { return N; }
    [[nodiscard]] length_t get_nx() const { return nx; }
    [[nodiscard]] length_t get_nu() const { return nu; }
    [[nodiscard]] length_t get_nh() const { return nh; }
    [[nodiscard]] length_t get_nh_N() const { return nh_N; }
    [[nodiscard]] length_t get_nc() const { return nc; }
    [[nodiscard]] length_t get_nc_N() const { return nc; } // TODO

  private:
    using Functions = casadi_loader::CasADiControlFunctionsWithParam<Conf>;
    util::copyable_unique_ptr<Functions> impl;
};

CASADI_OCP_LOADER_EXPORT_EXTERN_TEMPLATE(class, CasADiControlProblem,
                                         EigenConfigd);
CASADI_OCP_LOADER_EXPORT_EXTERN_TEMPLATE(class, CasADiControlProblem,
                                         DefaultConfig);

} // namespace alpaqa