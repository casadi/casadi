#pragma once

#include <alpaqa/casadi-ocp-loader-export.hpp>
#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/box.hpp>
#include <alpaqa/util/check-dim.hpp>
#include <alpaqa/util/copyable_unique_ptr.hpp>
#include <filesystem>

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
    length_t N, nx, nu, nh, nh_N, nc, nc_N;
    vec x_init;
    vec param;
    Box U, D, D_N;
    mutable vec work;

    /// Components of the constraint function with indices below this number are
    /// handled using a quadratic penalty method rather than using an
    /// augmented Lagrangian method. Specifically, the Lagrange multipliers for
    /// these components (which determine the shifts in ALM) are kept at zero.
    index_t penalty_alm_split = 0;
    /// Same as @ref penalty_alm_split, but for the terminal constraint.
    index_t penalty_alm_split_N = 0;

    CasADiControlProblem(const std::string &so_name, length_t N);
    ~CasADiControlProblem();

    CasADiControlProblem(const CasADiControlProblem &);
    CasADiControlProblem &operator=(const CasADiControlProblem &);
    CasADiControlProblem(CasADiControlProblem &&) noexcept;
    CasADiControlProblem &operator=(CasADiControlProblem &&) noexcept;

    /// Load the numerical problem data (bounds and parameters) from a CSV file.
    /// The file should contain 8 rows, with the following contents:
    ///   1. @ref U lower bound [nu]
    ///   2. @ref U upper bound [nu]
    ///   3. @ref D lower bound [nc]
    ///   4. @ref D upper bound [nc]
    ///   5. @ref D_N lower bound [nc_N]
    ///   6. @ref D_N upper bound [nc_N]
    ///   7. @ref x_init [nx]
    ///   8. @ref param [p]
    ///
    /// Line endings are encoded using a single line feed (`\n`), and the column
    /// separator can be specified using the @p sep argument.
    void load_numerical_data(const std::filesystem::path &filepath,
                             char sep = ',');

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
    void eval_grad_constr_prod(index_t timestep, crvec x, crvec p,
                               rvec grad_cx_p) const;
    void eval_add_gn_hess_constr(index_t timestep, crvec x, crvec M,
                                 rmat out) const;
    void eval_constr_N(crvec x, rvec c) const;
    void eval_grad_constr_prod_N(crvec x, crvec p, rvec grad_cx_p) const;
    void eval_add_gn_hess_constr_N(crvec x, crvec M, rmat out) const;

    void check() const {
        util::check_dim_msg<config_t>(U.lowerbound, nu,
                                      "Length of problem.U.lowerbound does not "
                                      "match problem size problem.nu");
        util::check_dim_msg<config_t>(U.upperbound, nu,
                                      "Length of problem.U.upperbound does not "
                                      "match problem size problem.nu");
        util::check_dim_msg<config_t>(D.lowerbound, nc,
                                      "Length of problem.D.lowerbound does not "
                                      "match problem size problem.nc");
        util::check_dim_msg<config_t>(D.upperbound, nc,
                                      "Length of problem.D.upperbound does not "
                                      "match problem size problem.nc");
        util::check_dim_msg<config_t>(D_N.lowerbound, nc_N,
                                      "Length of problem.D_N.lowerbound does "
                                      "not match problem size problem.nc_N");
        util::check_dim_msg<config_t>(D_N.upperbound, nc_N,
                                      "Length of problem.D_N.upperbound does "
                                      "not match problem size problem.nc_N");
        if (penalty_alm_split < 0 || penalty_alm_split > nc)
            throw std::invalid_argument("Invalid penalty_alm_split");
        if (penalty_alm_split_N < 0 || penalty_alm_split > nc_N)
            throw std::invalid_argument("Invalid penalty_alm_split_N");
    }

    [[nodiscard]] length_t get_N() const { return N; }
    [[nodiscard]] length_t get_nx() const { return nx; }
    [[nodiscard]] length_t get_nu() const { return nu; }
    [[nodiscard]] length_t get_nh() const { return nh; }
    [[nodiscard]] length_t get_nh_N() const { return nh_N; }
    [[nodiscard]] length_t get_nc() const { return nc; }
    [[nodiscard]] length_t get_nc_N() const { return nc_N; }

    /// @see @ref TypeErasedControlProblem::eval_proj_diff_g
    void eval_proj_diff_g(crvec z, rvec e) const {
        for (index_t t = 0; t < N; ++t)
            e.segment(t * nc, nc) =
                alpaqa::projecting_difference(z.segment(t * nc, nc), D);
        e.segment(N * nc, nc_N) =
            alpaqa::projecting_difference(z.segment(N * nc, nc_N), D_N);
    }
    /// @see @ref TypeErasedControlProblem::eval_proj_multipliers
    void eval_proj_multipliers(rvec y, real_t M) const {
        // If there's no lower bound, the multipliers can only be positive
        auto max_lb = [M](real_t y, real_t z_lb) {
            real_t y_lb = z_lb == -alpaqa::inf<config_t> ? 0 : -M;
            return std::max(y, y_lb);
        };
        // If there's no upper bound, the multipliers can only be negative
        auto min_ub = [M](real_t y, real_t z_ub) {
            real_t y_ub = z_ub == alpaqa::inf<config_t> ? 0 : M;
            return std::min(y, y_ub);
        };
        for (index_t t = 0; t < N; ++t) {
            auto num_alm    = nc - penalty_alm_split;
            auto &&yt       = y.segment(t * nc, nc);
            auto &&y_qpm    = yt.topRows(penalty_alm_split);
            auto &&y_alm    = yt.bottomRows(num_alm);
            auto &&z_alm_lb = D.lowerbound.bottomRows(num_alm);
            auto &&z_alm_ub = D.upperbound.bottomRows(num_alm);
            y_qpm.setZero();
            y_alm =
                y_alm.binaryExpr(z_alm_lb, max_lb).binaryExpr(z_alm_ub, min_ub);
        }
        {
            auto &&yt       = y.segment(N * nc, nc_N);
            auto num_alm    = nc_N - penalty_alm_split_N;
            auto &&y_qpm    = yt.topRows(penalty_alm_split_N);
            auto &&y_alm    = yt.bottomRows(num_alm);
            auto &&z_alm_lb = D.lowerbound.bottomRows(num_alm);
            auto &&z_alm_ub = D.upperbound.bottomRows(num_alm);
            y_qpm.setZero();
            y_alm =
                y_alm.binaryExpr(z_alm_lb, max_lb).binaryExpr(z_alm_ub, min_ub);
        }
    }

  private:
    using Functions = casadi_loader::CasADiControlFunctionsWithParam<Conf>;
    util::copyable_unique_ptr<Functions> impl;
};

CASADI_OCP_LOADER_EXPORT_EXTERN_TEMPLATE(class, CasADiControlProblem,
                                         EigenConfigd);
CASADI_OCP_LOADER_EXPORT_EXTERN_TEMPLATE(class, CasADiControlProblem,
                                         DefaultConfig);

} // namespace alpaqa