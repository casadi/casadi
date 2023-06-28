#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/box.hpp>
#include <alpaqa/problem/ocproblem-counters.hpp>
#include <alpaqa/problem/problem-counters.hpp>
#include <alpaqa/util/not-implemented.hpp>
#include <alpaqa/util/required-method.hpp>
#include <alpaqa/util/timed.hpp>
#include <alpaqa/util/type-erasure.hpp>
#include <array>
#include <concepts>
#include <stdexcept>
#include <type_traits>

#if !ALPAQA_WITH_OCP
#error "OCP support disabled"
#endif

#ifndef NDEBUG
#include <iostream>
#endif

namespace alpaqa {

template <Config Conf>
struct OCPDim {
    USING_ALPAQA_CONFIG(Conf);
    length_t N, nx, nu, nh, nc;
};

template <Config Conf>
struct ControlProblemVTable : util::BasicVTable {
    USING_ALPAQA_CONFIG(Conf);
    using Box = alpaqa::Box<config_t>;

    template <class F>
    using optional_function_t = util::BasicVTable::optional_function_t<F, ControlProblemVTable>;
    template <class F>
    using optional_const_function_t =
        util::BasicVTable::optional_const_function_t<F, ControlProblemVTable>;

    // clang-format off
    required_const_function_t<void(crvec z, rvec e)>
        eval_proj_diff_g;
    required_const_function_t<void(rvec y, real_t M)>
        eval_proj_multipliers;
    required_const_function_t<void(Box &U)>
        get_U;
    optional_const_function_t<void(Box &D)>
        get_D = nullptr;
    optional_const_function_t<void(Box &D)>
        get_D_N = &default_get_D_N;
    required_const_function_t<void(rvec x_init)>
        get_x_init;
    required_const_function_t<void(index_t timestep, crvec x, crvec u, rvec fxu)>
        eval_f;
    required_const_function_t<void(index_t timestep, crvec x, crvec u, rmat J_fxu)>
        eval_jac_f;
    required_const_function_t<void(index_t timestep, crvec x, crvec u, crvec p, rvec grad_fxu_p)>
        eval_grad_f_prod;
    optional_const_function_t<void(index_t timestep, crvec x, crvec u, rvec h)>
        eval_h = nullptr;
    optional_const_function_t<void(crvec x, rvec h)>
        eval_h_N = nullptr;
    required_const_function_t<real_t(index_t timestep, crvec h)>
        eval_l;
    required_const_function_t<real_t(crvec h)>
        eval_l_N;
    required_const_function_t<void(index_t timestep, crvec xu, crvec h, rvec qr)>
        eval_qr;
    required_const_function_t<void(crvec x, crvec h, rvec q)>
        eval_q_N;
    required_const_function_t<void(index_t timestep, crvec xu, crvec h, rmat Q)>
        eval_add_Q;
    optional_const_function_t<void(crvec x, crvec h, rmat Q)>
        eval_add_Q_N = &default_eval_add_Q_N;
    required_const_function_t<void(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat R, rvec work)>
        eval_add_R_masked;
    required_const_function_t<void(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat S, rvec work)>
        eval_add_S_masked;
    optional_const_function_t<void(index_t timestep, crvec xu, crvec h, crindexvec mask_J, crindexvec mask_K, crvec v, rvec out, rvec work)>
        eval_add_R_prod_masked = &default_eval_add_R_prod_masked;
    optional_const_function_t<void(index_t timestep, crvec xu, crvec h, crindexvec mask_K, crvec v, rvec out, rvec work)>
        eval_add_S_prod_masked = &default_eval_add_S_prod_masked;
    optional_const_function_t<length_t()>
        get_R_work_size = &default_get_R_work_size;
    optional_const_function_t<length_t()>
        get_S_work_size = &default_get_S_work_size;
    optional_const_function_t<void(index_t timestep, crvec x, rvec c)>
        eval_constr = nullptr;
    optional_const_function_t<void(crvec x, rvec c)>
        eval_constr_N = &default_eval_constr_N;
    optional_const_function_t<void(index_t timestep, crvec x, crvec p, rvec grad_cx_p)>
        eval_grad_constr_prod = nullptr;
    optional_const_function_t<void(crvec x, crvec p, rvec grad_cx_p)>
        eval_grad_constr_prod_N = &default_eval_grad_constr_prod_N;
    optional_const_function_t<void(index_t timestep, crvec x, crvec M, rmat out)>
        eval_add_gn_hess_constr = nullptr;
    optional_const_function_t<void(crvec x, crvec M, rmat out)>
        eval_add_gn_hess_constr_N = &default_eval_add_gn_hess_constr_N;
    required_const_function_t<void()>
        check;
    // clang-format on

    length_t N, nu, nx, nh, nh_N, nc, nc_N;

    template <class P>
    ControlProblemVTable(std::in_place_t, P &p) : util::BasicVTable{std::in_place, p} {
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_proj_diff_g);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_proj_multipliers);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, get_U);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, get_D, p);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, get_D_N, p);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, get_x_init);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_f);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_jac_f);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_grad_f_prod);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_h, p);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_h_N, p);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_l);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_l_N);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_qr);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_q_N);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_add_Q);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_add_Q_N, p);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_add_R_masked);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_add_S_masked);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_add_R_prod_masked, p);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_add_S_prod_masked, p);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, get_R_work_size, p);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, get_S_work_size, p);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_constr, p);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_constr_N, p);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_grad_constr_prod, p);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_grad_constr_prod_N, p);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_add_gn_hess_constr, p);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_add_gn_hess_constr_N, p);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, check);
        N    = p.get_N();
        nu   = p.get_nu();
        nx   = p.get_nx();
        nh   = p.get_nh();
        nh_N = p.get_nh_N();
        nc   = p.get_nc();
        nc_N = p.get_nc_N();
        if (nc > 0 && get_D == nullptr)
            throw std::runtime_error("ControlProblem: missing 'get_D'");
        if (nc > 0 && eval_constr == nullptr)
            throw std::runtime_error("ControlProblem: missing 'eval_constr'");
        if (nc > 0 && eval_grad_constr_prod == nullptr)
            throw std::runtime_error("ControlProblem: missing 'eval_grad_constr_prod'");
        if (nh > 0 && eval_h == nullptr)
            throw std::runtime_error("ControlProblem: missing 'eval_h'");
        if (nh_N > 0 && eval_h_N == nullptr)
            throw std::runtime_error("ControlProblem: missing 'eval_h_N'");
    }
    ControlProblemVTable() = default;

    ALPAQA_EXPORT static void default_get_D_N(const void *self, Box &D,
                                              const ControlProblemVTable &vtable);
    ALPAQA_EXPORT static void default_eval_add_Q_N(const void *self, crvec x, crvec h, rmat Q,
                                                   const ControlProblemVTable &vtable);
    ALPAQA_EXPORT static void default_eval_add_R_prod_masked(const void *, index_t, crvec, crvec,
                                                             crindexvec, crindexvec, crvec, rvec,
                                                             rvec, const ControlProblemVTable &);
    ALPAQA_EXPORT static void default_eval_add_S_prod_masked(const void *, index_t, crvec, crvec,
                                                             crindexvec, crvec, rvec, rvec,
                                                             const ControlProblemVTable &);
    [[nodiscard]] ALPAQA_EXPORT static length_t
    default_get_R_work_size(const void *, const ControlProblemVTable &);
    [[nodiscard]] ALPAQA_EXPORT static length_t
    default_get_S_work_size(const void *, const ControlProblemVTable &);
    ALPAQA_EXPORT static void default_eval_constr_N(const void *self, crvec x, rvec c,
                                                    const ControlProblemVTable &vtable);
    ALPAQA_EXPORT static void default_eval_grad_constr_prod_N(const void *self, crvec x, crvec p,
                                                              rvec grad_cx_p,
                                                              const ControlProblemVTable &vtable);
    ALPAQA_EXPORT static void default_eval_add_gn_hess_constr_N(const void *self, crvec x, crvec M,
                                                                rmat out,
                                                                const ControlProblemVTable &vtable);
};

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ControlProblemVTable, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ControlProblemVTable, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ControlProblemVTable, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ControlProblemVTable, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ControlProblemVTable, EigenConfigq);
#endif

/**
 * Nonlinear optimal control problem with finite horizon @f$ N @f$.
 * @f[
 * \newcommand\U{U}
 * \newcommand\D{D}
 * \newcommand\nnu{{n_u}}
 * \newcommand\nnx{{n_x}}
 * \newcommand\nny{{n_y}}
 * \newcommand\xinit{x_\text{init}}
 * \begin{equation}\label{eq:OCP} \tag{OCP}\hspace{-0.8em}
 *     \begin{aligned}
 *         &\minimize_{u,x} && \sum_{k=0}^{N-1} \ell_k\big(h_k(x^k, u^k)\big) + \ell_N\big(h_N(x^N)\big)\hspace{-0.8em} \\
 *         &\subjto && u^k \in \U \\
 *         &&& c_k(x^k) \in \D \\
 *         &&& c_N(x^N) \in \D_N \\
 *         &&& x^0 = \xinit \\
 *         &&& x^{k+1} = f(x^k, u^k) \quad\quad (0 \le k \lt N)
 *     \end{aligned}
 * \end{equation}
 * @f]
 * 
 * The function @f$ f : \R^\nnx \times \R^\nnu \to \R^\nnx @f$ models the 
 * discrete-time, nonlinear dynamics of the system, which starts from an initial
 * state @f$ \xinit @f$. 
 * The functions @f$ h_k : \R^\nnx \times \R^\nnu \to \R^{n_h} @f$ for
 * @f$ 0 \le k \lt N @f$ and @f$ h_N : \R^\nnx \to \R^{n_h^N} @f$ can be used to
 * represent the (possibly time-varying) output mapping of the system,
 * and the convex functions @f$ \ell_k : \R^{n_h} \to \R @f$ and
 * @f$ \ell_N : \R^{n_h^N} \to \R @f$ define the stage costs and the terminal
 * cost respectively. Stage constraints and terminal constraints are represented
 * by the functions @f$ c_k : \R^{n_x} \to \R^{n_c} @f$ and
 * @f$ c_N : \R^{n_x} \to \R^{n_c^N} @f$, and the boxes @f$ D @f$ and
 * @f$ D_N @f$.
 *
 * Additional functions for computing Gauss-Newton approximations of the cost
 * Hessian are included as well:
 * @f[ \begin{aligned}
 * q^k &\defeq \tp{\jac_{h_k}^x\!(\barxuk)} \nabla \ell_k(\hhbar^k) \\
 * r^k &\defeq \tp{\jac_{h_k}^u\!(\barxuk)} \nabla \ell_k(\hhbar^k) \\
 * \Lambda_k &\defeq \partial^2 \ell_k(\hhbar^k) \\
 * Q_k &\defeq \tp{\jac_{h_k}^x\!(\barxuk)} \Lambda_k\, \jac_{h_k}^x\!(\barxuk) \\
 * S_k &\defeq \tp{\jac_{h_k}^u\!(\barxuk)} \Lambda_k\, \jac_{h_k}^x\!(\barxuk) \\
 * R_k &\defeq \tp{\jac_{h_k}^u\!(\barxuk)} \Lambda_k\, \jac_{h_k}^u\!(\barxuk). \\
 * \end{aligned} @f]
 * See @cite pas2022gaussnewton for more details.
 *
 * @ingroup grp_Problems
 */
template <Config Conf = DefaultConfig, class Allocator = std::allocator<std::byte>>
class TypeErasedControlProblem : public util::TypeErased<ControlProblemVTable<Conf>, Allocator> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using VTable         = ControlProblemVTable<config_t>;
    using allocator_type = Allocator;
    using Box            = typename VTable::Box;
    using Dim            = OCPDim<config_t>;
    using TypeErased     = util::TypeErased<VTable, allocator_type>;
    using TypeErased::TypeErased;

  protected:
    using TypeErased::call;
    using TypeErased::self;
    using TypeErased::vtable;

  public:
    template <class T, class... Args>
    static TypeErasedControlProblem make(Args &&...args) {
        return TypeErased::template make<TypeErasedControlProblem, T>(std::forward<Args>(args)...);
    }

    /// @name Problem dimensions
    /// @{

    /// Horizon length.
    [[nodiscard]] length_t get_N() const { return vtable.N; }
    /// Number of inputs.
    [[nodiscard]] length_t get_nu() const { return vtable.nu; }
    /// Number of states.
    [[nodiscard]] length_t get_nx() const { return vtable.nx; }
    /// Number of outputs.
    [[nodiscard]] length_t get_nh() const { return vtable.nh; }
    [[nodiscard]] length_t get_nh_N() const { return vtable.nh_N; }
    /// Number of constraints.
    [[nodiscard]] length_t get_nc() const { return vtable.nc; }
    [[nodiscard]] length_t get_nc_N() const { return vtable.nc_N; }
    /// All dimensions
    [[nodiscard]] Dim get_dim() const {
        return {
            .N    = vtable.N,
            .nx   = vtable.nx,
            .nu   = vtable.nu,
            .nh   = vtable.nh,
            .nh_N = vtable.nh_N,
            .nc   = vtable.nc,
            .nc_N = vtable.nc_N,
        };
    }
    /// Total number of variables.
    [[nodiscard]] length_t get_n() const { return get_N() * get_nu(); }
    /// Total number of constraints.
    [[nodiscard]] length_t get_m() const { return get_N() * get_nc() + get_nc_N(); }

    /// @}

    /// @name Projections onto constraint sets
    /// @{

    /// **[Required]**
    /// Function that evaluates the difference between the given point @f$ z @f$
    /// and its projection onto the constraint set @f$ D @f$.
    /// @param  [in] z
    ///         Slack variable, @f$ z \in \R^m @f$
    /// @param  [out] e
    ///         The difference relative to its projection,
    ///         @f$ e = z - \Pi_D(z) \in \R^m @f$
    /// @note   @p z and @p e can refer to the same vector.
    void eval_proj_diff_g(crvec z, rvec e) const;
    /// **[Required]**
    /// Function that projects the Lagrange multipliers for ALM.
    /// @param  [inout] y
    ///         Multipliers, @f$ y \leftarrow \Pi_Y(y) \in \R^m @f$
    /// @param  [in] M
    ///         The radius/size of the set @f$ Y @f$.
    ///         See @ref ALMParams::max_multiplier.
    void eval_proj_multipliers(rvec y, real_t M) const;

    /// @}

    /// @name Constraint sets
    /// @{

    /// Input box constraints @f$ U @f$.
    void get_U(Box &U) const;
    /// Stage box constraints @f$ D @f$.
    void get_D(Box &D) const;
    /// Terminal box constraints @f$ D_N @f$.
    void get_D_N(Box &D) const;

    /// @}

    /// @name Dynamics and initial state
    /// @{

    /// Initial state @f$ x_\text{init} @f$.
    void get_x_init(rvec x_init) const;
    /// Discrete-time dynamics @f$ x^{k+1} = f_k(x^k, u^k) @f$.
    void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const;
    /// Jacobian of discrete-time dynamics @f$ \jac_f(x^k, u^k) @f$.
    void eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const;
    /// Gradient-vector product of discrete-time dynamics @f$ \nabla f(x^k, u^k)\,p @f$.
    void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p, rvec grad_fxu_p) const;

    /// @}

    /// @name Output mapping
    /// @{

    /// Stage output mapping @f$ h_k(x^k, u^k) @f$.
    void eval_h(index_t timestep, crvec x, crvec u, rvec h) const;
    /// Terminal output mapping @f$ h_N(x^N) @f$.
    void eval_h_N(crvec x, rvec h) const;

    /// @}

    /// @name Stage and terminal cost
    /// @{

    /// Stage cost @f$ \ell_k(\hbar^k) @f$.
    [[nodiscard]] real_t eval_l(index_t timestep, crvec h) const;
    /// Terminal cost @f$ \ell_N(\hbar^N) @f$.
    [[nodiscard]] real_t eval_l_N(crvec h) const;

    /// @}

    /// @name Gauss-Newton approximations
    /// @{

    /// Cost gradients w.r.t. states and inputs
    /// @f$ q^k = \tp{\jac_{h_k}^x\!(\barxuk)} \nabla \ell_k(\hbar^k) @f$ and
    /// @f$ r^k = \tp{\jac_{h_k}^u\!(\barxuk)} \nabla \ell_k(\hbar^k) @f$.
    void eval_qr(index_t timestep, crvec xu, crvec h, rvec qr) const;
    /// Terminal cost gradient w.r.t. states
    /// @f$ q^N = \tp{\jac_{h_N}(\bar x^N)} \nabla \ell_k(\hbar^N) @f$.
    void eval_q_N(crvec x, crvec h, rvec q) const;
    /// Cost Hessian w.r.t. states @f$ Q_k = \tp{\jac_{h_k}^x\!(\barxuk)}
    /// \partial^2\ell_k(\hbar^k)\, \jac_{h_k}^x\!(\barxuk) @f$,
    /// added to the given matrix @p Q.
    /// @f$ Q \leftarrow Q + Q_k @f$.
    void eval_add_Q(index_t timestep, crvec xu, crvec h, rmat Q) const;
    /// Terminal cost Hessian w.r.t. states @f$ Q_N = \tp{\jac_{h_N}(\bar x^N)}
    /// \partial^2\ell_N(\hbar^N)\, \jac_{h_N}(\bar x^N) @f$,
    /// added to the given matrix @p Q.
    /// @f$ Q \leftarrow Q + Q_N @f$.
    void eval_add_Q_N(crvec x, crvec h, rmat Q) const;
    /// Cost Hessian w.r.t. inputs @f$ R_k = \tp{\jac_{h_k}^u\!(\barxuk)}
    /// \partial^2\ell_k(\hbar^k)\, \jac_{h_k}^u\!(\barxuk) @f$, keeping only
    /// rows and columns in the mask @f$ \mathcal J @f$, added to the given
    /// matrix @p R.
    /// @f$ R \leftarrow R + R_k[\mathcal J, \mathcal J] @f$.
    /// The size of @p work should be @ref get_R_work_size().
    void eval_add_R_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat R,
                           rvec work) const;
    /// Cost Hessian w.r.t. inputs and states @f$ S_k = \tp{\jac_{h_k}^u\!(\barxuk)}
    /// \partial^2\ell_k(\hbar^k)\, \jac_{h_k}^x\!(\barxuk) @f$, keeping only
    /// rows in the mask @f$ \mathcal J @f$, added to the given matrix @p S.
    /// @f$ S \leftarrow S + S_k[\mathcal J, \cdot] @f$.
    /// The size of @p work should be @ref get_S_work_size().
    void eval_add_S_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat S,
                           rvec work) const;
    /// @f$ out \leftarrow out + R[\mathcal J, \mathcal K]\,v[\mathcal K] @f$.
    /// Work should contain the contents written to it by a prior call to
    /// @ref eval_add_R_masked() in the same point.
    void eval_add_R_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_J,
                                crindexvec mask_K, crvec v, rvec out, rvec work) const;
    /// @f$ out \leftarrow out + \tp{S[\mathcal K, \cdot]}\, v[\mathcal K] @f$.
    /// Work should contain the contents written to it by a prior call to
    /// @ref eval_add_S_masked() in the same point.
    void eval_add_S_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_K, crvec v,
                                rvec out, rvec work) const;
    /// Size of the workspace required by @ref eval_add_R_masked() and
    /// @ref eval_add_R_prod_masked().
    [[nodiscard]] length_t get_R_work_size() const;
    /// Size of the workspace required by @ref eval_add_S_masked() and
    /// @ref eval_add_S_prod_masked().
    [[nodiscard]] length_t get_S_work_size() const;

    /// @}

    /// @name Constraints
    /// @{

    /// Stage constraints @f$ c_k(x^k) @f$.
    void eval_constr(index_t timestep, crvec x, rvec c) const;
    /// Terminal constraints @f$ c_N(x^N) @f$.
    void eval_constr_N(crvec x, rvec c) const;
    /// Gradient-vector product of stage constraints @f$ \nabla c_k(x^k)\, p @f$.
    void eval_grad_constr_prod(index_t timestep, crvec x, crvec p, rvec grad_cx_p) const;
    /// Gradient-vector product of terminal constraints @f$ \nabla c_N(x^N)\, p @f$.
    void eval_grad_constr_prod_N(crvec x, crvec p, rvec grad_cx_p) const;
    /// Gauss-Newton Hessian of stage constraints @f$ \tp{\jac_{c_k}}(x^k)\,
    /// \operatorname{diag}(M)\; \jac_{c_k}(x^k) @f$.
    void eval_add_gn_hess_constr(index_t timestep, crvec x, crvec M, rmat out) const;
    /// Gauss-Newton Hessian of terminal constraints @f$ \tp{\jac_{c_N}}(x^N)\,
    /// \operatorname{diag}(M)\; \jac_{c_N}(x^N) @f$.
    void eval_add_gn_hess_constr_N(crvec x, crvec M, rmat out) const;

    /// @}

    /// @name Checks
    /// @{

    /// Check that the problem formulation is well-defined, the dimensions match,
    /// etc. Throws an exception if this is not the case.
    void check() const;

    /// @}

    /// @name Querying specialized implementations
    /// @{

    // clang-format off
    [[nodiscard]] bool provides_get_D() const { return vtable.get_D != nullptr; }
    [[nodiscard]] bool provides_get_D_N() const { return vtable.get_D_N != &vtable.default_get_D_N; }
    [[nodiscard]] bool provides_eval_h() const { return vtable.eval_h != nullptr; }
    [[nodiscard]] bool provides_eval_h_N() const { return vtable.eval_h_N != nullptr; }
    [[nodiscard]] bool provides_eval_add_Q_N() const { return vtable.eval_add_Q_N != &vtable.default_eval_add_Q_N; }
    [[nodiscard]] bool provides_eval_add_R_prod_masked() const { return vtable.eval_add_R_prod_masked != &vtable.default_eval_add_R_prod_masked; }
    [[nodiscard]] bool provides_eval_add_S_prod_masked() const { return vtable.eval_add_S_prod_masked != &vtable.default_eval_add_S_prod_masked; }
    [[nodiscard]] bool provides_get_R_work_size() const { return vtable.get_R_work_size != &vtable.default_get_R_work_size; }
    [[nodiscard]] bool provides_get_S_work_size() const { return vtable.get_S_work_size != &vtable.default_get_S_work_size; }
    [[nodiscard]] bool provides_eval_constr() const { return vtable.eval_constr != nullptr; }
    [[nodiscard]] bool provides_eval_constr_N() const { return vtable.eval_constr_N != &vtable.default_eval_constr_N; }
    [[nodiscard]] bool provides_eval_grad_constr_prod() const { return vtable.eval_grad_constr_prod != nullptr; }
    [[nodiscard]] bool provides_eval_grad_constr_prod_N() const { return vtable.eval_grad_constr_prod_N != &vtable.default_eval_grad_constr_prod_N; }
    [[nodiscard]] bool provides_eval_add_gn_hess_constr() const { return vtable.eval_add_gn_hess_constr != nullptr; }
    [[nodiscard]] bool provides_eval_add_gn_hess_constr_N() const { return vtable.eval_add_gn_hess_constr_N != &vtable.default_eval_add_gn_hess_constr_N; }
    // clang-format on

    /// @}
};

// clang-format off
#ifdef NDEBUG
[[gnu::always_inline]] inline void check_finiteness(auto &&, auto &&) {}
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_proj_diff_g(crvec z, rvec e) const { return call(vtable.eval_proj_diff_g, z, e); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_proj_multipliers(rvec y, real_t M) const { return call(vtable.eval_proj_multipliers, y, M); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::get_U(Box &U) const { return call(vtable.get_U, U); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::get_D(Box &D) const { return call(vtable.get_D, D); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::get_D_N(Box &D) const { return call(vtable.get_D_N, D); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::get_x_init(rvec x_init) const { return call(vtable.get_x_init, x_init); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const { return call(vtable.eval_f, timestep, x, u, fxu); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const { return call(vtable.eval_jac_f, timestep, x, u, J_fxu); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p, rvec grad_fxu_p) const { return call(vtable.eval_grad_f_prod, timestep, x, u, p, grad_fxu_p); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_h(index_t timestep, crvec x, crvec u, rvec h) const { return call(vtable.eval_h, timestep, x, u, h); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_h_N(crvec x, rvec h) const { return call(vtable.eval_h_N, x, h); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline auto TypeErasedControlProblem<Conf, Allocator>::eval_l(index_t timestep, crvec h) const -> real_t { return call(vtable.eval_l, timestep, h); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline auto TypeErasedControlProblem<Conf, Allocator>::eval_l_N(crvec h) const -> real_t { return call(vtable.eval_l_N, h); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_qr(index_t timestep, crvec xu, crvec h, rvec qr) const { return call(vtable.eval_qr, timestep, xu, h, qr); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_q_N(crvec x, crvec h, rvec q) const { return call(vtable.eval_q_N, x, h, q); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_Q(index_t timestep, crvec xu, crvec h, rmat Q) const { return call(vtable.eval_add_Q, timestep, xu, h, Q); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_Q_N(crvec x, crvec h, rmat Q) const { return call(vtable.eval_add_Q_N, x, h, Q); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_R_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat R, rvec work) const { return call(vtable.eval_add_R_masked, timestep, xu, h, mask, R, work); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_S_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat S, rvec work) const { return call(vtable.eval_add_S_masked, timestep, xu, h, mask, S, work); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_R_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_J, crindexvec mask_K, crvec v, rvec out, rvec work) const { return call(vtable.eval_add_R_prod_masked, timestep, xu, h, mask_J, mask_K, v, out, work); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_S_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_K, crvec v, rvec out, rvec work) const { return call(vtable.eval_add_S_prod_masked, timestep, xu, h, mask_K, v, out, work); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline auto TypeErasedControlProblem<Conf, Allocator>::get_R_work_size() const -> length_t { return call(vtable.get_R_work_size); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline auto TypeErasedControlProblem<Conf, Allocator>::get_S_work_size() const -> length_t { return call(vtable.get_S_work_size); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_constr(index_t timestep, crvec x, rvec c) const { return call(vtable.eval_constr, timestep, x, c); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_constr_N(crvec x, rvec c) const { return call(vtable.eval_constr_N, x, c); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_grad_constr_prod(index_t timestep, crvec x, crvec p, rvec grad_cx_p) const { return call(vtable.eval_grad_constr_prod, timestep, x, p, grad_cx_p); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_grad_constr_prod_N(crvec x, crvec p, rvec grad_cx_p) const { return call(vtable.eval_grad_constr_prod_N, x, p, grad_cx_p); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_gn_hess_constr(index_t timestep, crvec x, crvec M, rmat out) const { return call(vtable.eval_add_gn_hess_constr, timestep, x, M, out); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_gn_hess_constr_N(crvec x, crvec M, rmat out) const { return call(vtable.eval_add_gn_hess_constr_N, x, M, out); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::check() const { return call(vtable.check); }
#else
/// If the given vector @p v is not finite, break or throw an exception with the
/// given message @p msg.
inline void check_finiteness(const auto &v, std::string_view msg) {
    using std::begin;
    using std::end;
    if (!v.allFinite()) {
        std::cout << msg << std::endl;
        throw std::runtime_error(std::string(msg));
    }
}
inline void check_finiteness(const std::floating_point auto &v, std::string_view msg) {
    if (!std::isfinite(v)) {
        std::cout << msg << std::endl;
        throw std::runtime_error(std::string(msg));
    }
}
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_proj_diff_g(crvec z, rvec e) const { return call(vtable.eval_proj_diff_g, z, e); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_proj_multipliers(rvec y, real_t M) const { return call(vtable.eval_proj_multipliers, y, M); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::get_U(Box &U) const { return call(vtable.get_U, U); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::get_D(Box &D) const { return call(vtable.get_D, D); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::get_D_N(Box &D_N) const { return call(vtable.get_D_N, D_N); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::get_x_init(rvec x_init) const { call(vtable.get_x_init, x_init); check_finiteness(x_init, "Infinite output of get_x_init"); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const { check_finiteness(x, "Infinite input x of f"); check_finiteness(u, "Infinite input u of f"); call(vtable.eval_f, timestep, x, u, fxu); check_finiteness(fxu, "Infinite output of f"); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const { check_finiteness(x, "Infinite input x of jac_f"); check_finiteness(u, "Infinite input u of jac_f"); call(vtable.eval_jac_f, timestep, x, u, J_fxu); check_finiteness(J_fxu.reshaped(), "Infinite output of jac_f"); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p, rvec grad_fxu_p) const { check_finiteness(x, "Infinite input x of grad_f_prod"); check_finiteness(u, "Infinite input u of grad_f_prod"); check_finiteness(p, "Infinite input p of grad_f_prod"); call(vtable.eval_grad_f_prod, timestep, x, u, p, grad_fxu_p); check_finiteness(grad_fxu_p, "Infinite output of jac_f"); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_h(index_t timestep, crvec x, crvec u, rvec h) const { check_finiteness(x, "Infinite input x of h"); check_finiteness(u, "Infinite input u of h"); call(vtable.eval_h, timestep, x, u, h); check_finiteness(h, "Infinite output of h"); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_h_N(crvec x, rvec h) const { check_finiteness(x, "Infinite input x of h_N"); call(vtable.eval_h_N, x, h); check_finiteness(h, "Infinite output of h_N"); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline auto TypeErasedControlProblem<Conf, Allocator>::eval_l(index_t timestep, crvec h) const -> real_t { check_finiteness(h, "Infinite input h of l"); auto l = call(vtable.eval_l, timestep, h); check_finiteness(l, "Infinite output of l"); return l; }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline auto TypeErasedControlProblem<Conf, Allocator>::eval_l_N(crvec h) const -> real_t { check_finiteness(h, "Infinite input h of l_N"); auto l =  call(vtable.eval_l_N, h); check_finiteness(l, "Infinite output of l_N"); return l; }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_qr(index_t timestep, crvec xu, crvec h, rvec qr) const { return call(vtable.eval_qr, timestep, xu, h, qr); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_q_N(crvec x, crvec h, rvec q) const { check_finiteness(x, "Infinite input x of q_N"); check_finiteness(h, "Infinite input h of q_N"); call(vtable.eval_q_N, x, h, q); check_finiteness(q, "Infinite output of q_N"); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_Q(index_t timestep, crvec xu, crvec h, rmat Q) const { return call(vtable.eval_add_Q, timestep, xu, h, Q); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_Q_N(crvec x, crvec h, rmat Q) const { return call(vtable.eval_add_Q_N, x, h, Q); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_R_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat R, rvec work) const { return call(vtable.eval_add_R_masked, timestep, xu, h, mask, R, work); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_S_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat S, rvec work) const { return call(vtable.eval_add_S_masked, timestep, xu, h, mask, S, work); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_R_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_J, crindexvec mask_K, crvec v, rvec out, rvec work) const { return call(vtable.eval_add_R_prod_masked, timestep, xu, h, mask_J, mask_K, v, out, work); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_S_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_K, crvec v, rvec out, rvec work) const { return call(vtable.eval_add_S_prod_masked, timestep, xu, h, mask_K, v, out, work); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline auto TypeErasedControlProblem<Conf, Allocator>::get_R_work_size() const -> length_t { return call(vtable.get_R_work_size); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline auto TypeErasedControlProblem<Conf, Allocator>::get_S_work_size() const -> length_t { return call(vtable.get_S_work_size); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_constr(index_t timestep, crvec x, rvec c) const { return call(vtable.eval_constr, timestep, x, c); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_constr_N(crvec x, rvec c) const { return call(vtable.eval_constr_N, x, c); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_grad_constr_prod(index_t timestep, crvec x, crvec p, rvec grad_cx_p) const { return call(vtable.eval_grad_constr_prod, timestep, x, p, grad_cx_p); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_grad_constr_prod_N(crvec x, crvec p, rvec grad_cx_p) const { return call(vtable.eval_grad_constr_prod_N, x, p, grad_cx_p); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_gn_hess_constr(index_t timestep, crvec x, crvec M, rmat out) const { return call(vtable.eval_add_gn_hess_constr, timestep, x, M, out); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::eval_add_gn_hess_constr_N(crvec x, crvec M, rmat out) const { return call(vtable.eval_add_gn_hess_constr_N, x, M, out); }
template <Config Conf, class Allocator> [[gnu::always_inline]] inline void TypeErasedControlProblem<Conf, Allocator>::check() const { return call(vtable.check); }
#endif
// clang-format on

template <class Problem>
struct ControlProblemWithCounters {
    USING_ALPAQA_CONFIG_TEMPLATE(std::remove_cvref_t<Problem>::config_t);
    using Box = typename TypeErasedControlProblem<config_t>::Box;

    [[nodiscard, gnu::always_inline]] length_t get_N() const { return problem.get_N(); }
    [[nodiscard, gnu::always_inline]] length_t get_nu() const { return problem.get_nu(); }
    [[nodiscard, gnu::always_inline]] length_t get_nx() const { return problem.get_nx(); }
    [[nodiscard, gnu::always_inline]] length_t get_nh() const { return problem.get_nh(); }
    [[nodiscard, gnu::always_inline]] length_t get_nh_N() const { return problem.get_nh_N(); }
    [[nodiscard, gnu::always_inline]] length_t get_nc() const { return problem.get_nc(); }
    [[nodiscard, gnu::always_inline]] length_t get_nc_N() const { return problem.get_nc_N(); }

    // clang-format off
    [[gnu::always_inline]] void eval_proj_diff_g(crvec z, rvec e) const { return problem.eval_proj_diff_g(z, e); }
    [[gnu::always_inline]] void eval_proj_multipliers(rvec y, real_t M) const { return problem.eval_proj_multipliers(y, M); }
    [[gnu::always_inline]] void get_x_init(rvec x_init) const { return problem.get_x_init(x_init); }
    [[nodiscard, gnu::always_inline]] length_t get_R_work_size() const requires requires { &std::remove_cvref_t<Problem>::get_R_work_size; } { return problem.get_R_work_size(); }
    [[nodiscard, gnu::always_inline]] length_t get_S_work_size() const requires requires { &std::remove_cvref_t<Problem>::get_S_work_size; } { return problem.get_S_work_size(); }
    [[gnu::always_inline]] void get_U(Box &U) const requires requires { &std::remove_cvref_t<Problem>::get_U; } { return problem.get_U(U); }
    [[gnu::always_inline]] void get_D(Box &D) const requires requires { &std::remove_cvref_t<Problem>::get_D; } { return problem.get_D(D); }
    [[gnu::always_inline]] void get_D_N(Box &D) const requires requires { &std::remove_cvref_t<Problem>::get_D_N; } { return problem.get_D_N(D); }
    [[gnu::always_inline]] void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const { ++evaluations->f; return timed(evaluations->time.f, std::bind(&std::remove_cvref_t<Problem>::eval_f, &problem, timestep, x, u, fxu)); }
    [[gnu::always_inline]] void eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const { ++evaluations->jac_f; return timed(evaluations->time.jac_f, std::bind(&std::remove_cvref_t<Problem>::eval_jac_f, &problem, timestep, x, u, J_fxu)); }
    [[gnu::always_inline]] void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p, rvec grad_fxu_p) const { ++evaluations->grad_f_prod; return timed(evaluations->time.grad_f_prod, std::bind(&std::remove_cvref_t<Problem>::eval_grad_f_prod, &problem, timestep, x, u, p, grad_fxu_p)); }
    [[gnu::always_inline]] void eval_h(index_t timestep, crvec x, crvec u, rvec h) const { ++evaluations->h; return timed(evaluations->time.h, std::bind(&std::remove_cvref_t<Problem>::eval_h, &problem, timestep, x, u, h)); }
    [[gnu::always_inline]] void eval_h_N(crvec x, rvec h) const { ++evaluations->h_N; return timed(evaluations->time.h_N, std::bind(&std::remove_cvref_t<Problem>::eval_h_N, &problem, x, h)); }
    [[nodiscard, gnu::always_inline]] real_t eval_l(index_t timestep, crvec h) const { ++evaluations->l; return timed(evaluations->time.l, std::bind(&std::remove_cvref_t<Problem>::eval_l, &problem, timestep, h)); }
    [[nodiscard, gnu::always_inline]] real_t eval_l_N(crvec h) const { ++evaluations->l_N; return timed(evaluations->time.l_N, std::bind(&std::remove_cvref_t<Problem>::eval_l_N, &problem, h)); }
    [[gnu::always_inline]] void eval_qr(index_t timestep, crvec xu, crvec h, rvec qr) const { ++evaluations->qr; return timed(evaluations->time.qr, std::bind(&std::remove_cvref_t<Problem>::eval_qr, &problem, timestep, xu, h, qr)); }
    [[gnu::always_inline]] void eval_q_N(crvec x, crvec h, rvec q) const requires requires { &std::remove_cvref_t<Problem>::eval_q_N; } { ++evaluations->q_N; return timed(evaluations->time.q_N, std::bind(&std::remove_cvref_t<Problem>::eval_q_N, &problem, x, h, q)); }
    [[gnu::always_inline]] void eval_add_Q(index_t timestep, crvec xu, crvec h, rmat Q) const { ++evaluations->add_Q; return timed(evaluations->time.add_Q, std::bind(&std::remove_cvref_t<Problem>::eval_add_Q, &problem, timestep, xu, h, Q)); }
    [[gnu::always_inline]] void eval_add_Q_N(crvec x, crvec h, rmat Q) const requires requires { &std::remove_cvref_t<Problem>::eval_add_Q_N; } { ++evaluations->add_Q_N; return timed(evaluations->time.add_Q_N, std::bind(&std::remove_cvref_t<Problem>::eval_add_Q_N, &problem, x, h, Q)); }
    [[gnu::always_inline]] void eval_add_R_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat R, rvec work) const { ++evaluations->add_R_masked; return timed(evaluations->time.add_R_masked, std::bind(&std::remove_cvref_t<Problem>::eval_add_R_masked, &problem, timestep, xu, h, mask, R, work)); }
    [[gnu::always_inline]] void eval_add_S_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat S, rvec work) const { ++evaluations->add_S_masked; return timed(evaluations->time.add_S_masked, std::bind(&std::remove_cvref_t<Problem>::eval_add_S_masked, &problem, timestep, xu, h, mask, S, work)); }
    [[gnu::always_inline]] void eval_add_R_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_J, crindexvec mask_K, crvec v, rvec out, rvec work) const requires requires { &std::remove_cvref_t<Problem>::eval_add_R_prod_masked; } { ++evaluations->add_R_prod_masked; return timed(evaluations->time.add_R_prod_masked, std::bind(&std::remove_cvref_t<Problem>::eval_add_R_prod_masked, &problem, timestep, xu, h, mask_J, mask_K, v, out, work)); }
    [[gnu::always_inline]] void eval_add_S_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_K, crvec v, rvec out, rvec work) const requires requires { &std::remove_cvref_t<Problem>::eval_add_S_prod_masked; } { ++evaluations->add_S_prod_masked; return timed(evaluations->time.add_S_prod_masked, std::bind(&std::remove_cvref_t<Problem>::eval_add_S_prod_masked, &problem, timestep, xu, h, mask_K, v, out, work)); }
    [[gnu::always_inline]] void eval_constr(index_t timestep, crvec x, rvec c) const requires requires { &std::remove_cvref_t<Problem>::eval_constr; } { ++evaluations->constr; return timed(evaluations->time.constr, std::bind(&std::remove_cvref_t<Problem>::eval_constr, &problem, timestep, x, c)); }
    [[gnu::always_inline]] void eval_constr_N(crvec x, rvec c) const requires requires { &std::remove_cvref_t<Problem>::eval_constr_N; } { ++evaluations->constr_N; return timed(evaluations->time.constr_N, std::bind(&std::remove_cvref_t<Problem>::eval_constr_N, &problem, x, c)); }
    [[gnu::always_inline]] void eval_grad_constr_prod(index_t timestep, crvec x, crvec p, rvec grad_cx_p) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_constr_prod; } { ++evaluations->grad_constr_prod; return timed(evaluations->time.grad_constr_prod, std::bind(&std::remove_cvref_t<Problem>::eval_grad_constr_prod, &problem, timestep, x, p, grad_cx_p)); }
    [[gnu::always_inline]] void eval_grad_constr_prod_N(crvec x, crvec p, rvec grad_cx_p) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_constr_prod_N; } { ++evaluations->grad_constr_prod_N; return timed(evaluations->time.grad_constr_prod_N, std::bind(&std::remove_cvref_t<Problem>::eval_grad_constr_prod_N, &problem, x, p, grad_cx_p)); }
    [[gnu::always_inline]] void eval_add_gn_hess_constr(index_t timestep, crvec x, crvec M, rmat out) const requires requires { &std::remove_cvref_t<Problem>::eval_add_gn_hess_constr; } { ++evaluations->add_gn_hess_constr; return timed(evaluations->time.add_gn_hess_constr, std::bind(&std::remove_cvref_t<Problem>::eval_add_gn_hess_constr, &problem, timestep, x, M, out)); }
    [[gnu::always_inline]] void eval_add_gn_hess_constr_N(crvec x, crvec M, rmat out) const requires requires { &std::remove_cvref_t<Problem>::eval_add_gn_hess_constr_N; } { ++evaluations->add_gn_hess_constr_N; return timed(evaluations->time.add_gn_hess_constr_N, std::bind(&std::remove_cvref_t<Problem>::eval_add_gn_hess_constr_N, &problem, x, M, out)); }
    [[gnu::always_inline]] void check() const { problem.check(); }

    [[nodiscard]] bool provides_get_D() const requires requires (Problem p) { { p.provides_get_D() } -> std::convertible_to<bool>; } { return problem.provides_get_D(); }
    [[nodiscard]] bool provides_get_D_N() const requires requires (Problem p) { { p.provides_get_D_N() } -> std::convertible_to<bool>; } { return problem.provides_get_D_N(); }
    [[nodiscard]] bool provides_eval_add_Q_N() const requires requires (Problem p) { { p.provides_eval_add_Q_N() } -> std::convertible_to<bool>; } { return problem.provides_eval_add_Q_N(); }
    [[nodiscard]] bool provides_eval_add_R_prod_masked() const requires requires (Problem p) { { p.provides_eval_add_R_prod_masked() } -> std::convertible_to<bool>; } { return problem.provides_eval_add_R_prod_masked(); }
    [[nodiscard]] bool provides_eval_add_S_prod_masked() const requires requires (Problem p) { { p.provides_eval_add_S_prod_masked() } -> std::convertible_to<bool>; } { return problem.provides_eval_add_S_prod_masked(); }
    [[nodiscard]] bool provides_get_R_work_size() const requires requires (Problem p) { { p.provides_get_R_work_size() } -> std::convertible_to<bool>; } { return problem.provides_get_R_work_size(); }
    [[nodiscard]] bool provides_get_S_work_size() const requires requires (Problem p) { { p.provides_get_S_work_size() } -> std::convertible_to<bool>; } { return problem.provides_get_S_work_size(); }
    [[nodiscard]] bool provides_eval_constr() const requires requires (Problem p) { { p.provides_eval_constr() } -> std::convertible_to<bool>; } { return problem.provides_eval_constr(); }
    [[nodiscard]] bool provides_eval_constr_N() const requires requires (Problem p) { { p.provides_eval_constr_N() } -> std::convertible_to<bool>; } { return problem.provides_eval_constr_N(); }
    [[nodiscard]] bool provides_eval_grad_constr_prod() const requires requires (Problem p) { { p.provides_eval_grad_constr_prod() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_constr_prod(); }
    [[nodiscard]] bool provides_eval_grad_constr_prod_N() const requires requires (Problem p) { { p.provides_eval_grad_constr_prod_N() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_constr_prod_N(); }
    [[nodiscard]] bool provides_eval_add_gn_hess_constr() const requires requires (Problem p) { { p.provides_eval_add_gn_hess_constr() } -> std::convertible_to<bool>; } { return problem.provides_eval_add_gn_hess_constr(); }
    [[nodiscard]] bool provides_eval_add_gn_hess_constr_N() const requires requires (Problem p) { { p.provides_eval_add_gn_hess_constr_N() } -> std::convertible_to<bool>; } { return problem.provides_eval_add_gn_hess_constr_N(); }
    // clang-format on

    std::shared_ptr<OCPEvalCounter> evaluations = std::make_shared<OCPEvalCounter>();
    Problem problem;

    ControlProblemWithCounters()
        requires std::is_default_constructible_v<Problem>
    = default;
    template <class P>
    explicit ControlProblemWithCounters(P &&problem)
        requires std::is_same_v<std::remove_cvref_t<P>, std::remove_cvref_t<Problem>>
        : problem{std::forward<P>(problem)} {}
    template <class... Args>
    explicit ControlProblemWithCounters(std::in_place_t, Args &&...args)
        requires(!std::is_lvalue_reference_v<Problem>)
        : problem{std::forward<Args>(args)...} {}

    /// Reset all evaluation counters and timers to zero. Affects all instances
    /// that share the same evaluations. If you only want to reset the counters
    /// of this instance, use @ref decouple_evaluations first.
    void reset_evaluations() { evaluations.reset(); }
    /// Give this instance its own evaluation counters and timers, decoupling
    /// it from any other instances they might have previously been shared with.
    /// The evaluation counters and timers are preserved (a copy is made).
    void decouple_evaluations() { evaluations = std::make_shared<OCPEvalCounter>(*evaluations); }

  private:
    template <class TimeT, class FunT>
    [[gnu::always_inline]] static decltype(auto) timed(TimeT &time, FunT &&f) {
        alpaqa::util::Timed timed{time};
        return std::forward<FunT>(f)();
    }
};

template <class Problem>
[[nodiscard]] auto ocproblem_with_counters(Problem &&p) {
    using Prob        = std::remove_cvref_t<Problem>;
    using ProbWithCnt = ControlProblemWithCounters<Prob>;
    return ProbWithCnt{std::forward<Problem>(p)};
}

template <class Problem>
[[nodiscard]] auto ocproblem_with_counters_ref(Problem &p) {
    using Prob        = std::remove_cvref_t<Problem>;
    using ProbWithCnt = ControlProblemWithCounters<const Prob &>;
    return ProbWithCnt{p};
}

} // namespace alpaqa