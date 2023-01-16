#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/box.hpp>
#include <alpaqa/problem/ocproblem-counters.hpp>
#include <alpaqa/problem/problem-counters.hpp>
#include <alpaqa/util/required-method.hpp>
#include <alpaqa/util/type-erasure.hpp>
#include <array>
#include <stdexcept>

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
struct OCProblemVTable : util::BasicVTable {
    USING_ALPAQA_CONFIG(Conf);
    using Box = alpaqa::Box<config_t>;

    template <class F>
    using optional_function_t = util::BasicVTable::optional_function_t<F, OCProblemVTable>;
    template <class F>
    using optional_const_function_t =
        util::BasicVTable::optional_const_function_t<F, OCProblemVTable>;

    // clang-format off
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
    required_const_function_t<void(index_t timestep, crvec x, crvec u, rvec h)>
        eval_h;
    required_const_function_t<void(crvec x, rvec h)>
        eval_h_N;
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
    optional_const_function_t<void(index_t timestep, crvec x, rmat J_c)>
        eval_jac_constr = nullptr;
    optional_const_function_t<void(crvec x, rmat J_c)>
        eval_jac_constr_N = &default_eval_jac_constr_N;
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
    OCProblemVTable(util::VTableTypeTag<P> t) : util::BasicVTable{t} {
        ALPAQA_TE_REQUIRED_METHOD(*this, P, get_U);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, get_D, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, get_D_N, t.t);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, get_x_init);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_f);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_jac_f);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_grad_f_prod);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_h);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_h_N);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_l);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_l_N);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_qr);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_q_N);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_add_Q);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_add_Q_N, t.t);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_add_R_masked);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, eval_add_S_masked);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_add_R_prod_masked, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_add_S_prod_masked, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, get_R_work_size, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, get_S_work_size, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_constr, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_constr_N, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_jac_constr, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_jac_constr_N, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_grad_constr_prod, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_grad_constr_prod_N, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_add_gn_hess_constr, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(*this, P, eval_add_gn_hess_constr_N, t.t);
        ALPAQA_TE_REQUIRED_METHOD(*this, P, check);
        N    = t.t->get_N();
        nu   = t.t->get_nu();
        nx   = t.t->get_nx();
        nh   = t.t->get_nh();
        nh_N = t.t->get_nh_N();
        nc   = t.t->get_nc();
        nc_N = t.t->get_nc_N();
        if (nc > 0 && get_D == nullptr)
            throw std::runtime_error("OCProblem: missing 'get_D'");
        if (nc > 0 && eval_constr == nullptr)
            throw std::runtime_error("OCProblem: missing 'eval_constr'");
        if (nc > 0 && eval_jac_constr == nullptr)
            throw std::runtime_error("OCProblem: missing 'eval_jac_constr'");
        if (nc > 0 && eval_grad_constr_prod == nullptr)
            throw std::runtime_error("OCProblem: missing 'eval_grad_constr_prod'");
        if (nc > 0 && eval_add_gn_hess_constr == nullptr)
            throw std::runtime_error("OCProblem: missing 'eval_add_gn_hess_constr'");
    }
    OCProblemVTable() = default;

    static void default_get_D_N(const void *self, Box &D, const OCProblemVTable &vtable) {
        vtable.get_D(self, D, vtable);
    }
    static void default_eval_add_Q_N(const void *self, crvec x, crvec h, rmat Q,
                                     const OCProblemVTable &vtable) {
        vtable.eval_add_Q(self, vtable.N, x, h, Q);
    }
    static void default_eval_add_R_prod_masked(const void *, index_t, crvec, crvec, crindexvec,
                                               crindexvec, crvec, rvec, rvec,
                                               const OCProblemVTable &) {}
    static void default_eval_add_S_prod_masked(const void *, index_t, crvec, crvec, crindexvec,
                                               crvec, rvec, rvec, const OCProblemVTable &) {}
    [[nodiscard]] static length_t default_get_R_work_size(const void *, const OCProblemVTable &) {
        return 0;
    }
    [[nodiscard]] static length_t default_get_S_work_size(const void *, const OCProblemVTable &) {
        return 0;
    }
    static void default_eval_constr_N(const void *self, crvec x, rvec c,
                                      const OCProblemVTable &vtable) {
        vtable.eval_constr(self, vtable.N, x, c, vtable);
    }
    static void default_eval_jac_constr_N(const void *self, crvec x, rmat J_c,
                                          const OCProblemVTable &vtable) {
        vtable.eval_jac_constr(self, vtable.N, x, J_c, vtable);
    }
    static void default_eval_grad_constr_prod_N(const void *self, crvec x, crvec p, rvec grad_cx_p,
                                                const OCProblemVTable &vtable) {
        vtable.eval_grad_constr_prod(self, vtable.N, x, p, grad_cx_p, vtable);
    }
    static void default_eval_add_gn_hess_constr_N(const void *self, crvec x, crvec M, rmat out,
                                                  const OCProblemVTable &vtable) {
        vtable.eval_add_gn_hess_constr(self, vtable.N, x, M, out, vtable);
    }
};

template <Config Conf = DefaultConfig, class Allocator = std::allocator<std::byte>>
class TypeErasedOCProblem : public util::TypeErased<OCProblemVTable<Conf>, Allocator> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using VTable         = OCProblemVTable<config_t>;
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
    static TypeErasedOCProblem make(Args &&...args) {
        return TypeErased::template make<TypeErasedOCProblem, T>(std::forward<Args>(args)...);
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

    /// @}

    void get_U(Box &U) const;
    void get_D(Box &D) const;
    void get_D_N(Box &D) const;
    void get_x_init(rvec x_init) const;
    void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const;
    void eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const;
    void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p, rvec grad_fxu_p) const;
    void eval_h(index_t timestep, crvec x, crvec u, rvec h) const;
    void eval_h_N(crvec x, rvec h) const;
    [[nodiscard]] real_t eval_l(index_t timestep, crvec h) const;
    [[nodiscard]] real_t eval_l_N(crvec h) const;
    void eval_qr(index_t timestep, crvec xu, crvec h, rvec qr) const;
    void eval_q_N(crvec x, crvec h, rvec q) const;
    void eval_add_Q(index_t timestep, crvec xu, crvec h, rmat Q) const;
    void eval_add_Q_N(crvec x, crvec h, rmat Q) const;
    void eval_add_R_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat R,
                           rvec work) const;
    void eval_add_S_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat S,
                           rvec work) const;
    void eval_add_R_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_J,
                                crindexvec mask_K, crvec v, rvec out, rvec work) const;
    void eval_add_S_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_K, crvec v,
                                rvec out, rvec work) const;
    [[nodiscard]] length_t get_R_work_size() const;
    [[nodiscard]] length_t get_S_work_size() const;
    void eval_constr(index_t timestep, crvec x, rvec c) const;
    void eval_constr_N(crvec x, rvec c) const;
    void eval_jac_constr(index_t timestep, crvec x, rmat J_c) const;
    void eval_jac_constr_N(crvec x, rmat J_c) const;
    void eval_grad_constr_prod(index_t timestep, crvec x, crvec p, rvec grad_cx_p) const;
    void eval_grad_constr_prod_N(crvec x, crvec p, rvec grad_cx_p) const;
    void eval_add_gn_hess_constr(index_t timestep, crvec x, crvec M, rmat out) const;
    void eval_add_gn_hess_constr_N(crvec x, crvec M, rmat out) const;
    void check() const;
};

// clang-format off
#ifdef NDEBUG
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::get_U(Box &U) const { return call(vtable.get_U, U); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::get_D(Box &D) const { return call(vtable.get_D, D); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::get_D_N(Box &D) const { return call(vtable.get_D_N, D); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::get_x_init(rvec x_init) const { return call(vtable.get_x_init, x_init); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const { return call(vtable.eval_f, timestep, x, u, fxu); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const { return call(vtable.eval_jac_f, timestep, x, u, J_fxu); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p, rvec grad_fxu_p) const { return call(vtable.eval_grad_f_prod, timestep, x, u, p, grad_fxu_p); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_h(index_t timestep, crvec x, crvec u, rvec h) const { return call(vtable.eval_h, timestep, x, u, h); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_h_N(crvec x, rvec h) const { return call(vtable.eval_h_N, x, h); }
template <Config Conf, class Allocator> auto TypeErasedOCProblem<Conf, Allocator>::eval_l(index_t timestep, crvec h) const -> real_t { return call(vtable.eval_l, timestep, h); }
template <Config Conf, class Allocator> auto TypeErasedOCProblem<Conf, Allocator>::eval_l_N(crvec h) const -> real_t { return call(vtable.eval_l_N, h); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_qr(index_t timestep, crvec xu, crvec h, rvec qr) const { return call(vtable.eval_qr, timestep, xu, h, qr); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_q_N(crvec x, crvec h, rvec q) const { return call(vtable.eval_q_N, x, h, q); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_Q(index_t timestep, crvec xu, crvec h, rmat Q) const { return call(vtable.eval_add_Q, timestep, xu, h, Q); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_Q_N(crvec x, crvec h, rmat Q) const { return call(vtable.eval_add_Q_N, x, h, Q); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_R_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat R, rvec work) const { return call(vtable.eval_add_R_masked, timestep, xu, h, mask, R, work); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_S_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat S, rvec work) const { return call(vtable.eval_add_S_masked, timestep, xu, h, mask, S, work); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_R_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_J, crindexvec mask_K, crvec v, rvec out, rvec work) const { return call(vtable.eval_add_R_prod_masked, timestep, xu, h, mask_J, mask_K, v, out, work); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_S_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_K, crvec v, rvec out, rvec work) const { return call(vtable.eval_add_S_prod_masked, timestep, xu, h, mask_K, v, out, work); }
template <Config Conf, class Allocator> auto TypeErasedOCProblem<Conf, Allocator>::get_R_work_size() const -> length_t { return call(vtable.get_R_work_size); }
template <Config Conf, class Allocator> auto TypeErasedOCProblem<Conf, Allocator>::get_S_work_size() const -> length_t { return call(vtable.get_S_work_size); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_constr(index_t timestep, crvec x, rvec c) const { return call(vtable.eval_constr, timestep, x, c); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_constr_N(crvec x, rvec c) const { return call(vtable.eval_constr_N, x, c); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_jac_constr(index_t timestep, crvec x, rmat J_c) const { return call(vtable.eval_jac_constr, timestep, x, J_c); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_jac_constr_N(crvec x, rmat J_c) const { return call(vtable.eval_jac_constr_N, x, J_c); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_grad_constr_prod(index_t timestep, crvec x, crvec p, rvec grad_cx_p) const { return call(vtable.eval_grad_constr_prod, timestep, x, p, grad_cx_p); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_grad_constr_prod_N(crvec x, crvec p, rvec grad_cx_p) const { return call(vtable.eval_grad_constr_prod_N, x, p, grad_cx_p); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_gn_hess_constr(index_t timestep, crvec x, crvec M, rmat out) const { return call(vtable.eval_add_gn_hess_constr, timestep, x, M, out); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_gn_hess_constr_N(crvec x, crvec M, rmat out) const { return call(vtable.eval_add_gn_hess_constr_N, x, M, out); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::check() const { return call(vtable.check); }
#else
/// If the given vector @p v is not finite, break or throw an exception with the
/// given message @p msg.
inline void check_finiteness(const auto &v, std::string_view msg) {
    using std::begin;
    using std::end;
    if (std::any_of(begin(v), end(v),
                    [](double d) { return !std::isfinite(d); })) {
        std::cout << msg << std::endl;
        throw std::runtime_error(std::string(msg));
    }
}
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::get_U(Box &U) const { return call(vtable.get_U, U); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::get_D(Box &D) const { return call(vtable.get_D, D); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::get_x_init(rvec x_init) const { call(vtable.get_x_init, x_init); check_finiteness(x_init, "Infinite output of get_x_init"); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const { check_finiteness(x, "Infinite input x of f"); check_finiteness(u, "Infinite input u of f"); call(vtable.eval_f, timestep, x, u, fxu); check_finiteness(fxu, "Infinite output of f"); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const { check_finiteness(x, "Infinite input x of jac_f"); check_finiteness(u, "Infinite input u of jac_f"); call(vtable.eval_jac_f, timestep, x, u, J_fxu); check_finiteness(J_fxu.reshaped(), "Infinite output of jac_f"); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p, rvec grad_fxu_p) const { check_finiteness(x, "Infinite input x of grad_f_prod"); check_finiteness(u, "Infinite input u of grad_f_prod"); check_finiteness(p, "Infinite input p of grad_prod_f"); call(vtable.eval_grad_f_prod, timestep, x, u, p, grad_fxu_p); check_finiteness(grad_fxu_p, "Infinite output of jac_f"); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_h(index_t timestep, crvec x, crvec u, rvec h) const { check_finiteness(x, "Infinite input x of h"); check_finiteness(u, "Infinite input u of h"); call(vtable.eval_h, timestep, x, u, h); check_finiteness(h, "Infinite output of h"); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_h_N(crvec x, rvec h) const { check_finiteness(x, "Infinite input x of h_N"); call(vtable.eval_h_N, x, h); check_finiteness(h, "Infinite output of h_N"); }
template <Config Conf, class Allocator> auto TypeErasedOCProblem<Conf, Allocator>::eval_l(index_t timestep, crvec h) const -> real_t { check_finiteness(h, "Infinite input h of l"); auto l = call(vtable.eval_l, timestep, h); check_finiteness(std::array{l}, "Infinite output of l"); return l; }
template <Config Conf, class Allocator> auto TypeErasedOCProblem<Conf, Allocator>::eval_l_N(crvec h) const -> real_t { check_finiteness(h, "Infinite input h of l_N"); auto l =  call(vtable.eval_l_N, h); check_finiteness(std::array{l}, "Infinite output of l_N"); return l; }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_qr(index_t timestep, crvec xu, crvec h, rvec qr) const { return call(vtable.eval_qr, timestep, xu, h, qr); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_q_N(crvec x, crvec h, rvec q) const { check_finiteness(x, "Infinite input x of q_N"); check_finiteness(h, "Infinite input h of q_N"); call(vtable.eval_q_N, x, h, q); check_finiteness(q, "Infinite output of q_N"); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_Q(index_t timestep, crvec xu, crvec h, rmat Q) const { return call(vtable.eval_add_Q, timestep, xu, h, Q); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_Q_N(crvec x, crvec h, rmat Q) const { return call(vtable.eval_add_Q_N, x, h, Q); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_R_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat R, rvec work) const { return call(vtable.eval_add_R_masked, timestep, xu, h, mask, R, work); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_S_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat S, rvec work) const { return call(vtable.eval_add_S_masked, timestep, xu, h, mask, S, work); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_R_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_J, crindexvec mask_K, crvec v, rvec out, rvec work) const { return call(vtable.eval_add_R_prod_masked, timestep, xu, h, mask_J, mask_K, v, out, work); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_S_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_K, crvec v, rvec out, rvec work) const { return call(vtable.eval_add_S_prod_masked, timestep, xu, h, mask_K, v, out, work); }
template <Config Conf, class Allocator> auto TypeErasedOCProblem<Conf, Allocator>::get_R_work_size() const -> length_t { return call(vtable.get_R_work_size); }
template <Config Conf, class Allocator> auto TypeErasedOCProblem<Conf, Allocator>::get_S_work_size() const -> length_t { return call(vtable.get_S_work_size); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_constr(index_t timestep, crvec x, rvec c) const { return call(vtable.eval_constr, timestep, x, c); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_constr_N(crvec x, rvec c) const { return call(vtable.eval_constr_N, x, c); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_jac_constr(index_t timestep, crvec x, rmat J_c) const { return call(vtable.eval_jac_constr, timestep, x, J_c); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_jac_constr_N(crvec x, rmat J_c) const { return call(vtable.eval_jac_constr_N, x, J_c); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_grad_constr_prod(index_t timestep, crvec x, crvec p, rvec grad_cx_p) const { return call(vtable.eval_grad_constr_prod, timestep, x, p, grad_cx_p); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_grad_constr_prod_N(crvec x, crvec p, rvec grad_cx_p) const { return call(vtable.eval_grad_constr_prod_N, x, p, grad_cx_p); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_gn_hess_constr(index_t timestep, crvec x, crvec M, rmat out) const { return call(vtable.eval_add_gn_hess_constr, timestep, x, M, out); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::eval_add_gn_hess_constr_N(crvec x, crvec M, rmat out) const { return call(vtable.eval_add_gn_hess_constr_N, x, M, out); }
template <Config Conf, class Allocator> void TypeErasedOCProblem<Conf, Allocator>::check() const { return call(vtable.check); }
#endif
// clang-format on

template <class Problem>
struct OCProblemWithCounters {
    USING_ALPAQA_CONFIG_TEMPLATE(std::remove_cvref_t<Problem>::config_t);
    using Box = typename TypeErasedOCProblem<config_t>::Box;

    [[nodiscard]] length_t get_N() const { return problem.get_N(); }
    [[nodiscard]] length_t get_nu() const { return problem.get_nu(); }
    [[nodiscard]] length_t get_nx() const { return problem.get_nx(); }
    [[nodiscard]] length_t get_nh() const { return problem.get_nh(); }
    [[nodiscard]] length_t get_nh_N() const { return problem.get_nh_N(); }
    [[nodiscard]] length_t get_nc() const { return problem.get_nc(); }
    [[nodiscard]] length_t get_nc_N() const { return problem.get_nc_N(); }

    // clang-format off
    void get_x_init(rvec x_init) const { return problem.get_x_init(x_init); }
    [[nodiscard]] length_t get_R_work_size() const requires requires { &std::remove_cvref_t<Problem>::get_R_work_size; } { return problem.get_R_work_size(); }
    [[nodiscard]] length_t get_S_work_size() const requires requires { &std::remove_cvref_t<Problem>::get_S_work_size; } { return problem.get_S_work_size(); }
    void get_U(Box &U) const requires requires { &std::remove_cvref_t<Problem>::get_U; } { return problem.get_U(U); }
    void get_D(Box &D) const requires requires { &std::remove_cvref_t<Problem>::get_D; } { return problem.get_D(D); }
    void get_D_N(Box &D) const requires requires { &std::remove_cvref_t<Problem>::get_D_N; } { return problem.get_D_N(D); }
    void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const { ++evaluations->f; return timed(evaluations->time.f, std::bind(&std::remove_cvref_t<Problem>::eval_f, &problem, timestep, x, u, fxu)); }
    void eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const { ++evaluations->jac_f; return timed(evaluations->time.jac_f, std::bind(&std::remove_cvref_t<Problem>::eval_jac_f, &problem, timestep, x, u, J_fxu)); }
    void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p, rvec grad_fxu_p) const { ++evaluations->grad_f_prod; return timed(evaluations->time.grad_f_prod, std::bind(&std::remove_cvref_t<Problem>::eval_grad_f_prod, &problem, timestep, x, u, p, grad_fxu_p)); }
    void eval_h(index_t timestep, crvec x, crvec u, rvec h) const { ++evaluations->h; return timed(evaluations->time.h, std::bind(&std::remove_cvref_t<Problem>::eval_h, &problem, timestep, x, u, h)); }
    void eval_h_N(crvec x, rvec h) const { ++evaluations->h_N; return timed(evaluations->time.h_N, std::bind(&std::remove_cvref_t<Problem>::eval_h_N, &problem, x, h)); }
    [[nodiscard]] real_t eval_l(index_t timestep, crvec h) const { ++evaluations->l; return timed(evaluations->time.l, std::bind(&std::remove_cvref_t<Problem>::eval_l, &problem, timestep, h)); }
    [[nodiscard]] real_t eval_l_N(crvec h) const { ++evaluations->l_N; return timed(evaluations->time.l_N, std::bind(&std::remove_cvref_t<Problem>::eval_l_N, &problem, h)); }
    void eval_qr(index_t timestep, crvec xu, crvec h, rvec qr) const { ++evaluations->qr; return timed(evaluations->time.qr, std::bind(&std::remove_cvref_t<Problem>::eval_qr, &problem, timestep, xu, h, qr)); }
    void eval_q_N(crvec x, crvec h, rvec q) const requires requires { &std::remove_cvref_t<Problem>::eval_q_N; } { ++evaluations->q_N; return timed(evaluations->time.q_N, std::bind(&std::remove_cvref_t<Problem>::eval_q_N, &problem, x, h, q)); }
    void eval_add_Q(index_t timestep, crvec xu, crvec h, rmat Q) const { ++evaluations->add_Q; return timed(evaluations->time.add_Q, std::bind(&std::remove_cvref_t<Problem>::eval_add_Q, &problem, timestep, xu, h, Q)); }
    void eval_add_Q_N(crvec x, crvec h, rmat Q) const requires requires { &std::remove_cvref_t<Problem>::eval_add_Q_N; } { ++evaluations->add_Q_N; return timed(evaluations->time.add_Q_N, std::bind(&std::remove_cvref_t<Problem>::eval_add_Q_N, &problem, x, h, Q)); }
    void eval_add_R_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat R, rvec work) const { ++evaluations->add_R_masked; return timed(evaluations->time.add_R_masked, std::bind(&std::remove_cvref_t<Problem>::eval_add_R_masked, &problem, timestep, xu, h, mask, R, work)); }
    void eval_add_S_masked(index_t timestep, crvec xu, crvec h, crindexvec mask, rmat S, rvec work) const { ++evaluations->add_S_masked; return timed(evaluations->time.add_S_masked, std::bind(&std::remove_cvref_t<Problem>::eval_add_S_masked, &problem, timestep, xu, h, mask, S, work)); }
    void eval_add_R_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_J, crindexvec mask_K, crvec v, rvec out, rvec work) const requires requires { &std::remove_cvref_t<Problem>::eval_add_R_prod_masked; } { ++evaluations->add_R_prod_masked; return timed(evaluations->time.add_R_prod_masked, std::bind(&std::remove_cvref_t<Problem>::eval_add_R_prod_masked, &problem, timestep, xu, h, mask_J, mask_K, v, out, work)); }
    void eval_add_S_prod_masked(index_t timestep, crvec xu, crvec h, crindexvec mask_K, crvec v, rvec out, rvec work) const requires requires { &std::remove_cvref_t<Problem>::eval_add_S_prod_masked; } { ++evaluations->add_S_prod_masked; return timed(evaluations->time.add_S_prod_masked, std::bind(&std::remove_cvref_t<Problem>::eval_add_S_prod_masked, &problem, timestep, xu, h, mask_K, v, out, work)); }
    void eval_constr(index_t timestep, crvec x, rvec c) const requires requires { &std::remove_cvref_t<Problem>::eval_constr; } { ++evaluations->constr; return timed(evaluations->time.constr, std::bind(&std::remove_cvref_t<Problem>::eval_constr, &problem, timestep, x, c)); }
    void eval_constr_N(crvec x, rvec c) const requires requires { &std::remove_cvref_t<Problem>::eval_constr_N; } { ++evaluations->constr_N; return timed(evaluations->time.constr_N, std::bind(&std::remove_cvref_t<Problem>::eval_constr_N, &problem, x, c)); }
    void eval_jac_constr(index_t timestep, crvec x, rmat J_c) const { ++evaluations->jac_constr; return timed(evaluations->time.jac_constr, std::bind(&std::remove_cvref_t<Problem>::eval_jac_constr, &problem, timestep, x, J_c)); }
    void eval_jac_constr_N(crvec x, rmat J_c) const requires requires { &std::remove_cvref_t<Problem>::eval_jac_constr_N; } { ++evaluations->jac_constr_N; return timed(evaluations->time.jac_constr_N, std::bind(&std::remove_cvref_t<Problem>::eval_jac_constr_N, &problem, x, J_c)); }
    void eval_grad_constr_prod(index_t timestep, crvec x, crvec p, rvec grad_cx_p) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_constr_prod; } { ++evaluations->grad_constr_prod; return timed(evaluations->time.grad_constr_prod, std::bind(&std::remove_cvref_t<Problem>::eval_grad_constr_prod, &problem, timestep, x, p, grad_cx_p)); }
    void eval_grad_constr_prod_N(crvec x, crvec p, rvec grad_cx_p) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_constr_prod_N; } { ++evaluations->grad_constr_prod_N; return timed(evaluations->time.grad_constr_prod_N, std::bind(&std::remove_cvref_t<Problem>::eval_grad_constr_prod_N, &problem, x, p, grad_cx_p)); }
    void eval_add_gn_hess_constr(index_t timestep, crvec x, crvec M, rmat out) const requires requires { &std::remove_cvref_t<Problem>::eval_add_gn_hess_constr; } { ++evaluations->add_gn_hess_constr; return timed(evaluations->time.add_gn_hess_constr, std::bind(&std::remove_cvref_t<Problem>::eval_add_gn_hess_constr, &problem, timestep, x, M, out)); }
    void eval_add_gn_hess_constr_N(crvec x, crvec M, rmat out) const requires requires { &std::remove_cvref_t<Problem>::eval_add_gn_hess_constr_N; } { ++evaluations->add_gn_hess_constr_N; return timed(evaluations->time.add_gn_hess_constr_N, std::bind(&std::remove_cvref_t<Problem>::eval_add_gn_hess_constr_N, &problem, x, M, out)); }
    void check() const { problem.check(); }

    [[nodiscard]] bool provides_get_D() const requires requires (Problem p) { { p.provides_get_D() } -> std::convertible_to<bool>; } { return problem.provides_get_D(); }
    [[nodiscard]] bool provides_get_D_N() const requires requires (Problem p) { { p.provides_get_D_N() } -> std::convertible_to<bool>; } { return problem.provides_get_D_N(); }
    [[nodiscard]] bool provides_eval_add_Q_N() const requires requires (Problem p) { { p.provides_eval_add_Q_N() } -> std::convertible_to<bool>; } { return problem.provides_eval_add_Q_N(); }
    [[nodiscard]] bool provides_eval_add_R_prod_masked() const requires requires (Problem p) { { p.provides_eval_add_R_prod_masked() } -> std::convertible_to<bool>; } { return problem.provides_eval_add_R_prod_masked(); }
    [[nodiscard]] bool provides_eval_add_S_prod_masked() const requires requires (Problem p) { { p.provides_eval_add_S_prod_masked() } -> std::convertible_to<bool>; } { return problem.provides_eval_add_S_prod_masked(); }
    [[nodiscard]] bool provides_get_R_work_size() const requires requires (Problem p) { { p.provides_get_R_work_size() } -> std::convertible_to<bool>; } { return problem.provides_get_R_work_size(); }
    [[nodiscard]] bool provides_get_S_work_size() const requires requires (Problem p) { { p.provides_get_S_work_size() } -> std::convertible_to<bool>; } { return problem.provides_get_S_work_size(); }
    [[nodiscard]] bool provides_eval_constr() const requires requires (Problem p) { { p.provides_eval_constr() } -> std::convertible_to<bool>; } { return problem.provides_eval_constr(); }
    [[nodiscard]] bool provides_eval_constr_N() const requires requires (Problem p) { { p.provides_eval_constr_N() } -> std::convertible_to<bool>; } { return problem.provides_eval_constr_N(); }
    [[nodiscard]] bool provides_eval_jac_constr() const requires requires (Problem p) { { p.provides_eval_jac_constr() } -> std::convertible_to<bool>; } { return problem.provides_eval_jac_constr(); }
    [[nodiscard]] bool provides_eval_jac_constr_N() const requires requires (Problem p) { { p.provides_eval_jac_constr_N() } -> std::convertible_to<bool>; } { return problem.provides_eval_jac_constr_N(); }
    [[nodiscard]] bool provides_eval_grad_constr_prod() const requires requires (Problem p) { { p.provides_eval_grad_constr_prod() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_constr_prod(); }
    [[nodiscard]] bool provides_eval_grad_constr_prod_N() const requires requires (Problem p) { { p.provides_eval_grad_constr_prod_N() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_constr_prod_N(); }
    [[nodiscard]] bool provides_eval_add_gn_hess_constr() const requires requires (Problem p) { { p.provides_eval_add_gn_hess_constr() } -> std::convertible_to<bool>; } { return problem.provides_eval_add_gn_hess_constr(); }
    [[nodiscard]] bool provides_eval_add_gn_hess_constr_N() const requires requires (Problem p) { { p.provides_eval_add_gn_hess_constr_N() } -> std::convertible_to<bool>; } { return problem.provides_eval_add_gn_hess_constr_N(); }
    // clang-format on

    std::shared_ptr<OCPEvalCounter> evaluations = std::make_shared<OCPEvalCounter>();
    Problem problem;

    OCProblemWithCounters(const Problem &problem) : problem(problem) {}
    OCProblemWithCounters(Problem &&problem)
        requires(!std::is_lvalue_reference_v<Problem>)
        : problem(std::forward<Problem>(problem)) {}

  private:
    template <class TimeT, class FunT>
    static decltype(auto) timed(TimeT &time, FunT &&f) {
        alpaqa::detail::Timed timed{time};
        return std::forward<FunT>(f)();
    }
};

template <class Problem>
[[nodiscard]] auto ocproblem_with_counters(Problem &&p) {
    using Prob        = std::remove_cvref_t<Problem>;
    using ProbWithCnt = OCProblemWithCounters<Prob>;
    return ProbWithCnt{std::forward<Problem>(p)};
}

template <class Problem>
[[nodiscard]] auto ocproblem_with_counters_ref(Problem &p) {
    using Prob        = std::remove_cvref_t<Problem>;
    using ProbWithCnt = OCProblemWithCounters<const Prob &>;
    return ProbWithCnt{p};
}

} // namespace alpaqa