#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/box.hpp>
#include <alpaqa/problem/structure.hpp>
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
    static length_t default_get_R_work_size(const void *, const OCProblemVTable &) { return 0; }
    static length_t default_get_S_work_size(const void *, const OCProblemVTable &) { return 0; }
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
    length_t get_N() const { return vtable.N; }
    /// Number of inputs.
    length_t get_nu() const { return vtable.nu; }
    /// Number of states.
    length_t get_nx() const { return vtable.nx; }
    /// Number of outputs.
    length_t get_nh() const { return vtable.nh; }
    length_t get_nh_N() const { return vtable.nh_N; }
    /// Number of constraints.
    length_t get_nc() const { return vtable.nc; }
    length_t get_nc_N() const { return vtable.nc_N; }
    /// All dimensions
    Dim get_dim() const {
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
    void get_x_init(rvec x_init) const;
    void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const;
    void eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const;
    void eval_grad_f_prod(index_t timestep, crvec x, crvec u, crvec p, rvec grad_fxu_p) const;
    void eval_h(index_t timestep, crvec x, crvec u, rvec h) const;
    void eval_h_N(crvec x, rvec h) const;
    real_t eval_l(index_t timestep, crvec h) const;
    real_t eval_l_N(crvec h) const;
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
    length_t get_R_work_size() const;
    length_t get_S_work_size() const;
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

} // namespace alpaqa