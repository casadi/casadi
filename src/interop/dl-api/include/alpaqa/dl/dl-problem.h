#pragma once

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#define ALPAQA_DL_ABI_VERSION 0xA1A000000001

#ifdef __cplusplus
extern "C" {
#define ALPAQA_DL_ABI_VERSION_DEFAULT = ALPAQA_DL_ABI_VERSION
#else
#define ALPAQA_DL_ABI_VERSION_DEFAULT
#endif

typedef double alpaqa_real_t;
typedef ptrdiff_t alpaqa_length_t;
typedef alpaqa_length_t alpaqa_index_t;

typedef struct {
    uint64_t abi_version ALPAQA_DL_ABI_VERSION_DEFAULT;
    alpaqa_length_t n, m;

    // clang-format off
    alpaqa_real_t (*eval_f)(
        void *instance,
        const alpaqa_real_t *x);
    void (*eval_grad_f)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *grad_fx);
    void (*eval_g)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *gx);
    void (*eval_grad_g_prod)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        alpaqa_real_t *grad_gxy);
    void (*eval_jac_g)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_index_t *inner_idx,
        alpaqa_index_t *outer_ptr,
        alpaqa_real_t *J_values);
    alpaqa_length_t (*get_jac_g_num_nonzeros)(
        void *instance);
    void (*eval_grad_gi)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_index_t i,
        alpaqa_real_t *grad_gi);
    void (*eval_hess_L_prod)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        alpaqa_real_t scale,
        const alpaqa_real_t *v,
        alpaqa_real_t *Hv);
    void (*eval_hess_L)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        alpaqa_real_t scale,
        alpaqa_index_t *inner_idx,
        alpaqa_index_t *outer_ptr,
        alpaqa_real_t *H_values);
    alpaqa_length_t (*get_hess_L_num_nonzeros)(
        void *instance);
    void (*eval_hess_ψ_prod)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        const alpaqa_real_t *Σ,
        alpaqa_real_t scale,
        const alpaqa_real_t *zl,
        const alpaqa_real_t *zu,
        const alpaqa_real_t *v,
        alpaqa_real_t *Hv);
    void (*eval_hess_ψ)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        const alpaqa_real_t *Σ,
        alpaqa_real_t scale,
        const alpaqa_real_t *zl,
        const alpaqa_real_t *zu,
        alpaqa_index_t *inner_idx,
        alpaqa_index_t *outer_ptr,
        alpaqa_real_t *H_values);
    alpaqa_length_t (*get_hess_ψ_num_nonzeros)(
        void *instance);
    alpaqa_real_t (*eval_f_grad_f)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *grad_fx);
    alpaqa_real_t (*eval_f_g)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *g);
    void (*eval_grad_f_grad_g_prod)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        alpaqa_real_t *grad_f,
        alpaqa_real_t *grad_gxy);
    void (*eval_grad_L)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        alpaqa_real_t *grad_L,
        alpaqa_real_t *work_n);
    alpaqa_real_t (*eval_ψ)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        const alpaqa_real_t *Σ,
        const alpaqa_real_t *zl,
        const alpaqa_real_t *zu,
        alpaqa_real_t *ŷ);
    void (*eval_grad_ψ)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        const alpaqa_real_t *Σ,
        const alpaqa_real_t *zl,
        const alpaqa_real_t *zu,
        alpaqa_real_t *grad_ψ,
        alpaqa_real_t *work_n,
        alpaqa_real_t *work_m);
    alpaqa_real_t (*eval_ψ_grad_ψ)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        const alpaqa_real_t *Σ,
        const alpaqa_real_t *zl,
        const alpaqa_real_t *zu,
        alpaqa_real_t *grad_ψ,
        alpaqa_real_t *work_n,
        alpaqa_real_t *work_m);
    alpaqa_real_t (*eval_prox_grad_step)(
        void *instance,
        alpaqa_real_t γ,
        const alpaqa_real_t *x,
        const alpaqa_real_t *grad_ψ,
        alpaqa_real_t *x̂,
        alpaqa_real_t *p);
    void (*initialize_box_C)(
        void *instance,
        alpaqa_real_t *lb,
        alpaqa_real_t *ub);
    void (*initialize_box_D)(
        void *instance,
        alpaqa_real_t *lb,
        alpaqa_real_t *ub);
    void (*initialize_l1_reg)(
        void *instance,
        alpaqa_real_t *lambda,
        alpaqa_length_t *size);
    // clang-format on
} alpaqa_problem_functions_t;

/// Opaque type for a C++-only map of extra functions.
typedef struct alpaqa_function_dict_s alpaqa_function_dict_t;

typedef struct {
    /// Owning pointer.
    void *instance;
    /// Non-owning pointer, lifetime at least as long as @ref instance.
    alpaqa_problem_functions_t *functions;
    /// Pointer to the function to clean up @ref instance.
    void (*cleanup)(void *);
    /// Pointer to a map of extra functions (C++ only).
    alpaqa_function_dict_t *extra_functions;
} alpaqa_problem_register_t;

typedef struct {
    uint64_t abi_version ALPAQA_DL_ABI_VERSION_DEFAULT;
    alpaqa_length_t N, nx, nu, nh, nh_N, nc, nc_N;

    // clang-format off
    void (*get_U)(
        void *instance,
        alpaqa_real_t *lb,
        alpaqa_real_t *ub);
    void (*get_D)(
        void *instance,
        alpaqa_real_t *lb,
        alpaqa_real_t *ub);
    void (*get_D_N)(
        void *instance,
        alpaqa_real_t *lb,
        alpaqa_real_t *ub);
    void (*get_x_init)(
        void *instance,
        alpaqa_real_t *x_init);
    void (*eval_f)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *x,
        const alpaqa_real_t *u,
        alpaqa_real_t *fxu);
    void (*eval_jac_f)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *x,
        const alpaqa_real_t *u,
        alpaqa_real_t *J_fxu);
    void (*eval_grad_f_prod)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *x,
        const alpaqa_real_t *u,
        const alpaqa_real_t *p,
        alpaqa_real_t *grad_fxu_p);
    void (*eval_h)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *x,
        const alpaqa_real_t *u,
        alpaqa_real_t *h);
    void (*eval_h_N)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *h);
    alpaqa_real_t (*eval_l)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *h);
    alpaqa_real_t (*eval_l_N)(
        void *instance,
        const alpaqa_real_t *h);
    void (*eval_qr)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *xu,
        const alpaqa_real_t *h,
        alpaqa_real_t *qr);
    void (*eval_q_N)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *h,
        alpaqa_real_t *q);
    void (*eval_add_Q)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *xu,
        const alpaqa_real_t *h,
        alpaqa_real_t *Q);
    void (*eval_add_Q_N)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *h,
        alpaqa_real_t *Q);
    void (*eval_add_R_masked)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *xu,
        const alpaqa_real_t *h,
        const alpaqa_index_t *mask,
        alpaqa_real_t *R,
        alpaqa_real_t *work);
    void (*eval_add_S_masked)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *xu,
        const alpaqa_real_t *h,
        const alpaqa_index_t *mask,
        alpaqa_real_t *S,
        alpaqa_real_t *work);
    void (*eval_add_R_prod_masked)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *xu,
        const alpaqa_real_t *h,
        const alpaqa_index_t *mask_J,
        const alpaqa_index_t *mask_K,
        const alpaqa_real_t *v,
        alpaqa_real_t *out,
        alpaqa_real_t *work);
    void (*eval_add_S_prod_masked)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *xu,
        const alpaqa_real_t *h,
        const alpaqa_index_t *mask_K,
        const alpaqa_real_t *v,
        alpaqa_real_t *out,
        alpaqa_real_t *work);
    alpaqa_length_t (*get_R_work_size)(
        void *instance);
    alpaqa_length_t (*get_S_work_size)(
        void *instance);
    void (*eval_constr)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *x,
        alpaqa_real_t *c);
    void (*eval_constr_N)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *c);
    void (*eval_grad_constr_prod)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *x,
        const alpaqa_real_t *p,
        alpaqa_real_t *grad_cx_p);
    void (*eval_grad_constr_prod_N)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *p,
        alpaqa_real_t *grad_cx_p);
    void (*eval_add_gn_hess_constr)(
        void *instance,
        alpaqa_index_t timestep,
        const alpaqa_real_t *x,
        const alpaqa_real_t *M,
        alpaqa_real_t *out);
    void (*eval_add_gn_hess_constr_N)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *M,
        alpaqa_real_t *out);
    // clang-format on
} alpaqa_control_problem_functions_t;

typedef struct {
    /// Owning pointer.
    void *instance;
    /// Non-owning pointer, lifetime at least as long as @ref instance.
    alpaqa_control_problem_functions_t *functions;
    /// Pointer to the function to clean up @ref instance.
    void (*cleanup)(void *);
    /// Pointer to a map of extra functions (C++ only).
    alpaqa_function_dict_t *extra_functions;
} alpaqa_control_problem_register_t;

#ifdef __cplusplus
}
#endif

#ifndef __cplusplus
#define ALPAQA_PROBLEM_FUNCTIONS_INIT(self)                                    \
    do {                                                                       \
        memset((self), 0, sizeof(*(self)));                                    \
        (self)->abi_version = ALPAQA_DL_ABI_VERSION;                           \
    } while (0)
#endif

#if defined(__cplusplus) && __cplusplus > 201703L

#include <any>
#include <functional>
#include <map>
#include <string>

struct alpaqa_function_dict_s {
    std::map<std::string, std::any> dict{};
};

namespace alpaqa {

using function_dict_t             = alpaqa_function_dict_t;
using problem_register_t          = alpaqa_problem_register_t;
using control_problem_register_t  = alpaqa_control_problem_register_t;
using problem_functions_t         = alpaqa_problem_functions_t;
using control_problem_functions_t = alpaqa_control_problem_functions_t;

/// Make the given function available to alpaqa.
/// @see @ref alpaqa::dl::DLProblem::call_extra_func
/// @see @ref alpaqa::dl::DLControlProblem::call_extra_func
template <class Func>
void register_function(function_dict_t *&extra_functions, std::string name,
                       Func &&func) {
    if (extra_functions == nullptr)
        extra_functions = new function_dict_t{};
    extra_functions->dict.insert_or_assign(
        std::move(name), std::function{std::forward<Func>(func)});
}

template <class Func>
void register_function(problem_register_t &result, std::string name,
                       Func &&func) {
    register_function(result.extra_functions, std::move(name),
                      std::forward<Func>(func));
}

template <class Func>
void register_function(control_problem_register_t &result, std::string name,
                       Func &&func) {
    register_function(result.extra_functions, std::move(name),
                      std::forward<Func>(func));
}

template <class Result, class T, class Ret, class... Args>
void register_member_function(Result &result, std::string name,
                              Ret (T::*member)(Args...)) {
    register_function(result, std::move(name),
                      [member](void *self_, Args... args) -> Ret {
                          auto *self = reinterpret_cast<T *>(self_);
                          return (self->*member)(std::forward<Args>(args)...);
                      });
}

template <class Result, class T, class Ret, class... Args>
void register_member_function(Result &result, std::string name,
                              Ret (T::*member)(Args...) const) {
    register_function(result, std::move(name),
                      [member](const void *self_, Args... args) -> Ret {
                          const auto *self = reinterpret_cast<const T *>(self_);
                          return (self->*member)(std::forward<Args>(args)...);
                      });
}

namespace detail {
template <auto Member, class Class, class Ret, class... Args>
static auto member_caller(Ret (Class::*)(Args...)) {
    return [](void *self_, Args... args) -> Ret {
        auto *self = reinterpret_cast<Class *>(self_);
        return (self->*Member)(std::forward<Args>(args)...);
    };
}

template <auto Member, class Class, class Ret, class... Args>
static auto member_caller(Ret (Class::*)(Args...) const) {
    return []<class Self>(Self * self_, Args... args) -> Ret
               requires std::is_void_v<Self>
    {
        const auto *self = reinterpret_cast<const Class *>(self_);
        return (self->*Member)(std::forward<Args>(args)...);
    };
}

template <auto Member, class Class, class Ret>
static auto member_caller(Ret Class::*) {
    return []<class Self>(Self * self_) -> decltype(auto)
               requires std::is_void_v<Self>
    {
        using CClass = std::conditional_t<std::is_const_v<Self>,
                                          std::add_const_t<Class>, Class>;
        auto *self   = reinterpret_cast<CClass *>(self_);
        return self->*Member;
    };
}
} // namespace detail

/// Wrap the given member function of signature into a lambda function that
/// accepts the instance as a void pointer.
///
/// - `Ret Class::member(args...)` → `Ret(void *self, args...)`
/// - `Ret Class::member(args...) const` → `Ret(const void *self, args...)`
/// - `Type Class::member` → `Type &(void *self)`
/// - `Type Class::member` → `const Type &(const void *self)`
template <auto Member>
static auto member_caller() {
    return detail::member_caller<Member>(Member);
}

/// Cleans up the extra functions registered by @ref register_function.
/// @note   This does not need to be called for the functions returned by the
///         registration function, those functions will be cleaned up by alpaqa.
/// @note   The @ref alpaqa_problem_register_t and
///         @ref alpaqa_control_problem_register_t structs are part of the C API
///         and do not automatically clean up their resources when destroyed,
///         you have to do it manually by calling this function.
inline void unregister_functions(function_dict_t *&extra_functions) {
    delete extra_functions;
}

} // namespace alpaqa

#endif

#undef ALPAQA_DL_ABI_VERSION_DEFAULT
