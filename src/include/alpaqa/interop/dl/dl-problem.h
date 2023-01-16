#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef double alpaqa_real_t;
typedef long int alpaqa_length_t;
typedef alpaqa_length_t alpaqa_index_t;

typedef struct {
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
    void (*eval_grad_gi)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_index_t i,
        alpaqa_real_t *grad_gi);
    void (*eval_hess_L_prod)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        const alpaqa_real_t *v,
        alpaqa_real_t *Hv);
    void (*eval_hess_L)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        alpaqa_real_t *H);
    alpaqa_real_t (*eval_f_grad_f)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *grad_fx);
    alpaqa_real_t (*eval_f_g)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *g);
    alpaqa_real_t (*eval_f_grad_f_g)(
        void *instance,
        const alpaqa_real_t *x,
        alpaqa_real_t *grad_fx,
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
        alpaqa_real_t *ŷ);
    void (*eval_grad_ψ_from_ŷ)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *ŷ,
        alpaqa_real_t *grad_ψ,
        alpaqa_real_t *work_n);
    void (*eval_grad_ψ)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        const alpaqa_real_t *Σ,
        alpaqa_real_t *grad_ψ,
        alpaqa_real_t *work_n,
        alpaqa_real_t *work_m);
    alpaqa_real_t (*eval_ψ_grad_ψ)(
        void *instance,
        const alpaqa_real_t *x,
        const alpaqa_real_t *y,
        const alpaqa_real_t *Σ,
        alpaqa_real_t *grad_ψ,
        alpaqa_real_t *work_n,
        alpaqa_real_t *work_m);
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

#ifdef __cplusplus
}
#endif

#if __cplusplus > 201703L

#include <any>
#include <functional>
#include <map>
#include <string>

struct alpaqa_function_dict_s {
    std::map<std::string, std::any> dict{};
};

namespace alpaqa {

using function_dict_t     = alpaqa_function_dict_t;
using problem_register_t  = alpaqa_problem_register_t;
using problem_functions_t = alpaqa_problem_functions_t;

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

} // namespace alpaqa

#endif
