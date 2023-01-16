#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef double real_t;
typedef long int length_t;
typedef length_t index_t;

typedef struct {
    length_t n, m;

    // clang-format off
    real_t (*eval_f)(
        void *instance,
        const real_t *x);
    void (*eval_grad_f)(
        void *instance,
        const real_t *x,
        real_t *grad_fx);
    void (*eval_g)(
        void *instance,
        const real_t *x,
        real_t *gx);
    void (*eval_grad_g_prod)(
        void *instance,
        const real_t *x,
        const real_t *y,
        real_t *grad_gxy);
    void (*eval_grad_gi)(
        void *instance,
        const real_t *x,
        index_t i,
        real_t *grad_gi);
    void (*eval_hess_L_prod)(
        void *instance,
        const real_t *x,
        const real_t *y,
        const real_t *v,
        real_t *Hv);
    void (*eval_hess_L)(
        void *instance,
        const real_t *x,
        const real_t *y,
        real_t *H);
    real_t (*eval_f_grad_f)(
        void *instance,
        const real_t *x,
        real_t *grad_fx);
    real_t (*eval_f_g)(
        void *instance,
        const real_t *x,
        real_t *g);
    real_t (*eval_f_grad_f_g)(
        void *instance,
        const real_t *x,
        real_t *grad_fx,
        real_t *g);
    void (*eval_grad_f_grad_g_prod)(
        void *instance,
        const real_t *x,
        const real_t *y,
        real_t *grad_f,
        real_t *grad_gxy);
    void (*eval_grad_L)(
        void *instance,
        const real_t *x,
        const real_t *y,
        real_t *grad_L,
        real_t *work_n);
    real_t (*eval_ψ)(
        void *instance,
        const real_t *x,
        const real_t *y,
        const real_t *Σ,
        real_t *ŷ);
    void (*eval_grad_ψ_from_ŷ)(
        void *instance,
        const real_t *x,
        const real_t *ŷ,
        real_t *grad_ψ,
        real_t *work_n);
    void (*eval_grad_ψ)(
        void *instance,
        const real_t *x,
        const real_t *y,
        const real_t *Σ,
        real_t *grad_ψ,
        real_t *work_n,
        real_t *work_m);
    real_t (*eval_ψ_grad_ψ)(
        void *instance,
        const real_t *x,
        const real_t *y,
        const real_t *Σ,
        real_t *grad_ψ,
        real_t *work_n,
        real_t *work_m);
    // clang-format on
} alpaqa_problem_functions_t;

typedef struct {
    /// Owning pointer.
    void *instance;
    /// Non-owning pointer, lifetime at least as long as @ref instance.
    alpaqa_problem_functions_t *functions;
    /// Pointer to the function to clean up @ref instance.
    void (*cleanup)(void *);
} alpaqa_problem_register_t;

#ifdef __cplusplus
}
#endif
