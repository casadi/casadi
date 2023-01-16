#include <alpaqa/interop/dl/dl-problem.h>
#include <problem-c-export.h>

#include <stdlib.h>
#include <string.h>

#include "matmul.h"

typedef alpaqa_real_t real_t;
typedef alpaqa_length_t length_t;

struct ProblemData {
    real_t *Q;
    real_t *A;
    real_t *work;

    alpaqa_problem_functions_t functions;
};

static void eval_grad_f(void *instance, const real_t *x, real_t *grad_f) {
    struct ProblemData *problem = instance;
    length_t n                  = problem->functions.n;
    matmul(n, n, 1, problem->Q, x, grad_f); // grad_f = Q x
}

static real_t eval_f_grad_f(void *instance, const real_t *x, real_t *grad_f) {
    struct ProblemData *problem = instance;
    length_t n                  = problem->functions.n;
    real_t result;
    eval_grad_f(instance, x, grad_f);
    matmul(1, n, 1, x, grad_f, &result); // result = xᵀ grad_f
    return (real_t)0.5 * result;
}

static real_t eval_f(void *instance, const real_t *x) {
    struct ProblemData *problem = instance;
    return eval_f_grad_f(instance, x, problem->work);
}

static void eval_g(void *instance, const real_t *x, real_t *gx) {
    struct ProblemData *problem = instance;
    length_t n                  = problem->functions.n;
    length_t m                  = problem->functions.m;
    matmul(m, n, 1, problem->A, x, gx); // gx = A x
}

static void eval_grad_g_prod(void *instance, const real_t *x, const real_t *y,
                             real_t *grad_gxy) {
    (void)x;
    struct ProblemData *problem = instance;
    length_t n                  = problem->functions.n;
    length_t m                  = problem->functions.m;
    matvec_transp(m, n, problem->A, y, grad_gxy); // grad_gxy = Aᵀ y
}

static struct ProblemData *create_problem(void *user_data) {
    (void)user_data;
    struct ProblemData *problem = malloc(sizeof(*problem));
    memset(problem, 0, sizeof(*problem));
    size_t n = 2, m = 1;
    problem->functions.n                = (length_t)n;
    problem->functions.m                = (length_t)m;
    problem->functions.eval_f           = &eval_f;
    problem->functions.eval_grad_f      = &eval_grad_f;
    problem->functions.eval_f_grad_f    = &eval_f_grad_f;
    problem->functions.eval_g           = &eval_g;
    problem->functions.eval_grad_g_prod = &eval_grad_g_prod;
    problem->Q                          = malloc(sizeof(real_t) * n * n);
    problem->A                          = malloc(sizeof(real_t) * m * n);
    problem->work                       = malloc(sizeof(real_t) * n);
    // Note: column major order
    problem->Q[0] = 3;
    problem->Q[1] = -1;
    problem->Q[2] = -1;
    problem->Q[3] = 3;
    problem->A[0] = 2;
    problem->A[1] = 1;
    return problem;
}

static void cleanup_problem(void *instance) {
    struct ProblemData *problem = instance;
    free(problem->Q);
    free(problem->A);
    free(problem->work);
    free(problem);
}

PROBLEM_C_EXPORT alpaqa_problem_register_t
alpaqa_problem_register(void *user_data) {
    struct ProblemData *problem = create_problem(user_data);
    alpaqa_problem_register_t result;
    memset(&result, 0, sizeof(result));
    result.instance  = problem;
    result.cleanup   = &cleanup_problem;
    result.functions = &problem->functions;
    return result;
}
