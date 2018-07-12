#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "../globals/globals.h"
struct solver_parameters{
    unsigned int max_interations;
    real_t tolerance;
    unsigned int buffer_size;
};

struct optimizer_problem{
    unsigned int dimension;
    real_t (*proxg)(real_t* input, real_t gamma);
    real_t (*cost_gradient_function)(const real_t* input,real_t* output_gradient);
    struct solver_parameters solver_params;
};

struct optimizer_extended_problem{
    struct optimizer_problem problem;
    int (*constraints)(const real_t* x,real_t* out);
    int (*constraints_forwad_diff)(const real_t* x,const real_t* y,real_t* out); /* returns inner product -> out=constraint(x)^T*y */
    unsigned int number_of_constraints;
    real_t* lower_bounds_constraints;
    real_t* upper_bounds_constraints;
    unsigned int max_loops;
};

/*
 * Initializes a optimizer that makes use of the PANOC algorithm.
 */
int optimizer_init(struct optimizer_problem* problem_);

/*
 * This interfaces a optimization problem onto the panoc algorithm:
 *      min cost(x)
 *          subject to  L<x<H
 */
int optimizer_init_with_box(struct optimizer_problem* problem,real_t lower_bound,real_t upper_bound);

/*
 * This interfaces a optimization problem onto the panoc algorithm:
 *      min cost(x)
 *          subject to  g(x)=0
 */
int optimizer_init_with_costum_constraint(struct optimizer_problem* problem_,real_t (*proxg)(real_t* x, real_t gamma));

/*
 * This interfaces the format of the casadi library onto the panoc algorithm:
 *      min cost(x)
 *          subject to  hmin < h(x) < hmax
 *                      xmin <  x   < xmax
 *      weight : weight of the constraints h(x)
 */
int optimizer_init_extended_box(struct optimizer_extended_problem* extended_problem,real_t lower_bound,real_t upper_bound,real_t constraint_weight_);

int optimizer_cleanup(void);
int optimizer_cleanup_extended(void);

int solve_problem(real_t* solution);
int solve_extended_problem(real_t* solution);

#endif