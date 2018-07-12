#include "../globals/globals.h"
#include "panoc.h"
#include "function_evaluator.h"
#include "./matrix_operations.h"

#include"math.h"
#include<stdlib.h>
#include<stdio.h>

#include "../include/optimizer.h"

static real_t compute_constraint_part_LA(const real_t* x,real_t* output_gradient);
static int compute_lagrangian_multipliers(real_t* new_x);

static real_t* lagrangian_multipliers;
static real_t* slacked_constraint_violations; /* The violation of the constraint with it's slack variables */
static real_t* constraint_evaluations;
static real_t* buffer_part_constraint;

static real_t (*cost_gradient_function)(const real_t* input,real_t* output_gradient); /* old cost function/gradient */
static real_t cost_gradient_function_LA(const real_t* input,real_t* output_gradient); /* new cost function/gradient */

static struct optimizer_extended_problem* extended_problem;
static real_t constraint_weight;

int optimizer_init_extended_box(struct optimizer_extended_problem* extended_problem_,real_t lower_bound,real_t upper_bound,real_t constraint_weight_){
    extended_problem = extended_problem_;
    constraint_weight = constraint_weight_;
    cost_gradient_function = extended_problem->problem.cost_gradient_function;

    if(optimizer_init_with_box(&extended_problem->problem,lower_bound ,upper_bound)==FAILURE)
        goto fail_1;

    /* If the box constraint PANOC optimizer was succesfull created */
    cost_gradient_function  = extended_problem->problem.cost_gradient_function; /* save the cost/gradient function */
    extended_problem->problem.cost_gradient_function=cost_gradient_function_LA; /* set the cost/gradient to the Lagragian cost/gradient*/

    lagrangian_multipliers = malloc(sizeof(real_t)*extended_problem->problem.dimension);
    if(lagrangian_multipliers==NULL)goto fail_2;

    constraint_evaluations = malloc(sizeof(real_t)*extended_problem->number_of_constraints);
    if(constraint_evaluations==NULL) goto fail_3;

    slacked_constraint_violations = malloc(sizeof(real_t)*extended_problem->number_of_constraints);
    if(slacked_constraint_violations==NULL) goto fail_4;

    buffer_part_constraint = malloc(sizeof(real_t)*extended_problem->problem.dimension);
    if(buffer_part_constraint==NULL) goto fail_5;
    
    return SUCCESS;

    fail_5:
        free(slacked_constraint_violations);
    fail_4:
        free(constraint_evaluations);
    fail_3:
        free(lagrangian_multipliers);
    fail_2:
        optimizer_cleanup();
    fail_1:
        return FAILURE;
}

int optimizer_cleanup_extended(void){
    free(buffer_part_constraint);
    free(slacked_constraint_violations);
    free(constraint_evaluations);
    free(lagrangian_multipliers);
    optimizer_cleanup();
}

/*
 * Replacement function for the cost and gradient function
 */
static real_t cost_gradient_function_LA(const real_t* input,real_t* output_gradient){
    /* fill up output_gradient with the gradient of f */
    real_t buffer =  cost_gradient_function(input,output_gradient);

    /* compute constraint related stuff for this position, do this first output gradient is also used as buffer */
    buffer += compute_constraint_part_LA(input,output_gradient);
   
    return buffer;
}

/*
 * This function is called by the user externally to get the solution
 */
int solve_extended_problem(real_t* solution){
    unsigned i;

    /* init the slack variables to zero */
    for (i = 0; i < extended_problem->number_of_constraints; i++){
        lagrangian_multipliers[i]=0;
    }

    unsigned int interations_till_convergence = 0;
    for( i = 0; i < extended_problem->max_loops; i++){
        
        /* solve the problem with the current slack variables */
        interations_till_convergence += solve_problem(solution);

        compute_lagrangian_multipliers(solution);
    }

    return interations_till_convergence;
}

/*
 * Calculate reusable buffer buffer_part_constraint
 *      buffer_part_constraint = h(x)+(1/c)*lagrangian_multiplier - z_c(x,lagrangian_multiplier)
 */
static int compute_buffer_part_constraint(void){
    unsigned int index_constraint;real_t LA=0;
    for (index_constraint = 0; index_constraint < extended_problem->number_of_constraints; index_constraint++){
        buffer_part_constraint[index_constraint]= constraint_evaluations[index_constraint] + \
             (1/constraint_weight)*lagrangian_multipliers[index_constraint]- \
             slacked_constraint_violations[index_constraint];        
    }

    return SUCCESS;
}

static int compute_slacked_constraint_violations(void){
    /* calculate the violation */
    unsigned int index_constraint;
    for (index_constraint = 0; index_constraint < extended_problem->number_of_constraints; index_constraint++){
        slacked_constraint_violations[index_constraint] = MIN(extended_problem->upper_bounds_constraints[index_constraint],\
            MAX(constraint_evaluations[index_constraint],extended_problem->lower_bounds_constraints[index_constraint])\
        );
    }
    return SUCCESS;
}
/* LA = c/2 * || buffer_part_constraint ||^2 - 1/(2*c) || Lagrangian_multipliers||^2 */
static real_t compute_LA_cost(const real_t* output_gradient){
    return (constraint_weight/2)*inner_product(buffer_part_constraint,buffer_part_constraint,extended_problem->number_of_constraints)\
        - (1/(2*constraint_weight))*inner_product(lagrangian_multipliers,lagrangian_multipliers,extended_problem->number_of_constraints);
}
static int add_LA_gradient(real_t* output_gradient){
    extended_problem->constraints_forwad_diff(output_gradient,buffer_part_constraint,buffer_part_constraint);
    vector_add_ntimes(output_gradient,buffer_part_constraint,extended_problem->number_of_constraints,constraint_weight);

    return SUCCESS;
}

/*
 * Compute the constraint part of the lagrangian (cost and gradient)
 */
static real_t compute_constraint_part_LA(const real_t* x,real_t* output_gradient){
    /* calculate intermediate variables */
    extended_problem->constraints(x,constraint_evaluations);
    compute_slacked_constraint_violations();
    compute_buffer_part_constraint();
    
    /* calculate cost and gradient */
    real_t LA = compute_LA_cost(output_gradient);
    add_LA_gradient(output_gradient);

    return LA;
}

/*
 * Compute the lagrangian multipliers with the formula
 *     y = y + constraint_weight*(g(new_x)-z_c(new_x,y))
 */
static int compute_lagrangian_multipliers(real_t* new_x){
    extended_problem->constraints(new_x,constraint_evaluations);
    compute_slacked_constraint_violations();

    real_t* buffer = buffer_part_constraint; /* reuse this buffer for something else */
    vector_add_ntimes(constraint_evaluations,slacked_constraint_violations,extended_problem->number_of_constraints,-1.); /* g(x_c) - z_c(x_c,lagrangian_multipliers) */
    vector_add_ntimes(lagrangian_multipliers,constraint_evaluations,extended_problem->number_of_constraints,constraint_weight);

    return SUCCESS;
}