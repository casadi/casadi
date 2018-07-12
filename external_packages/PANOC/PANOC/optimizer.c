#include "../globals/globals.h"
#include "panoc.h"
#include "function_evaluator.h"

#include"math.h"
#include<stdlib.h>
#include<stdio.h>

#include "../include/optimizer.h"


static struct optimizer_problem* problem; /* data related to the problem */
static unsigned char initialized=FALSE;
static real_t* new_solution;

static void save_solution(real_t* solution);

int solve_problem(real_t* solution){
    /* 
     * take implicitly the previous inputs as the starting position for the algorithm 
     */
    unsigned int i_panoc;
    real_t current_residual=problem->solver_params.tolerance*10;
    for (i_panoc= 0; i_panoc < problem->solver_params.max_interations ; i_panoc++){
        if(current_residual<problem->solver_params.tolerance){
            /* if the residual is low enough stop iterating */
            break;
        }
        current_residual = panoc_get_new_location(solution,new_solution);

        /*
         * if the residual was larger then the machine accuracy
         * -> set the new_input as input for the next iteration 
         * WARNING: if the residual was smaller then the machine 
         *  accuracy you might get NaN thats why we won't use it
         */
        if(current_residual>MACHINE_ACCURACY)
            save_solution(solution);
            
    }

    panoc_reset_cycli(); /* reset all counters/buffers */

    return i_panoc;
}

static void save_solution(real_t* solution){
    unsigned int i;
    for (i = 0; i < problem->dimension; i++)
    {
        solution[i]=new_solution[i];
    }
}

int optimizer_init(struct optimizer_problem* problem_){
    if(initialized==TRUE) /* already initialized, remove old stuff first */
        optimizer_cleanup();
    
    problem = problem_;

    new_solution = malloc(sizeof(real_t)*problem->dimension);
    if(new_solution==NULL) goto fail_1;
    
    /* Prepare the solver */
    if(panoc_init(problem->dimension,problem->solver_params.buffer_size)==FAILURE) goto fail_2;

    /* Prepare the cost function */
    if(function_evaluator_init(problem)==FAILURE) goto fail_3;

    initialized=TRUE;
    return SUCCESS;

    fail_3:
        panoc_cleanup();
    fail_2:
        free(new_solution);
    fail_1:
        return FAILURE;
}

int optimizer_init_with_costum_constraint(struct optimizer_problem* problem_,real_t (*proxg)(real_t* x, real_t gamma)){

    /* prepare proxg(x) */
    problem_->proxg = proxg;

    if(optimizer_init(problem_)==FAILURE) return FAILURE;

    return SUCCESS;
}

int optimizer_cleanup(void){
    panoc_cleanup();
    function_evaluator_cleanup();
    free(new_solution);
    initialized=FALSE;
    return SUCCESS;
}