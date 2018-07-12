#include"lipschitz.h"
#include"matrix_operations.h"
#include"math.h"
#include"function_evaluator.h"
#include <stddef.h>
#include <stdlib.h>
#include "buffer.h"

unsigned int dimension;
real_t* df_current_position_delta;
real_t* current_position_delta;

static real_t get_lipschitz_get_delta(const real_t* current_position,real_t * current_position_delta);

/*
 * Get the step used to estimate the lipschitz constant
 * delta= max{small number,10^{-6}*u_0}
 */
static real_t get_lipschitz_get_delta(const real_t* current_position,real_t * current_position_delta){
    real_t norm_delta=0;

    size_t i;
    for (i = 0; i < dimension; i++)
    {
        if(DELTA_LIPSCHITZ_SAFETY_VALUE>current_position[i] ){
            current_position_delta[i]  = current_position[i] + DELTA_LIPSCHITZ;
            norm_delta += sq(DELTA_LIPSCHITZ);
        }else{
            current_position_delta[i] = current_position[i] + current_position[i]*DELTA_LIPSCHITZ_SAFETY_VALUE;
            norm_delta += sq(current_position[i]*DELTA_LIPSCHITZ_SAFETY_VALUE);
        }
    }

    return sqrt(norm_delta);
}
/*
 * Estimate the lipschitz constant by using the numerical hessian as an estimation
 * Theorem:
 *     ||gradient(x)|| < B
 *     f is B-lipschitz
 */
real_t get_lipschitz(void){
    /* get the curernt position an its gradient */
    const real_t* current_position = buffer_get_current_location();
    const real_t* df_current_position = buffer_get_current_df();

    /* get the shifted position  */
    const real_t denominator = get_lipschitz_get_delta(current_position,current_position_delta);
    
    /* get shifted gradient */
    function_evaluator_f_df(current_position_delta,df_current_position_delta);

    /* 
     * L = norm((df(x+delta)-df(x))/delta) 
     * reuse the current_position_delta values as buffer
     */
    vector_copy(df_current_position,current_position_delta,dimension); /* step1: df(x+delta)-df(x) */
    vector_add_ntimes(current_position_delta,df_current_position_delta,dimension,-1.); 
    
    const real_t numerator = vector_norm2(current_position_delta, dimension); /* step2: norm((df(x+delta)-df(x))) */

    return numerator/denominator;
}

int init_lipschitz_estimator(unsigned int dimension_problem){
    current_position_delta = malloc(dimension_problem*sizeof(real_t));
    if(current_position_delta==NULL) goto fail_1;

    df_current_position_delta = malloc(dimension_problem*sizeof(real_t));
    if(df_current_position_delta==NULL) goto fail_1;

    dimension = dimension_problem;

    return SUCCESS;

        free(current_position_delta);
    fail_1:
        return FAILURE;
}

int cleanup_lipschtiz_estimator(void){
    free(current_position_delta);
    free(df_current_position_delta);
    return SUCCESS;
}