#include "proximal_gradient_descent.h"
#include "function_evaluator.h"
#include "../globals/globals.h"
#include <stddef.h>
#include <stdlib.h>
#include "matrix_operations.h"
#include "lipschitz.h"
#include"buffer.h"
#include"lbfgs.h"

/* functions only used internally */
static int proximal_gradient_descent_check_linesearch(void);
static int proximal_gradient_descent_forward_backward_step(const real_t* location,const real_t* df_location);
static int proximal_gradient_descent_push(void);
static int proximal_gradient_descent_linesearch(void);
static real_t proximal_gradient_descent_forward_backward_envelop_precomputed_step(const real_t f_location,const real_t* df_location);

/* variables saved between direction calls */
static size_t iteration_index=0; /* is 0 at first and 1 after first time you call get direction */
static real_t linesearch_gamma=0; /* linesearch parameter */
static real_t last_current_residual_inf_norm=0;
static unsigned int dimension;

/* variables used by each iteration */
static real_t* new_location;
static real_t* direction;
static real_t* new_location_FBE;
static real_t* direction_FBE;
static real_t g_new_location;


int proximal_gradient_descent_init(unsigned int dimension_problem){
	dimension = dimension_problem;

    new_location = malloc(sizeof(real_t)*dimension);
    if(new_location==NULL)goto fail_1;

    direction = malloc(sizeof(real_t)*dimension);
    if(direction==NULL)goto fail_2;

    new_location_FBE = malloc(sizeof(real_t)*dimension);
    if(new_location_FBE==NULL)goto fail_3;

    direction_FBE = malloc(sizeof(real_t)*dimension);
    if(direction_FBE==NULL)goto fail_4;

    if(init_lipschitz_estimator(dimension_problem)==FAILURE) goto fail_5;

    return SUCCESS;

    /*
     * something went wrong when allocating memory, free up what was already taken
     */
    fail_5:
        free(direction_FBE);
    fail_4:
        free(new_location_FBE);
    fail_3:
        free(direction);
    fail_2:
        free(new_location);
    fail_1:
        return FAILURE;
}
int proximal_gradient_descent_cleanup(void){
    free(new_location);
    free(direction);
    free(new_location_FBE);
    free(direction_FBE);
    cleanup_lipschtiz_estimator();
    proximal_gradient_descent_reset_iteration_counters();
    return SUCCESS;
}

/*
 * Reset iteration index and gamma, call me if your starting with a new problem
 */
int proximal_gradient_descent_reset_iteration_counters(void){
    iteration_index=0;
    linesearch_gamma=0;
    return SUCCESS;
}
/*
 * Find the proximal gradient descent with linesearch
 */
const real_t* proximal_gradient_descent_get_direction(void){
   /* 
    * If this is the first time you call me, find the initial gamma value
    * by estimating the Lipschitz value of df
    */
    if(iteration_index==0){
        real_t lipschitz_value = get_lipschitz();

        linesearch_gamma = (1-PROXIMAL_GRAD_DESC_SAFETY_VALUE)/lipschitz_value;
        iteration_index++; /* index only needs to increment if it is 0 */
    }
    proximal_gradient_descent_linesearch();
    return direction;
}

const real_t* proximal_gradient_descent_get_buffered_direction(void){
    return direction;
}

/*
 * This function performs the linesearch
 */
static int proximal_gradient_descent_linesearch(void){
    const real_t* current_location = buffer_get_current_location();
    const real_t* df_current_location = buffer_get_current_df();

    proximal_gradient_descent_forward_backward_step(current_location,df_current_location);
    while(proximal_gradient_descent_check_linesearch()==FAILURE){
        linesearch_gamma=linesearch_gamma/2;
        proximal_gradient_descent_forward_backward_step(current_location,df_current_location);
        lbfgs_reset_iteration_counters();
    }
    return SUCCESS;
}

/*
 * check if the linesearch condition is satisfied
 */
static int proximal_gradient_descent_check_linesearch(void){
    const real_t* df_current_location=buffer_get_current_df();
    const real_t inner_product_df_direction = inner_product(df_current_location,direction,dimension);

    const real_t f_current_location=buffer_get_current_f();
    buffer_evaluate_pure_prox_location(new_location);
    const real_t f_new_location=buffer_get_pure_prox_f();

    const real_t norm_direction_gamma = inner_product(direction,direction,dimension); /* direction=gamma*r in paper */

    if(f_new_location>f_current_location - inner_product_df_direction + ( (1-PROXIMAL_GRAD_DESC_SAFETY_VALUE)/(2*linesearch_gamma) )*norm_direction_gamma + 1e-6*f_current_location)
        return FAILURE;
    else
        return SUCCESS;
}


/* 
 * This function performs an forward backward step. x=prox(x-gamma*df(x))
 */
static int proximal_gradient_descent_forward_backward_step(const real_t* location,const real_t* df_location){
    vector_copy(location,new_location,dimension);
    vector_add_ntimes(new_location,df_location,dimension,-1*linesearch_gamma); /* new_location = location - gamma * df_location */
    
    g_new_location = function_evaluator_proxg(new_location,linesearch_gamma); /* new_location = proxg(new_location) */
    
    vector_copy(location,direction,dimension); /* find the direction */
    vector_add_ntimes(direction,new_location,dimension,-1.);
     
    return SUCCESS;
}

/*
 * returns the residual, R(x) = 1/gamma[ x- proxg(x-df(x)*gamma)]
 */
int proximal_gradient_descent_get_new_residual(const real_t* location,real_t* residual){
    proximal_gradient_descent_push(); /* swap data fields to FBE calculation fields */

    /* 
     * use special buffer so values might be reused later
     */
    const real_t* df_location=buffer_get_pure_prox_df();

    proximal_gradient_descent_forward_backward_step(location,df_location);
    /* calculate the residual, as in normalize the current direction */
    vector_copy(direction,residual,dimension);
    vector_real_mul(residual,dimension,1/linesearch_gamma);

    proximal_gradient_descent_push(); /* swap data fields to FBE calculation fields */
    return SUCCESS;
}
/*
 * compute residual using direction used by fordward backward envelop
 */
int proximal_gradient_descent_get_new_residual_buffered(real_t* residual){
    proximal_gradient_descent_push(); /* swap data fields to FBE calculation fields */

    /* calculate the residual, as in normalize the current direction */
    vector_copy(direction,residual,dimension);
    vector_real_mul(residual,dimension,1/linesearch_gamma);


    proximal_gradient_descent_push(); /* swap data fields to FBE calculation fields */
    return SUCCESS;
}

/*
 * returns the residual using the previous forward backward step, R(x) = 1/gamma[ x- proxg(x-df(x)*gamma)]
 */
int proximal_gradient_descent_get_current_residual(real_t* residual){
    /* calculate the residual, as in normalize the current direction */
    vector_copy(direction,residual,dimension);
    vector_real_mul(residual,dimension,1/linesearch_gamma);


    /* calculate the inf-norm and safe it */
    last_current_residual_inf_norm=vector_norm_inf(residual,dimension);

    return SUCCESS;
}


real_t proximal_gradient_descent_get_current_residual_inf_norm(void){
    return last_current_residual_inf_norm;
}


/*
 * calculate the forward backward envelop using the internal gamma
 * Matlab cache.FBE = cache.fx + cache.gz - cache.gradfx(:)'*cache.FPR(:) + (0.5/gam)*(cache.normFPR^2);
 */
real_t proximal_gradient_descent_forward_backward_envelop(const real_t* location){
    proximal_gradient_descent_push(); /* swap data fields to FBE calculation fields */

    /* 
     * use special buffer so values might be reused later
     */
    buffer_evaluate_new_location(location);
    const real_t* df_location=buffer_get_new_location_df();
    real_t f_location=buffer_get_new_location_f();

    proximal_gradient_descent_forward_backward_step(location, df_location); /* this will fill the new_direction variable */
    real_t forward_backward_envelop = proximal_gradient_descent_forward_backward_envelop_precomputed_step(f_location,df_location);

    proximal_gradient_descent_push(); /* swap data fields to FBE calculation fields */

    return forward_backward_envelop;
}
real_t proximal_gradient_descent_get_current_forward_backward_envelop(void){
    const real_t f_location = buffer_get_current_f();
    const real_t* df_location = buffer_get_current_df();
    return proximal_gradient_descent_forward_backward_envelop_precomputed_step(f_location,df_location);
}

/*
 * calculate the forward backward envelop using internal forwardbackward step
 */
static real_t proximal_gradient_descent_forward_backward_envelop_precomputed_step(const real_t f_location,const real_t* df_location){
    const real_t norm_direction = inner_product(direction,direction,dimension);

    const real_t forward_backward_envelop = f_location + g_new_location \
     - inner_product(df_location,direction,dimension) \
     + (1/(linesearch_gamma*2))*norm_direction;

    return forward_backward_envelop;
}

real_t proximal_gradient_descent_get_gamma(void){
    return linesearch_gamma;
}

/*
 * Replace the new_location and direction arrays with 2 arrays used to compute the FBE
 */
static int proximal_gradient_descent_push(void){
    real_t* buffer;
    
    buffer=direction;
    direction=direction_FBE;
    direction_FBE=buffer;

    buffer=new_location;
    new_location=new_location_FBE;
    new_location_FBE=buffer;
    return SUCCESS;
}