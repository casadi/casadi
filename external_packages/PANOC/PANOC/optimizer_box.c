#include "../globals/globals.h"
#include "panoc.h"
#include "function_evaluator.h"

#include"math.h"
#include<stdlib.h>
#include<stdio.h>

#include "../include/optimizer.h"
#include "constraints/indicator_box.h"

static struct indicator_box* indicator_box_data; /* data related to the constraints used by the prox operator */

static real_t prox(real_t* input,real_t gamma){
    return prox_indicator_box(indicator_box_data,input,gamma);
}

int optimizer_init_with_box(struct optimizer_problem* problem,real_t lower_bound,real_t upper_bound){

    /* prepare proxg(x) */
    indicator_box_data = init_indicator_box(problem->dimension,lower_bound,upper_bound);
    if(indicator_box_data==NULL) goto fail_1;
    problem->proxg=prox;

    if(optimizer_init(problem)==FAILURE) goto fail_2;

    return SUCCESS;
     
    fail_2:
        free(indicator_box_data);
    fail_1:
        return FAILURE;
}