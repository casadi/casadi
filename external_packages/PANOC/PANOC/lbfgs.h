#ifndef LBFGS_H
#define LBFGS_H

#include "../globals/globals.h"
#include <stddef.h>
#include <stdlib.h>

int lbfgs_init(const unsigned int buffer_size_,const unsigned int dimension_problem);
int lbfgs_cleanup(void);
int lbfgs_update_hessian(real_t tau, const real_t* current_location, const real_t* new_location);

/*
 * returns the direction calculated with lbfgs
 */ 
const real_t* lbfgs_get_direction(void);
int lbfgs_reset_iteration_counters(void);
unsigned int lbfgs_get_active_buffer_size(void);

#endif 