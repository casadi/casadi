#ifndef BUFFER_H
#define BUFFER_H

#include "../globals/globals.h"
#include <stddef.h>
#include <stdlib.h>

int buffer_init(unsigned int dimension_problem);
int buffer_cleanup(void);

int buffer_renew(const real_t* current_location);

real_t buffer_get_current_f(void);
const real_t*  buffer_get_current_df(void);
const real_t* buffer_get_current_location(void);
int buffer_reset_cycle(void);

/*
 * 4 extra functions, these reuse df(x) and f(x) if it was a pure lbfgs step
 */
int buffer_evaluate_new_location(const real_t* lbfgs_new_location);
real_t buffer_get_new_location_f(void);
const real_t* buffer_get_new_location_df(void);
int buffer_set_new_location_as_current_location(void);


int buffer_evaluate_pure_prox_location(const real_t* pure_prox_location);
int buffer_set_pure_prox_location_as_current_location(void);
const real_t* buffer_get_pure_prox_df(void);
real_t buffer_get_pure_prox_f(void);

#endif