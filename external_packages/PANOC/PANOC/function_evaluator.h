#ifndef function_evaluator_H
#define function_evaluator_H

#include "../globals/globals.h"
#include "../include/optimizer.h"

#include <stddef.h>
#include <stdlib.h>

int function_evaluator_init(const struct optimizer_problem* problem);
int function_evaluator_cleanup(void);

/* 
 * returns the cost over the stack
 * returns the gradient in real_t* output
 */
real_t function_evaluator_f_df(const real_t* input,real_t* output);
real_t function_evaluator_proxg(real_t* input, real_t gamma);

#endif