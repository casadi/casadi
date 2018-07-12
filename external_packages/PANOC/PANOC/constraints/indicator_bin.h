#ifndef INDICATOR_BIN_H
#define INDICATOR_BIN_H

#include "../../globals/globals.h"

struct indicator_bin{
    real_t lower_bound;
    real_t upper_bound;
    unsigned int dimension;
};

struct indicator_bin* init_indicator_bin(const unsigned int dimension,const real_t lower_bound,const real_t upper_bound);
real_t prox_indicator_bin(const struct indicator_bin* data,real_t* input, real_t gamma);


#endif