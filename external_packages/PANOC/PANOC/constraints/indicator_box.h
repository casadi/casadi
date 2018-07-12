#ifndef INDICATOR_BOX_H
#define INDICATOR_BOX_H

#include "../../globals/globals.h"

struct indicator_box{
    unsigned int dimension;
    real_t lower_bound;
    real_t upper_bound;
    real_t inf;
};

real_t prox_indicator_box(const struct indicator_box* data, real_t* input, real_t gamma);

struct indicator_box* init_indicator_box(\
    const unsigned int dimension_,const real_t lower_bound,const real_t upper_bound);

#endif