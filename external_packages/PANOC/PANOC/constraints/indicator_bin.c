/*
 * Indicator binary function:
 *      - map any value x on the binary set x_0 and x_1
 */

#include "../../globals/globals.h"
#include "indicator_bin.h"
#include <stdlib.h>

struct indicator_bin* init_indicator_bin(\
    const unsigned int dimension,const real_t lower_bound,const real_t upper_bound){

    struct indicator_bin* data = malloc(sizeof(struct indicator_bin));
    if(data==NULL) return NULL;

    data->dimension=dimension;
    data->lower_bound=lower_bound;
    data->upper_bound=upper_bound;

    return data;
}

real_t prox_indicator_bin(const struct indicator_bin* data,real_t* input, real_t gamma){
    unsigned int i;
    real_t g=0; /* the new pont x is by definition allway's inside the box */

    real_t average = (data->upper_bound+data->lower_bound)/2;
    for (i = 0; i < data->dimension; i++)
    {
        if(input[i]<average)
            input[i]=data->lower_bound;
        else
            input[i]=data->upper_bound;
    }
    
    return g;
}