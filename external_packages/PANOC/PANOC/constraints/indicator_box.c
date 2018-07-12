#include "../../globals/globals.h"
#include "indicator_box.h"
#include <stdlib.h>

struct indicator_box* init_indicator_box(\
    const unsigned int dimension_,const real_t lower_bound,const real_t upper_bound){

    struct indicator_box* data = malloc(sizeof(struct indicator_box));
    if(data==NULL) return NULL;

    data->dimension=dimension_;
    data->lower_bound=lower_bound;
    data->upper_bound=upper_bound;

    return data;
}

real_t prox_indicator_box(const struct indicator_box* data,real_t* input, real_t gamma){
    unsigned int i;
    real_t g=0; /* the new pont x is by definition allway's inside the box */

    for (i = 0; i < data->dimension; i++)
    {
        if(input[i]<data->lower_bound)
        {
            input[i]=data->lower_bound;
        }
        else if (input[i]>data->upper_bound)
        {
            input[i]=data->upper_bound;
        }
    }
    
    return g;
}