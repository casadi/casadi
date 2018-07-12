#include "../globals/globals.h"
#include <stddef.h>
#include <stdlib.h>

#ifndef LIPSCHITZ_H
#define LIPSCHITZ_H

int init_lipschitz_estimator(unsigned int dimension);
real_t get_lipschitz(void);
int cleanup_lipschtiz_estimator(void);

#endif