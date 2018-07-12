#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include "../globals/globals.h"
#include "stddef.h"
#include "math.h"
 
/* squared with preprocessor */
#define sq(x) ((x)*(x))

/* copy vector1 into vector2 */
void vector_copy(const real_t* vector1,real_t* vector2,const int size_vector);

/* multiply vector with real */
void vector_real_mul(real_t* vector,const int size_vector,const real_t real);

/*
 * calculate the 2 norm of a vector defined as
 * sqrt(x[0]^2 + x[1]^2 + ... x[n]^2)
 */
real_t vector_norm2(const real_t* vector,const int vector_size);

real_t vector_norm1(const real_t* vector,const int vector_size);

real_t vector_norm_inf(const real_t* vector,const int vector_size);

/* inner product between vector 1 and 2 */
real_t inner_product(const real_t* vector1,const real_t* vector2,const int size_vector);

/* add vector2 n times to vector1 */
void vector_add_ntimes(real_t* vector1,const real_t* vector2,const int size_vector,const real_t n);

/* add vector2 a_vector2 times to vector1 and add vector3 a_vector3 times to vector1*/
void vector_add_2_vectors_a_times(const real_t* vector1,const real_t* vector2,const real_t* vector3,const int size_vector,
    const real_t a_vector2,const real_t a_vector3,real_t* result);

#endif