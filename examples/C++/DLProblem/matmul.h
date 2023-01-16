#pragma once

#include <alpaqa/interop/dl/dl-problem.h>

// A: m×n, B: n×p, C: m×p, column major
inline static void matmul(length_t m, length_t n, length_t p, const real_t *A,
                          const real_t *B, real_t *C) {
    for (index_t l = 0; l < m * p; ++l)
        C[l] = 0;
    for (index_t k = 0; k < p; ++k)
        for (index_t j = 0; j < n; ++j)
            for (index_t i = 0; i < m; ++i)
                C[i + k * m] += A[i + j * m] * B[j + k * n];
}

// A: m×n, b: m×1, c: n×1, column major
inline static void matvec_transp(length_t m, length_t n, const real_t *A,
                                 const real_t *b, real_t *c) {
    for (index_t l = 0; l < n; ++l)
        c[l] = 0;
    for (index_t j = 0; j < n; ++j)
        for (index_t i = 0; i < m; ++i)
            c[j] += A[i + j * m] * b[i];
}
