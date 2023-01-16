#pragma once

#include <alpaqa/interop/dl/dl-problem.h>

// A: m×n, B: n×p, C: m×p, column major
inline static void matmul(alpaqa_length_t m, alpaqa_length_t n,
                          alpaqa_length_t p, const alpaqa_real_t *A,
                          const alpaqa_real_t *B, alpaqa_real_t *C) {
    for (alpaqa_index_t l = 0; l < m * p; ++l)
        C[l] = 0;
    for (alpaqa_index_t k = 0; k < p; ++k)
        for (alpaqa_index_t j = 0; j < n; ++j)
            for (alpaqa_index_t i = 0; i < m; ++i)
                C[i + k * m] += A[i + j * m] * B[j + k * n];
}

// A: m×n, b: m×1, c: n×1, column major
inline static void matvec_transp(alpaqa_length_t m, alpaqa_length_t n,
                                 const alpaqa_real_t *A, const alpaqa_real_t *b,
                                 alpaqa_real_t *c) {
    for (alpaqa_index_t l = 0; l < n; ++l)
        c[l] = 0;
    for (alpaqa_index_t j = 0; j < n; ++j)
        for (alpaqa_index_t i = 0; i < m; ++i)
            c[j] += A[i + j * m] * b[i];
}
