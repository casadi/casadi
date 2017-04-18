// Sparse matrix-vector multiplication: z <- z + x*y
template<typename real_t>
void CASADI_PREFIX(mv)(const real_t* x, const int* sp_x, const real_t* y, real_t* z, int tr);
