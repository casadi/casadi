// Sparse matrix-matrix multiplication: z <- z + x*y
template<typename real_t>
void CASADI_PREFIX(mtimes)(const real_t* x, const int* sp_x, const real_t* y, const int* sp_y, real_t* z, const int* sp_z, real_t* w, int tr);
