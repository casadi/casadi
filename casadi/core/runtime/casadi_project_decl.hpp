// Sparse copy: y <- x, w work vector (length >= number of rows)
template<typename real_t>
void CASADI_PREFIX(project)(const real_t* x, const int* sp_x, real_t* y, const int* sp_y, real_t* w);
