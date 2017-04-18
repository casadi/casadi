// TRANS: y <- trans(x) , w work vector (length >= rows x)
template<typename real_t>
void CASADI_PREFIX(trans)(const real_t* x, const int* sp_x, real_t* y, const int* sp_y, int *tmp);
