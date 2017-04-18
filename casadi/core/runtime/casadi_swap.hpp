template<typename real_t>
void CASADI_PREFIX(swap)(int n, real_t* x, int inc_x, real_t* y, int inc_y) {
  real_t t;
  int i;
  for (i=0; i<n; ++i) {
    t = *x;
    *x = *y;
    *y = t;
    x += inc_x;
    y += inc_y;
  }
}
