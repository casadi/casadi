// NOLINT(legal/copyright)
// SYMBOL "swap"
template<typename T1>
void casadi_swap(int n, T1* x, int inc_x, T1* y, int inc_y) {
  T1 t;
  int i;
  for (i=0; i<n; ++i) {
    t = *x;
    *x = *y;
    *y = t;
    x += inc_x;
    y += inc_y;
  }
}
