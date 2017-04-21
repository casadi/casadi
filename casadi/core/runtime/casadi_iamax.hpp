// NOLINT(legal/copyright)
template<typename T1>
int CASADI_PREFIX(iamax)(int n, const T1* x, int inc_x) {
  T1 t;
  T1 largest_value = -1.0;
  int largest_index = -1;
  int i;
  for (i=0; i<n; ++i) {
    t = fabs(*x);
    x += inc_x;
    if (t>largest_value) {
      largest_value = t;
      largest_index = i;
    }
  }
  return largest_index;
}
