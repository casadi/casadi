// NOLINT(legal/copyright)
// SYMBOL "iamax"
template<typename T1>
casadi_int casadi_iamax(casadi_int n, const T1* x, casadi_int inc_x) {
  T1 t;
  T1 largest_value = -1.0;
  casadi_int largest_index = -1;
  casadi_int i;
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
