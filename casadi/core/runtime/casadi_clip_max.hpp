// NOLINT(legal/copyright)
// SYMBOL "clip_max"
template<typename T1>
void casadi_clip_max(T1* x, casadi_int n, T1 max, const casadi_int* mask) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i) {
      if (!mask || mask[i]) {
        if (x[i] > max) {
          x[i] = max;
        }
      }
    }
  }
}
