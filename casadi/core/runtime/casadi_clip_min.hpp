// NOLINT(legal/copyright)
// SYMBOL "clip_min"
template<typename T1>
void casadi_clip_min(T1* x, casadi_int n, T1 min, const casadi_int* mask) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i) {
      if (!mask || mask[i]) {
        if (x[i] < min) {
          x[i] = min;
        }
      }
    } 
  }
}
