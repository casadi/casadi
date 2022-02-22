// NOLINT(legal/copyright)
// SYMBOL "logsumexp"
template<typename T1>
T1 casadi_logsumexp(const T1* x, casadi_int n) {
  casadi_int i, max_ind;
  T1 max, r;
  if (n==1) return x[0];
  max_ind = 0;
  max = x[0];
  // Determine max, argmax
  for (i=1; i<n; ++i) {
    if (x[i]>x[0]) {
      max = x[i];
      max_ind = i;
    }
  }
  r = 0;
  for (i=0; i<n; ++i) {
    if (i!=max_ind) r += exp(x[i]-max);
  }
  // C-REPLACE "log1p" "casadi_log1p"
  return log1p(r)+max;
}
