// NOLINT(legal/copyright)
template<typename T1>
T1 CASADI_PREFIX(polyval)(const T1* p, int n, T1 x) {
  T1 r=p[0];
  int i;
  for (i=1; i<=n; ++i) {
    r = r*x + p[i];
  }
  return r;
}
