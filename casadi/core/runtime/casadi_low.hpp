// NOLINT(legal/copyright)
// SYMBOL "low"
template<typename T1>
int casadi_low(T1 x, const double* grid, int ng, int lookup_mode) {
  if (lookup_mode) {
    double g0 = grid[0];
    int ret = (int) ((x-g0)*(ng-1)/(grid[ng-1]-g0)); // NOLINT(readability/casting)
    if (ret<0) ret=0;
    if (ret>ng-2) ret=ng-2;
    return ret;
  } else {
    int i;
    for (i=0; i<ng-2; ++i) {
      if (x < grid[i+1]) break;
    }
    return i;
  }
}
