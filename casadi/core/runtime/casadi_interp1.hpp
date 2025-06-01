// NOLINT(legal/copyright)
// SYMBOL "interp1"
template<typename T1>
T1 casadi_interp1(const T1* grid, const casadi_int* offset, const T1* values, T1 x, casadi_int lookup_mode) { // NOLINT(whitespace/line_length)
  // Work vectors
  T1 alpha;
  casadi_int index;
  // Left index and fraction of interval
  index = casadi_low(x, grid, offset[1], lookup_mode);
  alpha = (x-grid[index])/(grid[index+1]-grid[index]);
  return (1-alpha)*values[index]+alpha*values[index+1];
}