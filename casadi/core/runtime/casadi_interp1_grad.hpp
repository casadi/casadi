// NOLINT(legal/copyright)
// SYMBOL "interp1_grad"
template<typename T1>
T1 casadi_interp1_grad(const T1* grid, const casadi_int* offset, const T1* values, T1 x, casadi_int lookup_mode) { // NOLINT(whitespace/line_length)
  casadi_int index;

  // Left index and fraction of interval
  index = casadi_low(x, grid, offset[1], lookup_mode);

  return (values[index+1]-values[index])/(grid[index+1]-grid[index]);
}