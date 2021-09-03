// NOLINT(legal/copyright)
// SYMBOL "interp1_weights"
template<typename T1>
void casadi_interp1_weights(const T1* grid, const casadi_int* offset, T1 x, T1* alpha, casadi_int* index, casadi_int lookup_mode) { // NOLINT(whitespace/line_length)
  // Left index and fraction of interval
    casadi_int j;
    // Find left index
    j = *index = casadi_low(x, grid, offset[1], lookup_mode);
    // Get interpolation/extrapolation alpha
    *alpha = (x-grid[j])/(grid[j+1]-grid[j]);
}
