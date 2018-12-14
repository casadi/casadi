// NOLINT(legal/copyright)

// C-REPLACE "fmax" "casadi_fmax"
// C-REPLACE "fmin" "casadi_fmin"
// C-REPLACE "std::isinf" "casadi_isinf"

// SYMBOL "bound_consistency"
template<typename T1>
void casadi_bound_consistency(casadi_int n, T1* x, T1* lam,
                                 const T1* lbx, const T1* ubx) {
    // Local variables
    casadi_int i;
    T1 lb, ub;
    // Loop over variables
    for (i=0; i<n; ++i) {
      // Get bounds
      lb = lbx ? lbx[i] : 0.;
      ub = ubx ? ubx[i] : 0.;
      // Make sure bounds are respected
      x[i] = fmin(fmax(x[i], lb), ub);
      // Adjust multipliers
      if (std::isinf(lb) && std::isinf(ub)) {
        // Both multipliers are infinite
        lam[i] = 0.;
      } else if (std::isinf(lb) || x[i] - lb > ub - x[i]) {
        // Infinite lower bound or closer to upper bound than lower bound
        lam[i] = fmax(0., lam[i]);
      } else if (std::isinf(ub) || x[i] - lb < ub - x[i]) {
        // Infinite upper bound or closer to lower bound than upper bound
        lam[i] = fmin(0., lam[i]);
      }
    }
}
