// NOLINT(legal/copyright)
// SYMBOL "de_boor"
template<typename T1>
void casadi_de_boor(T1 x, const T1* knots, int n_knots, int degree, T1* boor) {
  // length boor: n_knots-1
  for (int d=1;d<degree+1;++d) {
    for (int i=0;i<n_knots-d-1;++i) {
      T1 b = 0;
      T1 bottom = knots[i + d] - knots[i];
      if (bottom) b = (x - knots[i]) * boor[i] / bottom;
      bottom = knots[i + d + 1] - knots[i + 1];
      if (bottom) b += (knots[i + d + 1] - x) * boor[i + 1] / bottom;
      boor[i] = b;
    }
  }
}
