template<typename real_t>
void CASADI_PREFIX(de_boor)(real_t x, const real_t* knots, int n_knots, int degree, real_t* boor) {
  /* length boor: n_knots-1 */
  for (int d=1;d<degree+1;++d) {
    for (int i=0;i<n_knots-d-1;++i) {
      real_t b = 0;
      real_t bottom = knots[i + d] - knots[i];
      if (bottom) b = (x - knots[i]) * boor[i] / bottom;
      bottom = knots[i + d + 1] - knots[i + 1];
      if (bottom) b += (knots[i + d + 1] - x) * boor[i + 1] / bottom;
      boor[i] = b;
    }
  }
}
