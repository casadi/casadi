// NOLINT(legal/copyright)
// SYMBOL "flip"
inline
int casadi_flip(int* corner, int ndim) {
  int i;
  for (i=0; i<ndim; ++i) {
    if (corner[i]) {
      corner[i]=0;
    } else {
      corner[i]=1;
      return 1;
    }
  }
  return 0;
}
