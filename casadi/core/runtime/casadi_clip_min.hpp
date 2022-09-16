// NOLINT(legal/copyright)
// SYMBOL "fill"
template<typename T1>
void casadi_clip_min(T1* x, casadi_int n, T1 min) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i){
      if (*x < min){
        *x++ = min;
      } else{
        x++;
      }
    } 
  }
}
