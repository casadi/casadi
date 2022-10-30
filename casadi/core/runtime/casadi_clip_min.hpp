// NOLINT(legal/copyright)
// SYMBOL "fill"
template<typename T1>
void casadi_clip_min(T1* x, T1* binary_vector, casadi_int n, T1 min) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i){
      if (binary_vector[i] != 0){
        if (x[i] < min){
          x[i] = min;
        }
          // x++;
        // } else{
        //   x++;
        // }
      }
    } 
  }
}
