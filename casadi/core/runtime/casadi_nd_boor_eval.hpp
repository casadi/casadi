// NOLINT(legal/copyright)
// SYMBOL "nd_boor_eval"
template<typename T1>
void casadi_nd_boor_eval(T1* ret, int n_dims, const T1* all_knots, const int* offset, const int* all_degree, const int* strides, const T1* c, int m, const T1* all_x, const int* lookup_mode, int reverse, int* iw, T1* w) { // NOLINT(whitespace/line_length)
  int* boor_offset = iw; iw+=n_dims+1;
  int* starts = iw; iw+=n_dims;
  int* index = iw; iw+=n_dims;
  int* coeff_offset = iw;

  T1* cumprod = w; w+= n_dims+1;
  T1* all_boor = w;

  boor_offset[0] = 0;
  cumprod[n_dims] = 1;
  coeff_offset[n_dims] = 0;

  int n_iter = 1;
  for (int k=0;k<n_dims;++k) {
    T1* boor = all_boor+boor_offset[k];

    int degree = all_degree[k];
    const T1* knots = all_knots + offset[k];
    int n_knots = offset[k+1]-offset[k];
    int n_b = n_knots-degree-1;

    T1 x = all_x[k];
    int L = casadi_low(x, knots+degree, n_knots-2*degree, lookup_mode[k]);

    int start = L;
    if (start>n_b-degree-1) start = n_b-degree-1;

    starts[k] = start;

    casadi_fill(boor, 2*degree+1, 0.0);
    if (x>=knots[0] && x<=knots[n_knots-1]) {
      if (x==knots[1]) {
        casadi_fill(boor, degree+1, 1.0);
      } else if (x==knots[n_knots-1]) {
        boor[degree] = 1;
      } else if (knots[L+degree]==x) {
        boor[degree-1] = 1;
      } else {
        boor[degree] = 1;
      }
    }
    casadi_de_boor(x, knots+start, 2*degree+2, degree, boor);
    boor+= degree+1;
    n_iter*= degree+1;
    boor_offset[k+1] = boor_offset[k] + degree+1;
  }

  casadi_fill_int(index, n_dims, 0);

  // Prepare cumulative product
  for (int pivot=n_dims-1;pivot>=0;--pivot) {
    cumprod[pivot] = (*(all_boor+boor_offset[pivot]))*cumprod[pivot+1];
    coeff_offset[pivot] = starts[pivot]*strides[pivot]+coeff_offset[pivot+1];
  }

  for (int k=0;k<n_iter;++k) {

    // accumulate result
    for (int i=0;i<m;++i) {
      if (reverse) {
        ret[coeff_offset[0]+i] += c[i]*cumprod[0];
      } else {
        ret[i] += c[coeff_offset[0]+i]*cumprod[0];
      }
    }

    // Increment index
    index[0]++;
    int pivot = 0;

    // Handle index overflow
    {
      // increment next index (forward)
      while (index[pivot]==boor_offset[pivot+1]-boor_offset[pivot]) {
        index[pivot] = 0;
        if (pivot==n_dims-1) break;
        index[++pivot]++;
      }

      // update cumulative structures (reverse)
      while (pivot>0) {
        // Compute product
        cumprod[pivot] = (*(all_boor+boor_offset[pivot]+index[pivot]))*cumprod[pivot+1];
        // Compute offset
        coeff_offset[pivot] = (starts[pivot]+index[pivot])*strides[pivot]+coeff_offset[pivot+1];
        pivot--;
      }
    }

    // Compute product
    cumprod[0] = (*(all_boor+index[0]))*cumprod[1];

    // Compute offset
    coeff_offset[0] = (starts[0]+index[0])*m+coeff_offset[1];

  }
}
