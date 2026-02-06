//
//    MIT No Attribution
//
//    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

// Tensor times vector in all modes (ttv) - recursive multi-mode contraction.
//
// Given an N-dimensional coefficient tensor C of shape (m, n_0, n_1, ..., n_{N-1})
// stored in column-major order with m as the leading (fastest) dimension,
// and N weight vectors w_0, w_1, ..., w_{N-1}, compute:
//
//   ret_j += sum_{i_0, ..., i_{N-1}}  w_0[i_0] * w_1[i_1] * ... * w_{N-1}[i_{N-1}]
//                                      * C[j, s_0+i_0, s_1+i_1, ..., s_{N-1}+i_{N-1}]
//
// for j = 0, ..., m-1, where s_k = starts[k] selects the active sub-tensor.
//
// The result is accumulated into ret (not cleared).
//
// The recursion peels off one dimension at a time (outermost first):
//   - dim > 0:  loop over i_{dim}, multiply weight, add offset, recurse
//   - dim == 0: inner-product with the m-wide coefficient row
//
// Parameters:
//   ret       [m]          output vector (accumulated into, not cleared)
//   dim                    current dimension index (call with n_dims-1)
//   n_dims                 total number of dimensions N
//   all_w     [sum n_k]    packed weight vectors [w_0 | w_1 | ... | w_{N-1}]
//   w_offset  [N+1]        w_k starts at all_w[w_offset[k]], length w_offset[k+1]-w_offset[k]
//   starts    [N]          start index per dimension (selects sub-tensor)
//   strides   [N]          stride per dimension in the flat coefficient array
//                          strides[0] = m, strides[k+1] = strides[k] * n_k
//   c         [total]      flat coefficient tensor
//   m                      number of outputs (leading dimension)
//   weight                 accumulated product of weights from outer dims (call with 1)
//   offset                 accumulated flat offset from outer dims (call with 0)
//
// SYMBOL "tensor_ttv"
template<typename T1>
void casadi_tensor_ttv(T1* ret, casadi_int dim, casadi_int n_dims,
    const T1* all_w, const casadi_int* w_offset,
    const casadi_int* starts, const casadi_int* strides,
    const T1* c, casadi_int m,
    T1 weight, casadi_int offset) {
  casadi_int i, j, n_w;
  const T1* w;
  // Number of weights and pointer for this dimension
  n_w = w_offset[dim+1] - w_offset[dim];
  w = all_w + w_offset[dim];
  if (dim == 0) {
    // Base case: innermost dimension - contract with coefficient row
    const T1* coeff = c + offset + starts[0]*m;
    for (i = 0; i < n_w; ++i) {
      T1 ww = weight * w[i];
      for (j = 0; j < m; ++j)
        ret[j] += ww * coeff[i*m + j];
    }
  } else {
    // Recursive case: peel off dimension `dim`, accumulate weight and offset
    for (i = 0; i < n_w; ++i) {
      casadi_tensor_ttv(ret, dim-1, n_dims, all_w, w_offset,
        starts, strides, c, m,
        weight * w[i],
        offset + (starts[dim]+i)*strides[dim]);
    }
  }
}
