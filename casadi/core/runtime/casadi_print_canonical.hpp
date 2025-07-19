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

// SYMBOL "print_canonical"
template<typename T1>
void casadi_print_canonical(const casadi_int* sp, const T1* x) {
  // C-REPLACE "printf" "CASADI_PRINTF"
  // C-REPLACE "static_cast<int>" "(int) "
  if (x) {
    casadi_int nrow = sp[0];
    casadi_int ncol = sp[1];
    const casadi_int* colind = sp+2;
    const casadi_int* row = colind + ncol + 1;
    casadi_int nnz = sp[2+ncol];
    if (nrow==1 && ncol==1 && nnz==1) {
      casadi_print_scalar(x[0]);
    } else {
      casadi_int i;
      printf("%dx%d: ", static_cast<int>(nrow), static_cast<int>(ncol));
      casadi_print_vector(nnz, x);
      if (nnz!=nrow*ncol) {
        printf(", colind: [");
        for (i = 0; i < ncol; ++i) {
          if (i > 0) printf(", ");
          printf("%d", static_cast<int>(colind[i]));
        }
        printf("], row: [");
        for (i = 0; i < nnz; ++i) {
          if (i > 0) printf(", ");
          printf("%d", static_cast<int>(row[i]));
        }
        printf("]");
      }
    }
  } else {
    printf("NULL");
  }

}
