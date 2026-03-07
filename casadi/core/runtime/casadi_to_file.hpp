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

// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"
// C-REPLACE "static_cast<int>" "(int) "

// SYMBOL "to_file"
template<typename T1>
void casadi_to_file(FILE* f, const casadi_int* sp, const T1* x) {
  casadi_int nrow = sp[0];
  casadi_int ncol = sp[1];
  const casadi_int* colind = sp+2;
  const casadi_int* row = colind + ncol + 1;
  casadi_int nnz = colind[ncol];
  casadi_int cc, k;
  fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
  fprintf(f, "%d %d %d\n", static_cast<int>(nrow), static_cast<int>(ncol),
    static_cast<int>(nnz));
  for (cc=0; cc<ncol; ++cc) {
    for (k=colind[cc]; k<colind[cc+1]; ++k) {
      fprintf(f, "%d %d ", static_cast<int>(row[k]+1), static_cast<int>(cc+1));
      if (x) {
        T1 v = x[k];
        if (v!=v) {
          fprintf(f, "nan");
        } else if (v==std::numeric_limits<T1>::infinity()) {
          fprintf(f, "inf");
        } else if (v==-std::numeric_limits<T1>::infinity()) {
          fprintf(f, "-inf");
        } else {
          fprintf(f, "%.16e", v);
        }
      } else {
        fprintf(f, "nan");
      }
      fprintf(f, "\n");
    }
  }
}
