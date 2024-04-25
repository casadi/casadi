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

// C-REPLACE "OracleMemory* m;" ""


// SYMBOL "oracle_data"
template<typename T1>
struct casadi_oracle_data {
  const T1** arg;
  T1** res;
  casadi_int* iw;
  T1* w;

  void* m;
};

// C-REPLACE "casadi_oracle_data<T1>" "struct casadi_oracle_data"

// SYMBOL "oracle_init"
template<typename T1>
void casadi_oracle_init(casadi_oracle_data<T1>* d, const T1*** arg, T1*** res,
    casadi_int** iw, T1** w) {
  d->arg = *arg;
  d->res = *res;
  d->iw = *iw;
  d->w = *w;
}
