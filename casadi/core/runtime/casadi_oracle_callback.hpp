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

// SYMBOL "oracle_callback"
template<typename T1>
struct casadi_oracle_callback {
  int (*eval)(const T1** arg, T1** res, casadi_int* iw, T1* w, int mem);
  int (*checkout)(void);
  void (*release)(int);
};
// C-REPLACE "casadi_oracle_callback<T1>" "struct casadi_oracle_callback"
// C-REPLACE "casadi_oracle_data<T1>" "struct casadi_oracle_data"

// SYMBOL "call"
template<typename T1>
int casadi_oracle_call(const casadi_oracle_callback<T1>* cb, casadi_oracle_data<T1>* d) {
  int flag;
  int mem = 0;
  if (cb->checkout) mem = cb->checkout();
  flag = cb->eval(d->arg, d->res, d->iw, d->w, mem);
  if (cb->release) cb->release(mem);
  return flag;
}
