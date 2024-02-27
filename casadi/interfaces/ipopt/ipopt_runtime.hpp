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


// C-REPLACE "casadi_nlpsol_prob<T1>" "struct casadi_nlpsol_prob"
// C-REPLACE "casadi_nlpsol_data<T1>" "struct casadi_nlpsol_data"

// C-REPLACE "reinterpret_cast<int**>" "(int**) "
// C-REPLACE "reinterpret_cast<int*>" "(int*) "
// C-REPLACE "const_cast<int*>" "(int*) "

template<typename T1>
struct casadi_ipopt_prob {
  const casadi_nlpsol_prob<T1>* nlp;

  const int *colinda, *rowa;
  const int *colindh, *rowh;
  const int *integrality;
};
// C-REPLACE "casadi_ipopt_prob<T1>" "struct casadi_ipopt_prob"

// SYMBOL "ipopt_setup"
template<typename T1>
void casadi_ipopt_setup(casadi_ipopt_prob<T1>* p) {

}



// SYMBOL "ipopt_data"
template<typename T1>
struct casadi_ipopt_data {
  // Problem structure
  const casadi_ipopt_prob<T1>* prob;
  // Problem structure
  casadi_nlpsol_data<T1>* nlp;

  int return_status;

};
// C-REPLACE "casadi_ipopt_data<T1>" "struct casadi_ipopt_data"

// SYMBOL "ipopt_init_mem"
template<typename T1>
int ipopt_init_mem(casadi_ipopt_data<T1>* d) {
  //d->ipopt = Highs_create();
  return 0;
}

// SYMBOL "ipopt_free_mem"
template<typename T1>
void ipopt_free_mem(casadi_ipopt_data<T1>* d) {
  //Highs_destroy(d->ipopt);
}

// SYMBOL "ipopt_work"
template<typename T1>
void casadi_ipopt_work(const casadi_ipopt_prob<T1>* p, casadi_int* sz_arg, casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_nlpsol_work(p->nlp, sz_arg, sz_res, sz_iw, sz_w);
}

// SYMBOL "ipopt_init"
template<typename T1>
void casadi_ipopt_init(casadi_ipopt_data<T1>* d, const T1*** arg, T1*** res, casadi_int** iw, T1** w) {
  

}