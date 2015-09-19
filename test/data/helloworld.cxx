#include <cmath>
#include <iostream>
#include <vector>

extern "C"
int helloworld_cxx_init(int *f_type, int *n_in, int *n_out, int *sz_arg, int* sz_res) {
  *f_type = 0;
  *n_in = 0;
  *n_out = 1;
  return 0;
}

extern "C"
void helloworld_cxx(const double* arg, double* res) {
  std::cout << "hello, C++ world!" << std::endl;

  std::vector<double> v;
  for (size_t i=0; i<10; ++i) {
    v.push_back(i);
  }
  for (size_t i=0; i<10; ++i) {
    v[i] = sin(v[i])+cos(v[i]);
  }
  res[0] = 0;
  for (size_t i=0; i<10; ++i) {
    res[0] += v[i];
  }

  std::cout << "returning " << res[0] << std::endl;
}

