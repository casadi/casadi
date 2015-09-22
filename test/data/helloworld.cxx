#include <cmath>
#include <vector>
#ifdef WITH_IOSTREAM
#include <iostream>
#else
#include <stdio.h>
#endif

extern "C"
int helloworld_cxx_init(int *f_type, int *n_in, int *n_out, int *sz_arg, int* sz_res) {
  *f_type = 0;
  *n_in = 0;
  *n_out = 1;
  return 0;
}

extern "C"
void helloworld_cxx(const double* arg, double* res) {
#ifdef WITH_IOSTREAM
  std::cout << "hello, C++ world!" << std::endl;
#else
  printf("hello, C++ world!\n");
#endif

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

#ifdef WITH_IOSTREAM
  std::cout << "returning " << res[0] << std::endl;
#else
  printf("returning %f\n", res[0]);
#endif
}
