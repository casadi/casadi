#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdio.h>

int helloworld_c_init(int *f_type, int *n_in, int *n_out, int *sz_arg, int* sz_res) {
  *f_type = 0;
  *n_in = 0;
  *n_out = 1;
  return 0;
}

void helloworld_c(const double* arg, double* res) {
  printf("hello, C world!\n");

  double v[10];
  size_t i;
  for (i=0; i<10; ++i) {
    v[i] = i;
  }
  for (i=0; i<10; ++i) {
    v[i] = sin(v[i])+cos(v[i]);
  }
  res[0] = 0;
  for (i=0; i<10; ++i) {
    res[0] += v[i];
  }

  printf("returning %g\n", res[0]);
}

#ifdef __cplusplus
} /* extern "C" */
#endif
