/* File : example.c */
#include "test_swig_matlab.hpp"
#include <time.h>
     
int fact(int n) {
  if (n <= 1) return 1;
  else return n*fact(n-1);
}
     
int my_mod(int x, int y) {
  return (x%y);
}

char *get_time()
{
  time_t ltime;
  time(&ltime);
  return ctime(&ltime);
}
     

