#ifndef EXAMPLE_HPP
#define EXAMPLE_HPP

#include <time.h>
double My_variable = 3.0;

int fact(int n);
int my_mod(int x, int y);
char *get_time();

class MyClass{
public:
  MyClass(int n){
    n_fact_ = fact(n);
  }

  int get_n_fact(){
    return n_fact_;
  }
  
private:
  int n_fact_;
};

 

#endif // EXAMPLE_HPP
