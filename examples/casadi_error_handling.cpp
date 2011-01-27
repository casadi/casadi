// Uncomment this line to compile without error handling
//#define CASADI_NDEBUG

#include <casadi/casadi_exception.hpp>

bool bad_test(){
  return false;
}

bool bad_test2(){
  // This will fail
  casadi_assert(bad_test());
  
  // Returns true, but the code won't reach this place
  return true;
}

int main(){
  
  // Warning
  casadi_warning("This function will fail.");
  
  // No warning here
  casadi_assert_warning(0==0, "Not here.");
  
  // Warning due to failed assert
  casadi_assert_warning(1==0, "I am telling you, it WILL fail.");
  
  // Recursive error
  casadi_assert(bad_test2());

  return 0;
}
