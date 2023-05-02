#include "casadi/casadi.hpp"
#include <iostream>

using namespace casadi;

int main(int argc, char *argv[])
{
  SX x = SX::sym("x");
  SX y = SX::sym("y");
  
  std::cout << x+y << std::endl;
}
