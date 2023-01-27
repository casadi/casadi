#include "casadi/casadi.hpp"

using namespace casadi;

int main(int argc, char *argv[])
{
  SX x = SX::sym("x");
  SX y = SX::sym("y");
  
  cout << x+y << endl;
}
