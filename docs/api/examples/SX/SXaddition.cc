#include "casadi/symbolic/std_vector_tools.hpp"
#include "casadi/symbolic/sx/sx_tools.hpp"

using namespace casadi;
using namespace std;

main(int argc, char *argv[])
{
  SX x = SX::sym("x");
  SX y = SX::sym("y");
  
  cout << x+y << endl;
}
