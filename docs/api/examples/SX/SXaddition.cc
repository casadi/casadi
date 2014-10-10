#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/sx/sx_tools.hpp"

using namespace casadi;
using namespace std;

int main(int argc, char *argv[])
{
  SX x = SX::sym("x");
  SX y = SX::sym("y");
  
  cout << x+y << endl;
}
