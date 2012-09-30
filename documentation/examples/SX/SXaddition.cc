#include "symbolic/stl_vector_tools.hpp"
#include "symbolic/sx/sx.hpp"

using namespace CasADi;
using namespace std;

main(int argc, char *argv[])
{
  SX x("x");
  SX y("y");
  
  cout << x+y << endl;
}
