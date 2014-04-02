#include "symbolic/std_vector_tools.hpp"
#include "symbolic/sx/sx_element.hpp"

using namespace CasADi;
using namespace std;

main(int argc, char *argv[])
{
  SXElement x("x");
  SXElement y("y");
  
  cout << x+y << endl;
}
