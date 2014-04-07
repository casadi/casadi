#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/std_vector_tools.hpp"
#include <limits>

using namespace CasADi;

int main(int argc, char *argv[])
{
  DMatrix d(3,1,0);

  assert(!all(d));
  assert(!any(d));
  std::cout << d << "   all: " << all(d) << "  any: " <<  any(d) << std::endl;
  
  
  d(0) = 0;
  d(1) = 0;
  d(2) = 1;
  
  assert(!all(d));
  assert(any(d));
  std::cout << d << "   all: " << all(d) << "  any: " <<  any(d) << std::endl;
  
  d(0) = 1;
  d(1) = 1;
  d(2) = 1;
  
  assert(all(d));
  assert(any(d));
  std::cout << d << "   all: " << all(d) << "  any: " <<  any(d) << std::endl;
  
  return 0;
}
