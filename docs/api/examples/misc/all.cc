#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include <limits>

using namespace casadi;

int main(int argc, char *argv[])
{
  DMatrix d = DMatrix::zeros(3,1);

  casadi_assert(!all(d));
  casadi_assert(!any(d));
  std::cout << d << "   all: " << all(d) << "  any: " <<  any(d) << std::endl;
  
  
  d(0) = 0;
  d(1) = 0;
  d(2) = 1;
  
  casadi_assert(!all(d));
  casadi_assert(any(d));
  std::cout << d << "   all: " << all(d) << "  any: " <<  any(d) << std::endl;
  
  d(0) = 1;
  d(1) = 1;
  d(2) = 1;
  
  casadi_assert(all(d));
  casadi_assert(any(d));
  std::cout << d << "   all: " << all(d) << "  any: " <<  any(d) << std::endl;
  
  return 0;
}
