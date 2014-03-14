#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/std_vector_tools.hpp"
#include <limits>

using namespace CasADi;

int main(int argc, char *argv[])
{
  std::cout << range(10) << std::endl;
  std::cout << range(1,11) << std::endl;
  std::cout << range(1,11,2) << std::endl;
  std::cout << range(1,11,2,5) << std::endl;
  
  return 0;
}
