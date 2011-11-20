#include "casadi/sx/sx_tools.hpp"
#include <limits>

using namespace CasADi;

int main(int argc, char *argv[])
{

  SXMatrix x=ssym("x",4,1);
  SXMatrix s=sum_all(outer_prod((x-1),(x-1)));

  std::cout << "Default (10000)" << std::endl;
  std::cout << s << std::endl;
  
  std::cout << "Unlimited printing" << std::endl;
  SX::setMaxNumCallsInPrint(std::numeric_limits<long>::max());

  std::cout << s << std::endl;
  
  std::cout << "Limit to 10 calls" << std::endl;
  SX::setMaxNumCallsInPrint(10);
  std::cout << s << std::endl;
  
  std::cout << "Limit to 100 calls" << std::endl;
  SX::setMaxNumCallsInPrint(100);
  std::cout << s << std::endl;
  
  return 0;
}
