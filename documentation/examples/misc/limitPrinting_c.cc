#include "casadi/sx/sx_tools.hpp"

using namespace CasADi;

main(int argc, char *argv[])
{

  SXMatrix x=symbolic("x",4,1);
  SXMatrix s=sum_all(outer_prod((x-1),(x-1)));

  std::cout << "Unlimited printing" << std::endl;
  std::cout << s << std::endl;
  
  std::cout << "Limit to 10 characters" << std::endl;
  limitedTo(std::cout,10) << s << std::endl;
  
  std::cout << "Limit to 100 characters" << std::endl;
  limitedTo(std::cout,100) << s << std::endl;

}
