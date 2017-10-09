#include <casadi/casadi.hpp>
#include <limits>

using namespace casadi;

int main(int argc, char *argv[])
{
  DM d = DM::zeros(3,1);

  casadi_assert(!all(d).scalar(), "");
  casadi_assert(!any(d).scalar(), "");
  std::cout << d << "   all: " << all(d) << "  any: " <<  any(d) << std::endl;
  
  
  d(0) = 0;
  d(1) = 0;
  d(2) = 1;
  
  casadi_assert(!all(d).scalar(), "");
  casadi_assert(any(d).scalar(), "");
  std::cout << d << "   all: " << all(d) << "  any: " <<  any(d) << std::endl;
  
  d(0) = 1;
  d(1) = 1;
  d(2) = 1;
  
  casadi_assert(all(d).scalar(), "");
  casadi_assert(any(d).scalar(), "");
  std::cout << d << "   all: " << all(d) << "  any: " <<  any(d) << std::endl;
  
  return 0;
}
