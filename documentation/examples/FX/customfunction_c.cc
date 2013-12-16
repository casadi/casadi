#include "symbolic/fx/custom_function.hpp"
#include "symbolic/functor.hpp"

using namespace CasADi;

  /// Wrapper around functions
  typedef void (*CustomEvaluateCPtr)(CustomFunction &f, void* user_data);


void myEvaluate(CustomFunction &f, void* user_data) {
  DMatrix x = f.input(0);
  f.setOutput(sin(x));
}

int main(int argc, char *argv[])
{
  
  std::vector<CRSSparsity> ins;
  ins.push_back(sp_dense(1,1));
  
  std::vector<CRSSparsity> outs;
  outs.push_back(sp_dense(1,1));
  
  CustomFunction f(myEvaluate,ins,outs);
  f.init();
  
  f.setInput(2);
  f.evaluate();
  
  std::cout << "out:" << f.output() << std::endl;
  
  return 0;
}
