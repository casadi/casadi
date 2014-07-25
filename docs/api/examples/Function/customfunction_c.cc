#include "casadi/core/function/custom_function.hpp"
#include "casadi/core/functor.hpp"

using namespace casadi;

  /// Wrapper around functions
  typedef void (*CustomEvaluateCPtr)(CustomFunction &f, void* user_data);


void myEvaluate(CustomFunction &f, void* user_data) {
  DMatrix x = f.input(0);
  f.setOutput(sin(x));
}

int main(int argc, char *argv[])
{
  
  std::vector<Sparsity> ins;
  ins.push_back(Sparsity::dense(1,1));
  
  std::vector<Sparsity> outs;
  outs.push_back(Sparsity::dense(1,1));
  
  CustomFunction f(myEvaluate,ins,outs);
  f.init();
  
  f.setInput(2);
  f.evaluate();
  
  std::cout << "out:" << f.output() << std::endl;
  
  return 0;
}
