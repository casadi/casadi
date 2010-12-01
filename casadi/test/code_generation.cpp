#include "../fx/external_function.hpp"
#include "../fx/sx_function.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

using namespace std;
using namespace CasADi;

int main(){

  SX in("out");
  SX out("out");
  SXFunction fcn(in,out);

  // Generate c code for the function
  fcn->generateCode("test.c");

  // Compile the source into a dll
  int flag = system("gcc -fpic -shared test.c -o test.so");

  // delete the c source
  flag = system("rm test.c");

  // Create an external function
  ExternalFunction cg("test.so");

  // Call the external functions
  cg->evaluate();

//  cg.evaluateFwd();

  // delete the binary
  flag = system("rm test.so");


  return 0;
}
