#include <iostream>

#include "llvm/DerivedTypes.h"
#include "llvm/LLVMContext.h"
#include "llvm/Module.h"
#include "llvm/Analysis/Verifier.h"
#include "llvm/Support/IRBuilder.h"

using namespace std;
using namespace llvm;

int main(){
  LLVMContext &Context = getGlobalContext();

  
  
  cout << "hello work" << endl;

  return 0;
}

