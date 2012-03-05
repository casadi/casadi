#include "llvm/DerivedTypes.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/ExecutionEngine/JIT.h"
#include "llvm/LLVMContext.h"
#include "llvm/Module.h"
#include "llvm/PassManager.h"
#include "llvm/Analysis/Verifier.h"
#include "llvm/Target/TargetData.h"
#include "llvm/Target/TargetSelect.h"
#include "llvm/Transforms/Scalar.h"
#include "llvm/Support/IRBuilder.h"

#include <cstdio>
#include <string>
#include <map>
#include <vector>

static llvm::IRBuilder<> Builder(llvm::getGlobalContext());

int main() {
  llvm::InitializeNativeTarget();
  llvm::LLVMContext &Context = llvm::getGlobalContext();

  // Make the module, which holds all the code.
  llvm::Module *TheModule = new llvm::Module("my cool jit", Context);

  // Create the JIT.  This takes ownership of the module.
  std::string ErrStr;
  llvm::ExecutionEngine *TheExecutionEngine = llvm::EngineBuilder(TheModule).setErrorStr(&ErrStr).create();
  if (!TheExecutionEngine) {
    fprintf(stderr, "Could not create ExecutionEngine: %s\n", ErrStr.c_str());
    exit(1);
  }

  llvm::FunctionPassManager OurFPM(TheModule);

  // Set up the optimizer pipeline.  Start with registering info about how the
  // target lays out data structures.
  OurFPM.add(new llvm::TargetData(*TheExecutionEngine->getTargetData()));
  // Do simple "peephole" optimizations and bit-twiddling optzns.
  OurFPM.add(llvm::createInstructionCombiningPass());
  // Reassociate expressions.
  OurFPM.add(llvm::createReassociatePass());
  // Eliminate Common SubExpressions.
  OurFPM.add(llvm::createGVNPass());
  // Simplify the control flow graph (deleting unreachable blocks, etc).
  OurFPM.add(llvm::createCFGSimplificationPass());

  OurFPM.doInitialization();

  // Set the global so the code gen can use this.
  llvm::FunctionPassManager *TheFPM = &OurFPM;

  // Single argument
  std::vector<const llvm::Type*> unaryArg(1,llvm::Type::getDoubleTy(llvm::getGlobalContext()));

  // Two arguments
  std::vector<const llvm::Type*> binaryArg(2,llvm::Type::getDoubleTy(llvm::getGlobalContext()));
  
  // Unary operation
  llvm::FunctionType *unaryFun = llvm::FunctionType::get(llvm::Type::getDoubleTy(llvm::getGlobalContext()),unaryArg, false);

  // Binary operation
  llvm::FunctionType *binaryFun = llvm::FunctionType::get(llvm::Type::getDoubleTy(llvm::getGlobalContext()),binaryArg, false);

  // Declare sin
  llvm::Function *sin_ = llvm::Function::Create(unaryFun, llvm::Function::ExternalLinkage, "sin", TheModule);
  
  // Declare my function
  llvm::Function *myfun = llvm::Function::Create(binaryFun, llvm::Function::ExternalLinkage, "myfcn", TheModule);

  // Create a new basic block to start insertion into.
  llvm::BasicBlock *BB = llvm::BasicBlock::Create(llvm::getGlobalContext(), "entry", myfun);
  Builder.SetInsertPoint(BB);

  // Set names for all arguments.
  llvm::Function::arg_iterator AI = myfun->arg_begin();
  AI->setName("xx");
  llvm::Value *x1 = AI;
  AI++;
  AI->setName("yy");
  llvm::Value *x2 = AI;
  
  llvm::Value *five = llvm::ConstantFP::get(llvm::getGlobalContext(), llvm::APFloat(5.0));
  llvm::Value *x1_plus_5 = Builder.CreateAdd(x1, five, "x1_plus_5");
  
  // Call the sine function
  std::vector<llvm::Value*> sinarg(1,x2);
  llvm::Value* sin_x2 = Builder.CreateCall(sin_, sinarg.begin(), sinarg.end(), "callsin");
  
  // Finish off the function.
  Builder.CreateRet(x1_plus_5);

  // Validate the generated code, checking for consistency.
  verifyFunction(*myfun);

  // Optimize the function.
  TheFPM->run(*myfun);

  // JIT the function, returning a function pointer.
  void *FPtr = TheExecutionEngine->getPointerToFunction(myfun);
      
  // Cast it to the right type
  double (*FP)(double,double) = (double (*)(double,double))(intptr_t)FPtr;
  fprintf(stderr, "Evaluated to %f\n", FP(10.0,20.0));
  
  // Print out all of the generated code.
  TheModule->dump();

  return 0;
}
