/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

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

  // Make the module, which holds all the code.
  llvm::Module *TheModule = new llvm::Module("my cool jit", llvm::getGlobalContext());

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
  
  // Two arguments in and two references
  std::vector<const llvm::Type*> genArg(4);
  genArg[0] = genArg[1] = llvm::Type::getDoubleTy(llvm::getGlobalContext());
  genArg[2] = genArg[3] = llvm::Type::getDoublePtrTy(llvm::getGlobalContext());
  
  // Unary operation
  llvm::FunctionType *unaryFun = llvm::FunctionType::get(llvm::Type::getDoubleTy(llvm::getGlobalContext()),unaryArg, false);

  // Binary operation
  llvm::FunctionType *binaryFun = llvm::FunctionType::get(llvm::Type::getDoubleTy(llvm::getGlobalContext()),binaryArg, false);

  // More generic operation, return by reference
  llvm::FunctionType *genFun = llvm::FunctionType::get(llvm::Type::getVoidTy(llvm::getGlobalContext()),genArg, false);

  // Declare sin
  llvm::Function *sin_ = llvm::Function::Create(unaryFun, llvm::Function::ExternalLinkage, "sin", TheModule);
  
  // Declare my function
  llvm::Function *myfun = llvm::Function::Create(genFun, llvm::Function::ExternalLinkage, "myfcn", TheModule);

  // Create a new basic block to start insertion into.
  llvm::BasicBlock *BB = llvm::BasicBlock::Create(llvm::getGlobalContext(), "entry", myfun);
  Builder.SetInsertPoint(BB);

  // Set names for all arguments.
  llvm::Function::arg_iterator AI = myfun->arg_begin();
  AI->setName("x1");
  llvm::Value *x1 = AI;
  AI++;
  AI->setName("x2");
  llvm::Value *x2 = AI;
  AI++;
  AI->setName("r1");
  llvm::Value *r1 = AI;
  AI++;
  AI->setName("r2");
  llvm::Value *r2 = AI;
  
  llvm::Value *five = llvm::ConstantFP::get(llvm::getGlobalContext(), llvm::APFloat(5.0));
  llvm::Value *x1_plus_5 = Builder.CreateFAdd(x1, five, "x1_plus_5");
  
  // Call the sine function
  std::vector<llvm::Value*> sinarg(1,x2);
  llvm::Value* sin_x2 = Builder.CreateCall(sin_, sinarg.begin(), sinarg.end(), "callsin");
  
  // Set values
  llvm::StoreInst *what_is_this1 = Builder.CreateStore(sin_x2,r1);
  llvm::StoreInst *what_is_this2 = Builder.CreateStore(x1_plus_5,r2);

  // Finish off the function.
  Builder.CreateRetVoid();

  // Validate the generated code, checking for consistency.
  verifyFunction(*myfun);

  // Optimize the function.
  TheFPM->run(*myfun);

  // Print out all of the generated code.
  TheModule->dump();

  // JIT the function
  double x1_val = 10;
  double x2_val = 20;
  double r1_val = -1;
  double r2_val = -1;
  typedef void (*GenType)(double,double,double*,double*);
  GenType FP = GenType(intptr_t(TheExecutionEngine->getPointerToFunction(myfun)));

  FP(x1_val,x2_val,&r1_val,&r2_val);

  printf("r1 = %g\n", r1_val);
  printf("r2 = %g\n", r2_val);
  

  return 0;
}
