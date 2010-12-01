#include "c_function.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <dlfcn.h>
#include <cassert>

namespace CasADi{

using namespace std;

CFunction::CFunction(){
}

CFunction::CFunction(CFunctionWrapper c_fcn){
  assignNode(new CFunctionNode(c_fcn));
}

CFunctionNode* CFunction::operator->(){
  return (CFunctionNode*)(FX::operator->());
}

const CFunctionNode* CFunction::operator->() const{
   return (const CFunctionNode*)(FX::operator->()); 
}
  
void CFunction::assertNode() const{
  if(!dynamic_cast<const CFunctionNode*>(get()))
    throw CasadiException("CFunction::assertNode");
}

CFunctionNode::CFunctionNode(CFunctionWrapper c_fcn) : evaluate_(c_fcn){
  user_data_ = 0;
}

CFunctionNode::~CFunctionNode(){
  
}

void CFunctionNode::setUserData(void* user_data){
  user_data_ = user_data;
}

void CFunctionNode::evaluate(int fsens_order, int asens_order){
  if(evaluate_==0) throw CasadiException("CFunctionNode::evaluate: pointer is null");
  ref_.assignNode(this); // make the reference point to this object
  evaluate_(ref_,fsens_order, asens_order,user_data_);  
  ref_.assignNode(0); // make sure to remove the reference, otherwise the object will never be deleted
}

void CFunctionNode::init(){
  FXNode::init();
}


} // namespace CasADi

