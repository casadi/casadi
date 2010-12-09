/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "shared_object.hpp"
#include "casadi_exception.hpp"
#include <cassert>
#include <typeinfo>

using namespace std;
namespace CasADi{

SharedObject::SharedObject(){
  node = 0;
}

SharedObjectNode::SharedObjectNode(const SharedObjectNode& node){
  count = 0; // reference counter is _not_ copied
}

SharedObjectNode& SharedObjectNode::operator=(const SharedObjectNode& node){
  // do _not_ copy the reference counter
  return *this;
}

SharedObject::SharedObject(const SharedObject& ref){
  node = ref.node;
  count_up();
}

SharedObject::~SharedObject(){
  count_down();
}

void SharedObject::assignNode(SharedObjectNode* node_){
  count_down();
  node = node_;
  count_up();
}

SharedObject& SharedObject::operator=(const SharedObject& ref){
  // quick return if the old and new pointers point to the same object
  if(node == ref.node) return *this;

  // decrease the counter and delete if this was the last pointer       
  count_down();

  // save the new pointer
  node = ref.node;
  count_up();
  return *this;
}

const SharedObjectNode* SharedObject::get() const{
  return node;
}

SharedObjectNode* SharedObject::get(){
  return node;
}

bool SharedObject::isNull() const{
  return node==0;
}

void SharedObject::count_up(){
  if(node) node->count++;  
}

void SharedObject::count_down(){
  if(node && --node->count == 0){
    delete node;
    node = 0;
  }  
}

const SharedObjectNode* SharedObject::operator->() const{
  if(node==0) throw CasadiException("SharedObject::operator->(): Node is null"); 
  return node;
}

SharedObjectNode* SharedObject::operator->(){
  if(node==0) throw CasadiException("SharedObject::operator->(): Node is null"); 
  return node;
}
    
SharedObjectNode::SharedObjectNode(){
  count = 0;
}

SharedObjectNode::~SharedObjectNode(){
   assert(count==0);
}

void SharedObject::init(){
  (*this)->init();
}

void SharedObjectNode::init(){
}
    
void SharedObject::assertNode() const{
  if(isNull())
    throw CasadiException("SharedObject::assertNode(): Object is null");
}

void SharedObject::repr(std::ostream &stream) const{
  if(isNull())
    stream << 0;
  else
    (*this)->repr(stream);
}

void SharedObjectNode::repr(std::ostream &stream) const{
  // Print description by default
  print(stream);
}

void SharedObject::print(std::ostream &stream) const{
  if(isNull())    stream << "Null pointer of class \"" << typeid(this).name() << "\"";
  else           (*this)->print(stream);
}

void SharedObjectNode::print(std::ostream &stream) const{
  // Print the name of the object by default
  stream << typeid(this).name();
}

void SharedObject::makeUnique(){
  if(node && node->count>1){
    SharedObjectNode *newnode = node->clone();
    assignNode(newnode);
  }
}

SharedObjectNode* SharedObjectNode::clone() const{
  throw CasadiException(string("clone() has not defined for class") + typeid(this).name());
}
    
} // namespace CasADi
    
    