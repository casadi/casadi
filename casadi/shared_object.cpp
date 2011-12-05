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

#include <typeinfo>
#include <cassert>

using namespace std;
namespace CasADi{

SharedObject::SharedObject(){
  node = 0;
}

SharedObjectNode::SharedObjectNode(const SharedObjectNode& node){
  is_init_ = node.is_init_;
  count = 0; // reference counter is _not_ copied
}

SharedObjectNode& SharedObjectNode::operator=(const SharedObjectNode& node){
  is_init_ = node.is_init_;
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

void SharedObject::assignNodeNoCount(SharedObjectNode* node_){
  node = node_;
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
  casadi_assert(!isNull());
  return node;
}

SharedObjectNode* SharedObject::operator->(){
  casadi_assert(!isNull());
  return node;
}
    
SharedObjectNode::SharedObjectNode(){
  is_init_ = false;
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
    
bool SharedObject::checkNode() const{
  return dynamic_cast<const SharedObjectNode*>(get());
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

void SharedObject::makeUnique(bool clone_members){
  std::map<SharedObjectNode*,SharedObject> already_copied;
  makeUnique(already_copied,clone_members);
}

void SharedObject::makeUnique(std::map<SharedObjectNode*,SharedObject>& already_copied, bool clone_members){
  if(node && node->count>1){
    // First find out if the expression has already been copied
    std::map<SharedObjectNode*,SharedObject>::iterator it = already_copied.find(node);
    
    if(it==already_copied.end()){
      // If the expression has not yet been copied
      SharedObjectNode *newnode = node->clone();
      
      // Copy the data members
      if(clone_members) newnode->deepCopyMembers(already_copied);
      
      // Initialize object if parent was initialized
      if(isInit() && !newnode->is_init_){
        newnode->init();
      }
      
      // Assign cloned node to object
      assignNode(newnode);
    } else {
      // Use an existing copy
      assignNode(it->second.get());
    }
  }
}

SharedObject SharedObject::clone() const{
  SharedObject ret;
  if(!isNull()){
    ret.assignNode((*this)->clone());
  }
  return ret;
}

void SharedObject::swap(SharedObject& other){
  SharedObject temp = *this;
  *this = other;
  other = temp;
}

int SharedObject::getCount() const{
  return (*this)->getCount();
}

int SharedObjectNode::getCount() const{
  return count;
}

void SharedObjectNode::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
}

bool SharedObject::isInit() const{
  return (*this)->isInit();
}

void SharedObject::assertInit() const{
  (*this)->assertInit();
}

bool SharedObjectNode::isInit() const{
  return is_init_;
}

void SharedObjectNode::assertInit() const{
  casadi_assert_message(isInit(),"You must first initialize a Shared Object before you can use it." << std::endl <<  "Use something like f.init()");
}


} // namespace CasADi
    
    
