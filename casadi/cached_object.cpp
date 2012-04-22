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

#include "cached_object.hpp"

using namespace std;

namespace CasADi{
  CachedObject::CachedObject(){
  }

  CachedObject::~CachedObject(){
  }
  
  const CachedObjectNode* CachedObject::operator->() const{
    return static_cast<const CachedObjectNode*>(SharedObject::operator->());
  }

  CachedObjectNode* CachedObject::operator->(){
    return static_cast<CachedObjectNode*>(SharedObject::operator->());
  }
  
  bool CachedObject::checkNode() const{
    return dynamic_cast<const CachedObjectNode*>(get())!=0;
  }  
  
  CachedObjectNode::CachedObjectNode(){
  }
  
  CachedObjectNode::~CachedObjectNode(){
    // Owner from cache references
    for(vector<WeakRef*>::iterator it=weak_.begin(); it!=weak_.end(); ++it){
      if(*it){
	(*it)->clear();
      }
    }
  }
  
  int CachedObjectNode::regRef(WeakRef* ref){
    // Check if there are unused locations
    if(!unused_.empty()){
      // Use an existing cache location
      int weak_loc = unused_.top();
      unused_.pop();
      
      // Save to cache
      weak_[weak_loc] = ref;
      return weak_loc;
    } else {
      // Create a new cache location at the end
      int weak_loc = weak_.size();
      weak_.push_back(ref);
      return weak_loc;
    }
  }
  
  void CachedObjectNode::unregRef(int weak_loc){
    // Remove from cache
    weak_[weak_loc] = 0;
    
    // Save to list of unused locations
    unused_.push(weak_loc);
  }

  WeakRef::WeakRef() : owner_(0){
  }
  
  WeakRef::~WeakRef(){
    if(owner_){
      owner_->unregRef(weak_loc_);
    }
  }

  WeakRef::WeakRef(const WeakRef& obj) : owner_(obj.owner_){
    // Register the reference
    if(owner_){
      weak_loc_ = owner_->regRef(this);
    }
  }
  
  WeakRef& WeakRef::operator=(const WeakRef& obj){
    // if same owner, do nothing
    if(obj.owner_ != owner_){
      // Unregister the existing owner, if any
      if(owner_){
	owner_->unregRef(weak_loc_);
      }
      
      // Give a new owner and register
      owner_ = obj.owner_;
      if(owner_){
	weak_loc_ = owner_->regRef(this);
      }
    }
    return *this;
  }
  
  WeakRef::WeakRef(CachedObject& obj){
    owner_ = static_cast<CachedObjectNode*>(obj.get());
    if(owner_){
      weak_loc_ = owner_->regRef(this);
    }
  }

  WeakRef& WeakRef::operator=(CachedObject& obj){
    // if same owner, do nothing
    if(obj.get() != owner_){
      // Unregister the existing owner, if any
      if(owner_){
	owner_->unregRef(weak_loc_);
      }
      
      // Give a new owner and register
      owner_ = static_cast<CachedObjectNode*>(obj.get());
      if(owner_){
	weak_loc_ = owner_->regRef(this);
      }
    }
    return *this;
  }
  
  void WeakRef::clear(){
    owner_ = 0;
  }
  
  void WeakRef::repr(std::ostream &stream) const{
    stream << "WeakRef(";
    if(owner_==0){
      stream << "NULL";
    } else {
      owner_->repr(stream);
    }
    stream << ")";
  }

  void WeakRef::print(std::ostream &stream) const{
    if(owner_==0){
      stream << "NULL";
    } else {
      owner_->print(stream);
    }
  }
  
} // namespace CasADi

