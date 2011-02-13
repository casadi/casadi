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

#include "multiple_output.hpp"
#include "../fx/fx_internal.hpp"
#include "../stl_vector_tools.hpp"

using namespace std;

namespace CasADi{

MultipleOutput::MultipleOutput(){
}

MultipleOutput::~MultipleOutput(){
}

OutputNode::OutputNode(const MX& parent) : parent_(parent){
  // Get the node of the parent
  MultipleOutput* p = dynamic_cast<MultipleOutput*>(parent_.get());
  casadi_assert(p!=0);
  
  // Save a pointer to the object in the parent class
  p->children_.insert(this);
}

OutputNode::~OutputNode(){
  // Get the node of the parent
  MultipleOutput* p = static_cast<MultipleOutput*>(parent_.get());
  
  // Remove the parent from the parent
  p->children_.erase(this);
}
  

} // namespace CasADi
