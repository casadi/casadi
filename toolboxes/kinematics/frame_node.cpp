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

#include "frame_node.hpp"
#include <cassert>

namespace KINEMATICS{
using namespace std;
using namespace CasADi;

FrameNode::FrameNode(const std::string& name_,const SXMatrix& q_,const SXMatrix& dq_,const SXMatrix& ddq_) : name(name_),q(q_),dq(dq_),ddq(ddq_) {
  count = 0;
}

FrameNode::FrameNode(const std::string& name_, const Frame &ref_,const SXMatrix & T) : name(name_), ref(ref_){
  R=SXMatrix(3,3);
  p=SXMatrix(3,1);
  R(0,0)=T(0,0);R(0,1)=T(0,1);R(0,2)=T(0,2);
  R(1,0)=T(1,0);R(1,1)=T(1,1);R(1,2)=T(1,2);
  R(2,0)=T(2,0);R(2,1)=T(2,1);R(2,2)=T(2,2);
  p[0]=T(0,3);
  p[1]=T(1,3);
  p[2]=T(2,3);
  q = ref.node->q;
  dq = ref.node->dq;
  ddq = ref.node->ddq;
  count = 0;
}
  
FrameNode::~FrameNode(){
  assert(count==0);
}

void FrameNode::print(std::ostream &stream){
  stream << ref;
  stream << " -> ";
  stream << name;
}

} // namespace KINEMATICS

