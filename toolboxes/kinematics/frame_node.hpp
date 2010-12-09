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

#ifndef FRAME_NODE_HPP
#define FRAME_NODE_HPP

#include "frame.hpp"

namespace KINEMATICS{
using namespace CasADi;

/** \brief Internal class to make Frame trees reference-safe

*/
class FrameNode{
  public:
    friend class Frame;
      
/** \brief  Constructor */
    explicit FrameNode(const std::string& name, const SXMatrix & q,const SXMatrix & dq,const SXMatrix & ddq);
    FrameNode(const std::string& name, const Frame &ref,const SXMatrix & T);
    ~FrameNode();

/** \brief  Print */
    void print(std::ostream &stream);

    
  protected:
  std::string name;
  int count;
  Frame ref; 
  SXMatrix q;
  SXMatrix dq;
  SXMatrix ddq;
  SXMatrix R; // the 3x3 rotation matrix
  SXMatrix p; // 3 vectors
};

} // namespace KINEMATICS

#endif //FRAME_NODE_HPP

