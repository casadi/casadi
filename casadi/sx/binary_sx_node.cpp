/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#include "binary_sx_node.hpp"
#include <cassert>

using namespace std;
namespace CasADi{

void BinarySXNode::print(ostream &stream) const{
  stringstream s0,s1;
  s0 << child[0];
  s1 << child[1];
  print_c[op](stream,s0.str(),s1.str());
}

bool BinarySXNode::isSmooth() const{
  if(op == STEP_NODE || op == FLOOR_NODE)
    return false;
  else
    return true;
}

const SX& BinarySXNode::dependent(int i) const{
  assert(i==0 || i==1);
  return child[i];
}

} // namespace CasADi
