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

#include "qp_stabilizer_internal.hpp"
#include "qp_stabilizer.hpp"

using namespace std;
namespace CasADi{

QPStabilizer::QPStabilizer(){ 
}


QPStabilizer::QPStabilizer(const QPStructure & st)  {
  assignNode(new QPStabilizerInternal(st));
}

QPStabilizerInternal* QPStabilizer::operator->(){
  return (QPStabilizerInternal*)(Function::operator->());
}

const QPStabilizerInternal* QPStabilizer::operator->() const{
  return (const QPStabilizerInternal*)(Function::operator->());

}

bool QPStabilizer::checkNode() const{
  return dynamic_cast<const QPStabilizerInternal*>(get());
}

QPSolver & QPStabilizer::getSolver() {
  return (*this)->qp_solver_;
}

} // namespace CasADi
