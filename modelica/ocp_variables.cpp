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

#include "ocp_variables.hpp"
#include <algorithm>
#include <set>
#include <cassert>
#include "../casadi/casadi_exception.hpp"
#include "../casadi/stl_vector_tools.hpp"

using namespace std;
namespace CasADi{
  namespace Modelica{

OCPVariables::OCPVariables(const Variable& var){
    
  // Get all the variables
  vector<Variable> v = var;
  
  // Loop over variables
  for(vector<Variable>::iterator it=v.begin(); it!=v.end(); ++it){
    // Make sure that the variable is initialized
    switch(it->getType()){
      case TYPE_INDEPENDENT:        assert(t.isNull());     t = *it;  break;
      case TYPE_STATE:              x.push_back(*it);  break;
      case TYPE_ALGEBRAIC:          z.push_back(*it);  break;
      case TYPE_CONTROL:            u.push_back(*it);  break;
      case TYPE_PARAMETER:          p.push_back(*it);  break;
      case TYPE_CONSTANT:           c.push_back(*it);  break;
      case TYPE_DEPENDENT:          d.push_back(*it);  break;
      default: throw CasadiException("OCP::sortVariables: unknown type for " + it->getName());
    }
  }
}
    
void OCPVariables::print(ostream &stream) const{
  stream << "{" << endl;
  stream << "  t = " << t << endl;
  stream << "  x =  " << x << endl;
  stream << "  z =  " << z << endl;
  stream << "  u =  " << u << endl;
  stream << "  p =  " << p << endl;
  stream << "  c =  " << c << endl;
  stream << "  d =  " << d << endl;
  stream << "}" << endl;
}

  } // namespace Modelica
} // namespace CasADi

