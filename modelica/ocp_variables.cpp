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

void OCPVariables::print(ostream &stream) const{
  stream << "Time:                           " << t << endl;
  stream << "Differential states:            " << x << endl;
  stream << "Algebraic states:               " << z << endl;
  stream << "Controls:                       " << u << endl;
  stream << "Parameters:                     " << p << endl;
  stream << "Constants:                      " << c << endl;
  stream << "Dependent:                      " << d << endl;
  stream << endl;
}

  } // namespace Modelica
} // namespace CasADi

