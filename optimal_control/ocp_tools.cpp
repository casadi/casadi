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

#include "ocp_tools.hpp"
#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/sx_function.hpp"
#include "symbolic/stl_vector_tools.hpp"
#include "variable_tools.hpp"

#include <algorithm>
#include <set>
#include <queue>
#include <stack>

using namespace std;

namespace CasADi{
  
void updateDependent(SymbolicOCP& ocp){
  // Quick return if no dependent parameters
  if(ocp.pd.empty()) return;
  
  // Get the binding equations
  SX pd_def = binding(ocp.pd);
  
  // Get expressions for the variables
  SX pd = var(ocp.pd);
  
  // Sort out interdependencies
  bool reverse=false;
  substituteInPlace(pd, pd_def, reverse);

  // Create a function which evaluates the binding equations numerically
  SX ci = var(ocp.ci);
  SX pi = var(ocp.pi);
  vector<SX> f_in(2);
  f_in[0] = ci;
  f_in[1] = pi;
  SXFunction f(f_in,pd_def);
  f.init();
  
  // Evaluate the start attribute
  f.setInput(getStart(ocp.ci),0);
  f.setInput(getStart(ocp.pi),1);
  f.evaluate();
  const vector<double>& res = f.output().data();
  
  // Save to variables
  for(int k=0; k<ocp.pd.size(); ++k){
    ocp.pd[k].setStart(res[k]);
  }
}

} // namespace CasADi

