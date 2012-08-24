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

#ifndef GENERIC_MATRIX_TOOLS_HPP
#define GENERIC_MATRIX_TOOLS_HPP

#include "slice.hpp"
#include "submatrix.hpp"
#include "nonzeros.hpp"
#include "crs_sparsity.hpp"
#include "../casadi_math.hpp"

namespace CasADi{


/** \brief Matlab's linspace command
*/
template<typename T>
T linspace(const GenericMatrix<T> &a, const GenericMatrix<T> &b, int nsteps);

#ifndef SWIG
template<typename T>
T linspace(const GenericMatrix<T> &a_, const GenericMatrix<T> &b_, int nsteps){
  const T& a = static_cast<const T&>(a_);
  const T& b = static_cast<const T&>(b_);
  std::vector<T> ret(nsteps);
  ret[0] = a;
  T step = (b-a)/(nsteps-1);

  for(int i=1; i<nsteps-1; ++i)
    ret[i] = ret[i-1] + step;
  
  ret[nsteps-1] = b;
  return vertcat(ret);
}
#endif // SWIG

} // namespace CasADi

#ifdef SWIG

// map the template name to the instantiated name
#define GMTT_INST(T,function_name) \
%template(function_name) CasADi::function_name< T >;

// Define template instanciations
#define GENERIC_MATRIX_TOOLS_TEMPLATES(T) \
GMTT_INST(T,linspace) \

#endif //SWIG



#endif // GENERIC_MATRIX_TOOLS_HPP

