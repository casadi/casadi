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

#ifdef WITH_SWIG_SPLIT
%module casadi_primitive_tools

%include "common.i"

%import "casadi_primitive.i"
#endif //WITH_SWIG_SPLIT

#ifdef SWIGPYTHON
%pythoncode%{
try:
  from numpy import pi, inf
except:
  pass

try:
  from numpy import sin, cos, tan, sqrt, log, exp, floor, ceil, fmin, fmax, sinh, cosh, tanh, arcsin, arccos, arctan, arctan2, fabs, sign, arctanh, arcsinh, arccosh, copysign
except:
  sin = lambda x: x.sin()
  cos = lambda x: x.cos()
  tan = lambda x: x.tan()
  arcsin = lambda x: x.arcsin()
  arccos = lambda x: x.arccos()
  arctan = lambda x: x.arctan()
  sqrt = lambda x: x.sqrt()
  log = lambda x: x.log()
  exp = lambda x: x.exp()
  floor = lambda x: x.floor()
  ceil = lambda x: x.ceil()
  fmin = lambda x,y: x.fmin(y)
  fmax = lambda x,y: x.fmax(y)
  sinh = lambda x: x.sinh()
  cosh = lambda x: x.cosh()
  tanh = lambda x: x.tanh()
  fabs = lambda x: x.fabs()
  sign = lambda x: x.sign()
  arctan2 = lambda x,y: x.arctan2(y)
  arctanh = lambda x: x.arctanh()
  arcsinh = lambda x: x.arcsinh()
  arccosh = lambda x: x.arccosh()
  copysign = lambda x,y: x.copysign(y)
%}
#endif // SWIGPYTHON

// Matrix tools
%include "symbolic/matrix/matrix_tools.hpp"
%include "symbolic/matrix/generic_matrix_tools.hpp"

// General tools
%include "symbolic/matrix/generic_expression_tools.hpp"

// Instantiate the functions
MATRIX_TOOLS_TEMPLATES(int)
MATRIX_TOOLS_TEMPLATES(double)
MATRIX_TOOLS_TEMPLATES(CasADi::SXElement)

GENERIC_MATRIX_TOOLS_TEMPLATES(CasADi::Matrix<int>)
GENERIC_MATRIX_TOOLS_TEMPLATES(CasADi::Matrix<double>)
GENERIC_MATRIX_TOOLS_TEMPLATES(CasADi::Matrix<CasADi::SXElement>)
GENERIC_MATRIX_TOOLS_TEMPLATES(CasADi::MX)

GENERIC_MATRIX_TOOLS_TEMPLATES_REAL_ONLY(CasADi::Matrix<double>)
GENERIC_MATRIX_TOOLS_TEMPLATES_REAL_ONLY(CasADi::Matrix<CasADi::SXElement>)
GENERIC_MATRIX_TOOLS_TEMPLATES_REAL_ONLY(CasADi::MX)


GENERIC_EXPRESSION_TOOLS_TEMPLATES(int)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(double)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(CasADi::Matrix<int>)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(CasADi::Matrix<double>)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(CasADi::Matrix<CasADi::SXElement>)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(CasADi::MX)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(CasADi::SXElement)

// Sparsity tools
%{
#include "symbolic/matrix/sparsity_tools.hpp"
%}
%include "symbolic/matrix/sparsity_tools.hpp"

// SXElement tools
%include "symbolic/sx/sx_tools.hpp"


%include "symbolic/mx/mx_tools.hpp"
