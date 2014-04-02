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
%module casadi_primitive

%include "common.i"

%import "casadi_main.i"
#endif //WITH_SWIG_SPLIT


#ifdef SWIGPYTHON
#ifdef WITH_SWIG_SPLIT
%pythoncode %{
import _casadi_primitive_tools as _casadi_global
%}
#endif // WITH_SWIG_SPLIT
#ifndef WITH_SWIG_SPLIT
%pythoncode %{
_casadi_global = _casadi
%}
#endif // WITH_SWIG_SPLIT
#endif // SWIGPYTHON

// Matrix typemaps class
%include "matrix.i"

// SXElement class
%include "sx.i"

// MX class
%include "mx.i"
