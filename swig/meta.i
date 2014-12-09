/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


%define %casadi_in_typemap_genericmatrix(xName, xType...) 
%typemap(in) const casadi::GenericMatrix< xType > & (xType m) {
  if (is_a($input, $descriptor(xType *))) {
    if (SWIG_ConvertPtr($input, (void **) &$1, $descriptor(xType *), 0) == -1) {
      SWIG_exception_fail(SWIG_TypeError,"xType cast failed ($1_type)");
    }
  } else {
    if (!to_##xName($input, &m))
      SWIG_exception_fail(SWIG_TypeError, "Input type conversion failure ($1_type)");
    $1 = &m;
  }
 }
%enddef

%define %casadi_typecheck_typemap_genericmatrix(xName, xPrec, xType...)
%typemap(typecheck,precedence=xPrec) const casadi::GenericMatrix< xType > & {
  $1 = is_a($input, $descriptor(xType *)) || to_##xName($input, 0);
 }
%enddef

%define %casadi_typemaps_genericmatrix(xName, xPrec, xType...)
%casadi_in_typemap_genericmatrix(xName, xType)
%casadi_freearg_typemap(const casadi::GenericMatrix< xType >  &)
%casadi_typecheck_typemap_genericmatrix(xName, xPrec, xType)
%enddef

#ifdef SWIGPYTHON
%include "python/meta_python.i"
#endif

#ifdef SWIGMATLAB
%include "matlab/meta_matlab.i"
#endif
