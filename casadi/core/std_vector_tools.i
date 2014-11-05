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
#ifndef CASADI_STD_VECTOR_TOOLS_I
#define CASADI_STD_VECTOR_TOOLS_I

%include <casadi/core/shared_object.i>

%{
#include <casadi/core/std_vector_tools.hpp>
%}
%include <casadi/core/std_vector_tools.hpp>

VECTOR_TOOLS_TEMPLATES(int)
VECTOR_TOOLS_TEMPLATES(double)
%define VECTOR_REPR(type)
%extend std::vector< type >{
  std::string SWIG_REPR(){ return casadi::getRepresentation(*$self); }
  std::string SWIG_STR(){ return casadi::getDescription(*$self); }
};
%enddef

#ifdef SWIGPYTHON
%rename(SWIG_STR) getDescription;
#endif // SWIGPYTHON

#endif // CASADI_STD_VECTOR_TOOLS_I
