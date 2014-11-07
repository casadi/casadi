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
#ifndef CASADI_IO_SCHEME_I
#define CASADI_IO_SCHEME_I

%include <casadi/core/shared_object.i>

%rename(__call__original__) casadi::IOScheme::operator();

%include <casadi/core/function/io_scheme.hpp>

#ifdef SWIGPYTHON
%extend casadi::IOScheme {
%template(__call__original__) operator()< casadi::Sparsity >;
%template(__call__original__) operator()< casadi::MX> ;
%template(__call__original__) operator()< casadi::Matrix<casadi::SXElement> >;
%template(__call__original__) operator()< casadi::Matrix<double> >;

%pythoncode %{
  def __call__(self,*dummy,**kwargs):
    if len(dummy)>1: return self.__call__original__(dummy[0],dummy[1:])
    return self.__call__original__(kwargs.keys(),kwargs.values())
%}
}
#endif

#endif // CASADI_IO_SCHEME_I
