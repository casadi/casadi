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


%{
#include <sstream>
#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/printable_object.hpp"
#include "casadi/core/shared_object.hpp"
#include "casadi/core/generic_type.hpp"
#include "casadi/core/casadi_types.hpp"
#include "casadi/core/options_functionality.hpp"
#include "casadi/core/matrix/sparsity.hpp"
#include "casadi/core/matrix/slice.hpp"
#include "casadi/core/matrix/matrix.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/sx/sx_element.hpp"
#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/mx/mx.hpp"
#include "casadi/core/mx/mx_tools.hpp"
#include "casadi/core/function/function.hpp"
%}

%define %my_genericmatrix_const_typemap(Precedence,Type...) 
%typemap(in) const casadi::GenericMatrix< Type > & (Type m) {
  if (meta< Type >::isa($input)) { // Type object get passed on as-is, and fast.
    $1 = meta< Type >::get_ptr($input);
    if ($1==0) {
      SWIG_exception_fail(SWIG_TypeError,"Type cast failed");
    }
  } else {
    bool result=meta< Type >::as($input,m);
    if (!result)
      SWIG_exception_fail(SWIG_TypeError,meta< Type >::expected_message);
    $1 = &m;
  }
}

%typemap(typecheck,precedence=Precedence) const casadi::GenericMatrix< Type > & { $1 = meta< Type >::isa($input) || meta< Type >::couldbe($input); }
%typemap(freearg) const casadi::GenericMatrix< Type >  & {}

%enddef


%inline %{
template<> swig_type_info** meta< double >::name = &SWIGTYPE_p_double;
template<> swig_type_info** meta< int >::name = &SWIGTYPE_p_int;
template<> swig_type_info** meta< std::string >::name = &SWIGTYPE_p_std__string;
template<> swig_type_info** meta< std::vector<double> >::name = &SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t;
template<> swig_type_info** meta< std::vector<int> >::name = &SWIGTYPE_p_std__vectorT_int_std__allocatorT_int_t_t;
template<> swig_type_info** meta< casadi::GenericType >::name = &SWIGTYPE_p_casadi__GenericType;
template<> swig_type_info** meta< casadi::GenericType::Dictionary >::name = &SWIGTYPE_p_Dictionary;
template<> swig_type_info** meta< casadi::SXElement >::name = &SWIGTYPE_p_casadi__SXElement;
template<> swig_type_info** meta< casadi::SX >::name = &SWIGTYPE_p_casadi__MatrixT_casadi__SXElement_t;
template<> swig_type_info** meta< std::vector< casadi::SX > >::name = &SWIGTYPE_p_std__vectorT_casadi__MatrixT_casadi__SXElement_t_std__allocatorT_casadi__MatrixT_casadi__SXElement_t_t_t;

template<> swig_type_info** meta< std::vector< std::vector< casadi::SX > > >::name = &SWIGTYPE_p_std__vectorT_std__vectorT_casadi__MatrixT_casadi__SXElement_t_std__allocatorT_casadi__MatrixT_casadi__SXElement_t_t_t_std__allocatorT_std__vectorT_casadi__MatrixT_casadi__SXElement_t_std__allocatorT_casadi__MatrixT_casadi__SXElement_t_t_t_t_t;
template<> swig_type_info** meta< std::vector< casadi::SXElement > >::name = &SWIGTYPE_p_std__vectorT_casadi__SXElement_std__allocatorT_casadi__SXElement_t_t;
template<> swig_type_info** meta< std::vector< std::vector< casadi::SXElement > > >::name = &SWIGTYPE_p_std__vectorT_std__vectorT_casadi__SXElement_std__allocatorT_casadi__SXElement_t_t_std__allocatorT_std__vectorT_casadi__SXElement_std__allocatorT_casadi__SXElement_t_t_t_t;

template<> swig_type_info** meta< std::vector< casadi::Matrix<double> > >::name = &SWIGTYPE_p_std__vectorT_casadi__MatrixT_double_t_std__allocatorT_casadi__MatrixT_double_t_t_t;
template<> swig_type_info** meta< std::vector< std::vector< casadi::Matrix<double> > > >::name = &SWIGTYPE_p_std__vectorT_std__vectorT_casadi__MatrixT_double_t_std__allocatorT_casadi__MatrixT_double_t_t_t_std__allocatorT_std__vectorT_casadi__MatrixT_double_t_std__allocatorT_casadi__MatrixT_double_t_t_t_t_t;
template<> swig_type_info** meta< std::vector< casadi::Matrix<int> > >::name = &SWIGTYPE_p_std__vectorT_casadi__MatrixT_int_t_std__allocatorT_casadi__MatrixT_int_t_t_t;
template<> swig_type_info** meta< std::vector< std::vector< casadi::Matrix<int> > >  >::name = &SWIGTYPE_p_std__vectorT_std__vectorT_casadi__MatrixT_int_t_std__allocatorT_casadi__MatrixT_int_t_t_t_std__allocatorT_std__vectorT_casadi__MatrixT_int_t_std__allocatorT_casadi__MatrixT_int_t_t_t_t_t;
template<> swig_type_info** meta< casadi::Sparsity >::name = &SWIGTYPE_p_casadi__Sparsity;
template<> swig_type_info** meta< casadi::Matrix<double> >::name = &SWIGTYPE_p_casadi__MatrixT_double_t;
template<> swig_type_info** meta< casadi::Matrix<int> >::name = &SWIGTYPE_p_casadi__MatrixT_int_t;
template<> swig_type_info** meta< casadi::MX >::name = &SWIGTYPE_p_casadi__MX;
template<> swig_type_info** meta< std::vector< casadi::MX> >::name = &SWIGTYPE_p_std__vectorT_casadi__MX_std__allocatorT_casadi__MX_t_t;
	template<> swig_type_info** meta< std::vector< std::vector< casadi::MX> > >::name = &SWIGTYPE_p_std__vectorT_std__vectorT_casadi__MX_p_std__allocatorT_casadi__MX_p_t_t_std__allocatorT_std__vectorT_casadi__MX_p_std__allocatorT_casadi__MX_p_t_t_t_t;
template<> swig_type_info** meta< std::vector<std::string> >::name = &SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t;
%}

#ifdef SWIGPYTHON
%inline %{
template<> swig_type_info** meta< casadi::DerivativeGenerator >::name = & SWIGTYPE_p_casadi__DerivativeGenerator;
template<> swig_type_info** meta< casadi::CustomEvaluate >::name = & SWIGTYPE_p_casadi__CustomEvaluate;
template<> swig_type_info** meta< casadi::Callback >::name = & SWIGTYPE_p_casadi__Callback;
%}
%include "meta_python.i"
#endif

#ifdef SWIGMATLAB
%include "meta_matlab.i"
#endif
