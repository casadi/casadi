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

%{
#include <sstream>
#include "symbolic/stl_vector_tools.hpp"
#include "symbolic/printable_object.hpp"
#include "symbolic/shared_object.hpp"
#include "symbolic/generic_type.hpp"
#include "symbolic/casadi_types.hpp"
#include "symbolic/options_functionality.hpp"
#include "symbolic/matrix/sparsity.hpp"
#include "symbolic/matrix/slice.hpp"
#include "symbolic/matrix/matrix.hpp"
#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/sx/sx_element.hpp"
#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/mx/mx.hpp"
#include "symbolic/mx/mx_tools.hpp"
#include "symbolic/fx/fx.hpp"
%}

%define %my_genericmatrix_const_typemap(Precedence,Type...) 
%typemap(in) const CasADi::GenericMatrix< Type > & (Type m) {
  if (meta< Type >::isa($input)) { // Type object get passed on as-is, and fast.
    Type* temp = static_cast< Type* >($1);
    int result = meta< Type >::get_ptr($input,temp);
    $1 = temp;
    if (!result)
      SWIG_exception_fail(SWIG_TypeError,"Type cast failed");
  } else {
    bool result=meta< Type >::as($input,m);
    if (!result)
      SWIG_exception_fail(SWIG_TypeError,meta< Type >::expected_message);
    $1 = &m;
  }
}

%typemap(typecheck,precedence=Precedence) const CasADi::GenericMatrix< Type > & { $1 = meta< Type >::isa($input) || meta< Type >::couldbe($input); }
%typemap(freearg) const CasADi::GenericMatrix< Type >  & {}

%enddef


%inline %{
template<> swig_type_info** meta< double >::name = &SWIGTYPE_p_double;
template<> swig_type_info** meta< int >::name = &SWIGTYPE_p_int;
template<> swig_type_info** meta< std::vector<double> >::name = &SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t;
template<> swig_type_info** meta< std::vector<int> >::name = &SWIGTYPE_p_std__vectorT_int_std__allocatorT_int_t_t;
template<> swig_type_info** meta< CasADi::GenericType >::name = &SWIGTYPE_p_CasADi__GenericType;
template<> swig_type_info** meta< CasADi::GenericType::Dictionary >::name = &SWIGTYPE_p_std__mapT_std__string_CasADi__GenericType_t;
template<> swig_type_info** meta< CasADi::SXElement >::name = &SWIGTYPE_p_CasADi__SXElement;
template<> swig_type_info** meta< CasADi::Matrix<CasADi::SXElement> >::name = &SWIGTYPE_p_CasADi__MatrixT_CasADi__SXElement_t;
template<> swig_type_info** meta< std::vector< CasADi::Matrix<CasADi::SXElement> > >::name = &SWIGTYPE_p_std__vectorT_CasADi__MatrixT_CasADi__SXElement_t_std__allocatorT_CasADi__MatrixT_CasADi__SXElement_t_t_t;

template<> swig_type_info** meta< std::vector< std::vector< CasADi::Matrix<CasADi::SXElement> > > >::name = &SWIGTYPE_p_std__vectorT_std__vectorT_CasADi__MatrixT_CasADi__SXElement_t_std__allocatorT_CasADi__MatrixT_CasADi__SXElement_t_t_t_std__allocatorT_std__vectorT_CasADi__MatrixT_CasADi__SXElement_t_std__allocatorT_CasADi__MatrixT_CasADi__SXElement_t_t_t_t_t;
template<> swig_type_info** meta< std::vector< CasADi::SXElement > >::name = &SWIGTYPE_p_std__vectorT_CasADi__SXElement_std__allocatorT_CasADi__SXElement_t_t;
template<> swig_type_info** meta< std::vector< std::vector< CasADi::SXElement > > >::name = &SWIGTYPE_p_std__vectorT_std__vectorT_CasADi__SXElement_std__allocatorT_CasADi__SXElement_t_t_std__allocatorT_std__vectorT_CasADi__SXElement_std__allocatorT_CasADi__SXElement_t_t_t_t;

template<> swig_type_info** meta< std::vector< CasADi::Matrix<double> > >::name = &SWIGTYPE_p_std__vectorT_CasADi__MatrixT_double_t_std__allocatorT_CasADi__MatrixT_double_t_t_t;
template<> swig_type_info** meta< std::vector< std::vector< CasADi::Matrix<double> > > >::name = &SWIGTYPE_p_std__vectorT_std__vectorT_CasADi__MatrixT_double_t_std__allocatorT_CasADi__MatrixT_double_t_t_t_std__allocatorT_std__vectorT_CasADi__MatrixT_double_t_std__allocatorT_CasADi__MatrixT_double_t_t_t_t_t;
template<> swig_type_info** meta< std::vector< CasADi::Matrix<int> > >::name = &SWIGTYPE_p_std__vectorT_CasADi__MatrixT_int_t_std__allocatorT_CasADi__MatrixT_int_t_t_t;
template<> swig_type_info** meta< std::vector< std::vector< CasADi::Matrix<int> > >  >::name = &SWIGTYPE_p_std__vectorT_std__vectorT_CasADi__MatrixT_int_t_std__allocatorT_CasADi__MatrixT_int_t_t_t_std__allocatorT_std__vectorT_CasADi__MatrixT_int_t_std__allocatorT_CasADi__MatrixT_int_t_t_t_t_t;
template<> swig_type_info** meta< CasADi::Sparsity >::name = &SWIGTYPE_p_CasADi__Sparsity;
template<> swig_type_info** meta< CasADi::Matrix<double> >::name = &SWIGTYPE_p_CasADi__MatrixT_double_t;
template<> swig_type_info** meta< CasADi::Matrix<int> >::name = &SWIGTYPE_p_CasADi__MatrixT_int_t;
template<> swig_type_info** meta< CasADi::MX >::name = &SWIGTYPE_p_CasADi__MX;
template<> swig_type_info** meta< std::vector< CasADi::MX> >::name = &SWIGTYPE_p_std__vectorT_CasADi__MX_std__allocatorT_CasADi__MX_t_t;
	template<> swig_type_info** meta< std::vector< std::vector< CasADi::MX> > >::name = &SWIGTYPE_p_std__vectorT_std__vectorT_CasADi__MX_p_std__allocatorT_CasADi__MX_p_t_t_std__allocatorT_std__vectorT_CasADi__MX_p_std__allocatorT_CasADi__MX_p_t_t_t_t;
template<> swig_type_info** meta< std::vector<std::string> >::name = &SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t;
%}



#ifdef SWIGOCTAVE
%include "meta_octave.i"
#endif

#ifdef SWIGPYTHON
%inline %{
template<> swig_type_info** meta< CasADi::DerivativeGenerator >::name = & SWIGTYPE_p_CasADi__DerivativeGenerator;
template<> swig_type_info** meta< CasADi::CustomEvaluate >::name = & SWIGTYPE_p_CasADi__CustomEvaluate;
template<> swig_type_info** meta< CasADi::Callback >::name = & SWIGTYPE_p_CasADi__Callback;
%}
%include "meta_python.i"
#endif
