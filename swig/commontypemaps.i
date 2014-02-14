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

// Lower value means wil be checked first
#define PRECEDENCE_IVector 92
#define PRECEDENCE_PAIR_SLICE_SLICE 93
// Why are SLICE and IndexVector the same precedence?
// To circumvent an issue with typemap precedence.
// Originally, we had slice precedence < IndexList precedence, but this caused the following order:
//    indexed(Slice,Slice)
//    indexed(Slice,Martix<int>)
//    indexed(IndexList,IndexList)
//    indexed(IndexList,Martix<int>)
// While we intend it to be:
//    indexed(Slice,Slice)
//    indexed(IndexList,IndexList)
//    indexed(Slice,Martix<int>)
//    indexed(IndexList,Martix<int>)
#define PRECEDENCE_SLICE 94
#define PRECEDENCE_IndexVector 94
#define PRECEDENCE_PAIR_IVector_IVector 96

#define PRECEDENCE_IMatrix 97
#define PRECEDENCE_IMatrixVector 98
#define PRECEDENCE_IMatrixVectorVector 98

#define PRECEDENCE_DVector 99

#define PRECEDENCE_DMatrix 100
#define PRECEDENCE_DMatrixVector 101
#define PRECEDENCE_DMatrixVectorVector 101

#define PRECEDENCE_SXMatrix 103
#define PRECEDENCE_SX 102
#define PRECEDENCE_SXMatrixVector 103
#define PRECEDENCE_SXMatrixVectorVector 103
#define PRECEDENCE_SXVector 102
#define PRECEDENCE_SXVectorVector 102
#define PRECEDENCE_MX 104
#define PRECEDENCE_MXVector 105
#define PRECEDENCE_MXVectorVector 106

#define PRECEDENCE_CREATOR 150
#define PRECEDENCE_DERIVATIVEGENERATOR 21
#define PRECEDENCE_CUSTOMEVALUATE 21
#define PRECEDENCE_CALLBACK 21

#define PRECEDENCE_GENERICTYPE 22
#define PRECEDENCE_DICTIONARY 21

#ifdef SWIG_MAIN_MODULE
%template(SXVector) std::vector< CasADi::SX > ;
%template(SXMatrixVector) std::vector<CasADi::Matrix<CasADi::SX> > ;
%template(SXMatrixVectorVector) std::vector< std::vector<CasADi::Matrix<CasADi::SX> > > ;
%template(MXVector) std::vector<CasADi::MX>;
%template(MXVectorVector) std::vector< std::vector<CasADi::MX> >;
%template(IMatrixVector) std::vector<CasADi::Matrix<int> > ;
%template(DMatrixVector) std::vector<CasADi::Matrix<double> > ;
%template(DMatrixVectorVector) std::vector< std::vector<CasADi::Matrix<double> > > ;
%template(IMatrixVectorVector) std::vector< std::vector<CasADi::Matrix<int> > > ;
%template(SXVectorVector)       std::vector<std::vector<CasADi::SX> > ;
%template(SXVectorVectorVector) std::vector< std::vector<std::vector<CasADi::SX> > > ;
#endif //SWIG_MAIN_MODULE
#ifndef SWIG_MAIN_MODULE
%template() std::vector<CasADi::Matrix<CasADi::SX> > ;
%template() std::vector< std::vector<CasADi::Matrix<CasADi::SX> > > ;
%template() std::vector<CasADi::MX>;
%template() std::vector<CasADi::SX>;
%template() std::vector< std::vector<CasADi::MX> >;
%template() std::vector<CasADi::Matrix<int> > ;
%template() std::vector<CasADi::Matrix<double> > ;
%template() std::vector< std::vector<CasADi::Matrix<double> > > ;
%template() std::vector< std::vector<CasADi::Matrix<int> > > ;
%template() std::vector<std::vector<CasADi::SX> > ;
%template() std::vector< std::vector<std::vector<CasADi::SX> > > ;
#endif //SWIG_MAIN_MODULE

#ifdef CASADI_MODULE
#ifdef SWIGPYTHON
%typemap(in) int (int m) {
  bool result=meta< int >::as($input,m);
  if (!result)
    SWIG_exception_fail(SWIG_TypeError,meta< int >::expected_message);
  $1 = m;
}

%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) int { $1 = meta< int >::isa($input) || meta< int >::couldbe($input); }
%typemap(freearg) int {}

#endif //SWIGPYTHON

#ifdef SWIGPYTHON
%typemap(in) double (double m) {
  bool result=meta< double >::as($input,m);
  if (!result)
    SWIG_exception_fail(SWIG_TypeError,meta< double >::expected_message);
  $1 = m;
}

%typemap(typecheck,precedence=SWIG_TYPECHECK_DOUBLE) double { $1 = meta< double >::isa($input) || meta< double >::couldbe($input); }
%typemap(freearg) double {}

#endif //SWIGPYTHON

#ifdef SWIGPYTHON
%typemap(out) CasADi::GenericType {
bool ret=meta<  CasADi::GenericType >::toPython($1,$result);
if (!ret) {
  SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
}

%typemap(out) const CasADi::GenericType::Dictionary&  {
bool ret=meta<  CasADi::GenericType::Dictionary >::toPython(*$1,$result);
if (!ret) {
  SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
}
#endif // SWIGPYTHON

%my_generic_const_typemap(PRECEDENCE_GENERICTYPE,CasADi::GenericType)
#ifdef SWIGPYTHON
%my_generic_const_typemap(PRECEDENCE_DICTIONARY ,CasADi::GenericType::Dictionary)
#endif

%my_creator_typemap(PRECEDENCE_CREATOR, CasADi::implicitFunctionCreator);
%my_creator_typemap(PRECEDENCE_CREATOR, CasADi::linearSolverCreator);

#ifdef SWIGPYTHON
%my_generic_const_typemap(PRECEDENCE_DERIVATIVEGENERATOR,CasADi::DerivativeGenerator);
%my_generic_const_typemap(PRECEDENCE_CUSTOMEVALUATE,CasADi::CustomEvaluate);
%my_generic_const_typemap(PRECEDENCE_CALLBACK,CasADi::Callback);
#endif

#ifdef SWIGPYTHON
%my_generic_const_typemap(SWIG_TYPECHECK_DOUBLE,double);
#endif // SWIGPYTHON

%my_generic_const_typemap(PRECEDENCE_DVector,std::vector<double>);
%my_generic_const_typemap(PRECEDENCE_IVector,std::vector<int>);

%my_generic_const_typemap(PRECEDENCE_SX,CasADi::SX);
%my_generic_const_typemap(PRECEDENCE_SXVector,std::vector< CasADi::SX >);
%my_generic_const_typemap(PRECEDENCE_SXVectorVector,std::vector< std::vector< CasADi::SX > >);

%my_generic_const_typemap(PRECEDENCE_SXMatrix,CasADi::Matrix<CasADi::SX>);
%my_genericmatrix_const_typemap(PRECEDENCE_SXMatrix,CasADi::Matrix<CasADi::SX>);

%my_generic_const_typemap(PRECEDENCE_SXMatrixVector,std::vector< CasADi::Matrix<CasADi::SX> >);
%my_generic_const_typemap(PRECEDENCE_SXMatrixVectorVector,std::vector< std::vector< CasADi::Matrix<CasADi::SX> > >);


%my_generic_const_typemap(PRECEDENCE_MX,CasADi::MX);
%my_genericmatrix_const_typemap(PRECEDENCE_MX,CasADi::MX);

%my_generic_const_typemap(PRECEDENCE_MXVector,std::vector< CasADi::MX >);



%my_generic_const_typemap(PRECEDENCE_DMatrix,CasADi::Matrix<double>);
%my_genericmatrix_const_typemap(PRECEDENCE_DMatrix,CasADi::Matrix<double>);


#ifdef SWIGPYTHON
%my_generic_const_typemap(PRECEDENCE_IMatrix,CasADi::Matrix<int>);
%my_genericmatrix_const_typemap(PRECEDENCE_IMatrix,CasADi::Matrix<int>);
%my_generic_const_typemap(PRECEDENCE_MXVectorVector,std::vector< std::vector< CasADi::MX > >);
%my_generic_const_typemap(PRECEDENCE_DMatrixVector,std::vector< CasADi::Matrix<double> >);
%my_generic_const_typemap(PRECEDENCE_IMatrixVector,std::vector< CasADi::Matrix<int> >);
%my_generic_const_typemap(PRECEDENCE_DMatrixVectorVector,std::vector< std::vector< CasADi::Matrix<double> > >);
%my_generic_const_typemap(PRECEDENCE_IMatrixVectorVector,std::vector< std::vector< CasADi::Matrix<int> > >);
#endif // SWIGPYTHON

%define %my_value_output_typemaps(Type,...)
%value_output_typemap(%arg(swig::from), %arg(SWIG_Traits_frag(Type)), %arg(Type));
%enddef

// These make OUTPUT behave like expected for non std container types
%my_value_output_typemaps(CasADi::Matrix< CasADi::SX >);
%my_value_output_typemaps(CasADi::Matrix< double >);
%my_value_output_typemaps(CasADi::Matrix< int >);
%my_value_output_typemaps(CasADi::MX);
%my_value_output_typemaps(CasADi::CRSSparsity);
//%my_value_output_typemaps(CasADi::MXFunction);

#ifdef SWIGPYTHON
%outputRefOwn(CasADi::CRSSparsity)
%outputRefOwn(std::vector< CasADi::SX >)

%outputRefOwn(std::vector< int >)
%outputRefOwn(std::vector< double >)

%outputRefOwn(CasADi::Matrix< double >)
%outputRefOwn(CasADi::Matrix< CasADi::SX >)
#endif // CASADI_MODULE
#endif // SWIGPYTHON
