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

#define PRECEDENCE_SX 103
#define PRECEDENCE_SXElement 102
#define PRECEDENCE_SXVector 103
#define PRECEDENCE_SXVectorVector 103
#define PRECEDENCE_SXElementVector 102
#define PRECEDENCE_SXElementVectorVector 102
#define PRECEDENCE_MX 104
#define PRECEDENCE_MXVector 105
#define PRECEDENCE_MXVectorVector 106

#define PRECEDENCE_CREATOR 150
#define PRECEDENCE_DERIVATIVEGENERATOR 21
#define PRECEDENCE_CUSTOMEVALUATE 21
#define PRECEDENCE_CALLBACK 21

#define PRECEDENCE_GENERICTYPE 22
#define PRECEDENCE_DICTIONARY 21

//#ifdef SWIG_MAIN_MODULE
%template(SXElementVector) std::vector< casadi::SXElement > ;
%template(SparsityVector) std::vector< casadi::Sparsity > ;
%template(SparsityVectorVector) std::vector< std::vector< casadi::Sparsity> > ;
%template(SXVector) std::vector<casadi::Matrix<casadi::SXElement> > ;
%template(SXVectorVector) std::vector< std::vector<casadi::Matrix<casadi::SXElement> > > ;
%template(MXVector) std::vector<casadi::MX>;
%template(MXVectorVector) std::vector< std::vector<casadi::MX> >;
%template(IMatrixVector) std::vector<casadi::Matrix<int> > ;
%template(DMatrixVector) std::vector<casadi::Matrix<double> > ;
%template(DMatrixVectorVector) std::vector< std::vector<casadi::Matrix<double> > > ;
%template(IMatrixVectorVector) std::vector< std::vector<casadi::Matrix<int> > > ;
%template(SXElementVectorVector)       std::vector<std::vector<casadi::SXElement> > ;
%template(SXElementVectorVectorVector) std::vector< std::vector<std::vector<casadi::SXElement> > > ;
//#endif //SWIG_MAIN_MODULE
/**#ifndef SWIG_MAIN_MODULE
%template() std::vector<casadi::Matrix<casadi::SXElement> > ;
%template() std::vector< std::vector<casadi::Matrix<casadi::SXElement> > > ;
%template() std::vector<casadi::MX>;
%template() std::vector<casadi::SXElement>;
%template() std::vector< std::vector<casadi::MX> >;
%template() std::vector<casadi::Matrix<int> > ;
%template() std::vector<casadi::Matrix<double> > ;
%template() std::vector< std::vector<casadi::Matrix<double> > > ;
%template() std::vector< std::vector<casadi::Matrix<int> > > ;
%template() std::vector<std::vector<casadi::SXElement> > ;
%template() std::vector< std::vector<std::vector<casadi::SXElement> > > ;
%template() std::vector< casadi::Sparsity > ;
%template() std::vector< std::vector< casadi::Sparsity> > ;
#endif //SWIG_MAIN_MODULE*/

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
%typemap(out) casadi::GenericType {
bool ret=meta<  casadi::GenericType >::toPython($1,$result);
if (!ret) {
  SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
}

%typemap(out) std::vector< casadi::GenericType > {
  PyObject* ret = PyList_New(0);
  std::vector< casadi::GenericType > & in = $1;
  for (int k=0 ; k < in.size(); ++k) {
    PyObject* rete;
    bool retv = meta< casadi::GenericType >::toPython(in[k],rete);
    if (!retv) {
      SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
    }
    PyList_Append(ret, rete);
  }
  $result = ret;
}

%typemap(out) const casadi::GenericType::Dictionary&  {
bool ret=meta<  casadi::GenericType::Dictionary >::toPython(*$1,$result);
if (!ret) {
  SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
}
#endif // SWIGPYTHON

%my_generic_const_typemap(PRECEDENCE_GENERICTYPE,casadi::GenericType)
%my_generic_const_typemap(PRECEDENCE_GENERICTYPE,std::vector< casadi::GenericType >)
#ifdef SWIGPYTHON
%my_generic_const_typemap(PRECEDENCE_DICTIONARY ,casadi::GenericType::Dictionary)
#endif

%my_creator_typemap(PRECEDENCE_CREATOR, casadi::implicitFunctionCreator);
%my_creator_typemap(PRECEDENCE_CREATOR, casadi::linearSolverCreator);

#ifdef SWIGPYTHON
%my_generic_const_typemap(PRECEDENCE_DERIVATIVEGENERATOR,casadi::DerivativeGenerator);
%my_generic_const_typemap(PRECEDENCE_CUSTOMEVALUATE,casadi::CustomEvaluate);
%my_generic_const_typemap(PRECEDENCE_CALLBACK,casadi::Callback);
#endif

#ifdef SWIGPYTHON
%my_generic_const_typemap(SWIG_TYPECHECK_DOUBLE,double);
#endif // SWIGPYTHON

%my_generic_const_typemap(PRECEDENCE_DVector,std::vector<double>);
%my_generic_const_typemap(PRECEDENCE_IVector,std::vector<int>);

%my_generic_const_typemap(PRECEDENCE_SXElement,casadi::SXElement);
%my_generic_const_typemap(PRECEDENCE_SXElementVector,std::vector< casadi::SXElement >);
%my_generic_const_typemap(PRECEDENCE_SXElementVectorVector,std::vector< std::vector< casadi::SXElement > >);

%my_generic_const_typemap(PRECEDENCE_SX,casadi::Matrix<casadi::SXElement>);
%my_genericmatrix_const_typemap(PRECEDENCE_SX,casadi::Matrix<casadi::SXElement>);

%my_generic_const_typemap(PRECEDENCE_SXVector,std::vector< casadi::Matrix<casadi::SXElement> >);
%my_generic_const_typemap(PRECEDENCE_SXVectorVector,std::vector< std::vector< casadi::Matrix<casadi::SXElement> > >);


%my_generic_const_typemap(PRECEDENCE_MX,casadi::MX);
%my_genericmatrix_const_typemap(PRECEDENCE_MX,casadi::MX);

%my_generic_const_typemap(PRECEDENCE_MXVector,std::vector< casadi::MX >);



%my_generic_const_typemap(PRECEDENCE_DMatrix,casadi::Matrix<double>);
%my_genericmatrix_const_typemap(PRECEDENCE_DMatrix,casadi::Matrix<double>);


#ifdef SWIGPYTHON
%my_generic_const_typemap(PRECEDENCE_IMatrix,casadi::Matrix<int>);
%my_genericmatrix_const_typemap(PRECEDENCE_IMatrix,casadi::Matrix<int>);
%my_generic_const_typemap(PRECEDENCE_MXVectorVector,std::vector< std::vector< casadi::MX > >);
%my_generic_const_typemap(PRECEDENCE_DMatrixVector,std::vector< casadi::Matrix<double> >);
%my_generic_const_typemap(PRECEDENCE_IMatrixVector,std::vector< casadi::Matrix<int> >);
%my_generic_const_typemap(PRECEDENCE_DMatrixVectorVector,std::vector< std::vector< casadi::Matrix<double> > >);
%my_generic_const_typemap(PRECEDENCE_IMatrixVectorVector,std::vector< std::vector< casadi::Matrix<int> > >);
#endif // SWIGPYTHON

%define %my_value_output_typemaps(Type,...)
%value_output_typemap(%arg(swig::from), %arg(SWIG_Traits_frag(Type)), %arg(Type));
%enddef

// These make OUTPUT behave like expected for non std container types
%my_value_output_typemaps(casadi::Matrix< casadi::SXElement >);
%my_value_output_typemaps(casadi::Matrix< double >);
%my_value_output_typemaps(casadi::Matrix< int >);
%my_value_output_typemaps(casadi::MX);
%my_value_output_typemaps(casadi::Sparsity);
//%my_value_output_typemaps(casadi::MXFunction);

#ifdef SWIGPYTHON
%outputRefOwn(casadi::Sparsity)

%outputRefNew(std::vector< casadi::SXElement >)
%outputRefNew(std::vector< int >)
%outputRefNew(std::vector< double >)

%outputRefOwn(casadi::Matrix< double >)
%outputRefOwn(casadi::Matrix< casadi::SXElement >)
#endif // CASADI_MODULE
#endif // SWIGPYTHON
