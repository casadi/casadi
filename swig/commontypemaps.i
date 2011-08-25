// Lower value means wil be checked first
#define PRECEDENCE_DVector 98
#define PRECEDENCE_IVector 99

#define PRECEDENCE_DMatrix 100
#define PRECEDENCE_DMatrixVector 101
#define PRECEDENCE_SXMatrix 102
#define PRECEDENCE_SX 103
#define PRECEDENCE_SXMatrixVector 103
#define PRECEDENCE_MX 104
#define PRECEDENCE_MXVector 105
#define PRECEDENCE_PAIR_SLICE_SLICE 204
#define PRECEDENCE_SLICE 205
#define PRECEDENCE_IndexVector 210
#define PRECEDENCE_PAIR_IVector_IVector 206
#define PRECEDENCE_GENERICTYPE 22
#define PRECEDENCE_DICTIONARY 21


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

%my_generic_const_typemap(PRECEDENCE_DVector,std::vector<double>);
%my_generic_const_typemap(PRECEDENCE_IVector,std::vector<int>);

%my_generic_const_typemap(PRECEDENCE_SX,CasADi::SX);

%my_generic_const_typemap(PRECEDENCE_SXMatrix,CasADi::Matrix<CasADi::SX>);

%template(SXMatrixVector)       std::vector<CasADi::Matrix<CasADi::SX> > ;
%template(SXMatrixVectorVector) std::vector< std::vector<CasADi::Matrix<CasADi::SX> > > ;

%my_generic_const_typemap(PRECEDENCE_SXMatrixVector,std::vector< CasADi::Matrix<CasADi::SX> >);



%my_generic_const_typemap(PRECEDENCE_MX,CasADi::MX);

%template(MXVector) std::vector<CasADi::MX>;

%my_generic_const_typemap(PRECEDENCE_MXVector,std::vector< CasADi::MX >);



%my_generic_const_typemap(PRECEDENCE_DMatrix,CasADi::Matrix<double>);

%template(DMatrixVector)       std::vector<CasADi::Matrix<double> > ;
%template(DMatrixVectorVector) std::vector< std::vector<CasADi::Matrix<double> > > ;

#ifdef SWIGPYTHON
%my_generic_const_typemap(PRECEDENCE_DMatrixVector,std::vector< CasADi::Matrix<double> >);
#endif // SWIGPYTHON


#ifdef SWIGPYTHON
%outputRefOwn(CasADi::CRSSparsity)
%outputRefOwn(std::vector< CasADi::SX >)
%outputRefOwn(std::vector< int >)
%outputRefOwn(std::vector< double >)
%outputRefOwn(CasADi::Matrix< double >)
%outputRefOwn(CasADi::Matrix< CasADi::SX >)
#endif // SWIGPYTHON
