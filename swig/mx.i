%{
#include "casadi/mx/mx.hpp"
#include "casadi/mx/mx_tools.hpp"
%}

%my_generic_const_typemap(PRECEDENCE_MX,CasADi::MX);

%include "casadi/mx/mx.hpp"

%inline %{
template<> swig_type_info** meta< CasADi::MX >::name = &SWIGTYPE_p_CasADi__MX;
template<> swig_type_info** meta< std::vector< CasADi::MX> >::name = &SWIGTYPE_p_std__vectorT_CasADi__MX_std__allocatorT_CasADi__MX_t_t;

template<> swig_type_info** meta< std::pair< CasADi::MX, std::vector< CasADi::MX> > >::name = &SWIGTYPE_p_std__pairT_CasADi__MX_std__vectorT_CasADi__MX_std__allocatorT_CasADi__MX_t_t_t;
%}

%inline %{
template<> char meta< CasADi::MX >::expected_message[] = "Expecting (MX, numberarray)";
%}

/// CasADi::MX
#ifdef SWIGPYTHON
%inline %{

template <>
bool meta< CasADi::MX >::couldbe(PyObject * p) {
  return (meta< CasADi::MX >::isa(p) || meta< CasADi::Matrix<double> >::couldbe(p) );
}

template <>
int meta< CasADi::MX >::as(PyObject * p,CasADi::MX &m) {
  NATIVERETURN(CasADi::MX,m)
  if(meta< CasADi::Matrix<double> >::couldbe(p)) {
    CasADi::DMatrix mt;
    bool result=meta< CasADi::Matrix<double> >::as(p,mt);
    if (!result)
      return false;
    m = CasADi::MX(mt);
    return true;
  }
  return false;
}
%}
#endif //SWIGPYTHON

/// CasADi::MX
#ifdef SWIGOCTAVE
%inline %{

template <>
bool meta< CasADi::MX >::couldbe(const octave_value& p) {
  return (meta< CasADi::MX >::isa(p) || meta< CasADi::Matrix<double> >::couldbe(p) );
}

template <>
int meta< CasADi::MX >::as(const octave_value& p,CasADi::MX &m) {
  NATIVERETURN(CasADi::MX,m)
  if(meta< CasADi::Matrix<double> >::couldbe(p)) {
    CasADi::DMatrix mt;
    bool result=meta< CasADi::Matrix<double> >::as(p,mt);
    if (!result)
      return false;
    m = CasADi::MX(mt);
    return true;
  }
  return false;
}
%}
#endif //SWIGPOCTAVE

%inline %{
template<> char meta< std::vector< CasADi::MX > >::expected_message[] = "Expecting sequence(MX, number)";
%}

/// std::vector< CasADi::MX >
#ifdef SWIGPYTHON
%inline %{
template <>
int meta< std::vector< CasADi::MX > >::as(PyObject * p,std::vector< CasADi::MX > &m) {
  NATIVERETURN(std::vector< CasADi::MX >,m)
  if( PyIsSequence(p) ) {
    return meta< CasADi::MX >::as_vector(p,m);
  } else if (PyDict_Check(p)) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    int num=0;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (PyInt_Check(key)) {
        if (PyInt_AsLong(key)+1>num)
          num=PyInt_AsLong(key)+1;
      } else if (PyString_Check(key) && strcmp(PyString_AsString(key),"NUM")==0) {
        if (PyInt_Check(value))
          num=PyInt_AsLong(value);
      } else {
        return false;
      }
    }
    m.resize(num);
    pos=0;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (PyInt_Check(key)) {
        bool result=meta< CasADi::MX >::as(value,m[PyInt_AsLong(key)]);
        if (!result)
          return false;
      }
    }
    return true;
  } else {
    return false;
  }
return true;
}

template <>
bool meta< std::vector< CasADi::MX > >::couldbe(PyObject * p) {
  if(meta< CasADi::MX >::couldbe_sequence(p)) {
    return true;
  } else if (PyDict_Check(p)) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (!((PyInt_Check(key) || (PyString_Check(key) && strcmp(PyString_AsString(key),"NUM")==0)) && meta< CasADi::MX >::couldbe(value)))
        return false;
    }
    return true;
  }
  return meta< std::vector<CasADi::MX> >::isa(p);
}

%}
#endif //SWIGPYTHON

/// std::vector< CasADi::MX >
#ifdef SWIGOCTAVE
%inline %{
template <>
int meta< std::vector< CasADi::MX > >::as(const octave_value& p,std::vector< CasADi::MX > &m) {
  int nrow = p.rows();
  int ncol = p.columns();
  if(nrow != 1) return false;
  m.resize(ncol);
  
  for(int i=0; i<ncol; ++i){
    // Get the octave object
    const octave_value& obj_i = p.cell_value()(i);
    if (!(obj_i.is_real_matrix() && obj_i.is_empty())) {
      bool ret = meta< CasADi::MX >::as(obj_i,m[i]);
      if(!ret) return false;
    }
  }
  return true;
}

template <> bool meta< std::vector< CasADi::MX > >::couldbe(const octave_value& p) {return p.is_cell();}

%}
#endif //SWIGOCTAVE

%template(sparsity_vector) std::vector<CasADi::CRSSparsity>;

#ifdef SWIGPYTHON
%meta_pair(CasADi::MX, std::vector< CasADi::MX >)
%typemap(out) std::pair< CasADi::MX, std::vector< CasADi::MX >  > {
    bool ret = meta< std::pair< CasADi::MX, std::vector< CasADi::MX >  > >::toPython($1,$result);
    if (!ret) SWIG_exception_fail(SWIG_TypeError,"Could not convert to (MX,std::vector<MX>)");
}
#endif //SWIGPYTHON

%extend CasADi::MX{
  %python_matrix_helpers
};

#ifdef SWIGPYTHON
%extend CasADi::MX{ 
  #binopsFull(double b,CasADi::MX,,CasADi::MX)
  #binopsFull(const CasADi::Matrix<double>& b,CasADi::MX,,CasADi::MX)
};
#endif // SWIGPYTHON

%my_generic_const_typemap(PRECEDENCE_MXVector,std::vector< CasADi::MX >);
#ifdef SWIGPYTHON
%my_generic_const_typemap(PRECEDENCE_DMatrixVector,std::vector< CasADi::Matrix<double> >);
#endif // SWIGPYTHON

%include "casadi/mx/mx_tools.hpp"

