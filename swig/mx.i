%{
#include "casadi/mx/mx.hpp"
#include "casadi/mx/mx_tools.hpp"


%}

%include "casadi/mx/mx.hpp"

%inline %{
template<> swig_type_info** meta< CasADi::MX >::name = &SWIGTYPE_p_CasADi__MX;
template<> swig_type_info** meta< std::vector< CasADi::MX> >::name = &SWIGTYPE_p_std__vectorT_CasADi__MX_std__allocatorT_CasADi__MX_t_t;
template<> swig_type_info** meta< std::vector< CasADi::Matrix<double> > >::name = &SWIGTYPE_p_std__vectorT_CasADi__MatrixT_double_t_std__allocatorT_CasADi__MatrixT_double_t_t_t;
%}


/// std::vector< CasADi::Matrix<double> >
#ifdef SWIGPYTHON
%inline %{
template<> char meta< std::vector< CasADi::Matrix<double> > >::expected_message[] = "Expecting sequence(DMatrix)";

template <>
bool meta< std::vector< CasADi::Matrix<double> > >::couldbe(PyObject * p) {
  return meta< std::vector<  CasADi::Matrix<double> > >::isa(p) || meta< CasADi::Matrix<double> >::couldbe_sequence(p);
}

template <>
int meta< std::vector<  CasADi::Matrix<double> > >::as(PyObject * p,std::vector< CasADi::Matrix<double> > &m) {
  return meta< CasADi::Matrix<double> >::as_vector(p,m);
}

%}
#endif

/// CasADi::MX
#ifdef SWIGPYTHON
%inline %{
template<> char meta< CasADi::MX >::expected_message[] = "Expecting (MX, number)";

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


/// std::vector< CasADi::MX >
#ifdef SWIGPYTHON
%inline %{
template<> char meta< std::vector< CasADi::MX > >::expected_message[] = "Expecting sequence(MX, number)";

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


#ifdef SWIGPYTHON
%my_generic_const_typemap(PRECEDENCE_MXVector,std::vector< CasADi::MX >);
%my_generic_const_typemap(PRECEDENCE_DMatrixVector,std::vector< CasADi::Matrix<double> >);
#endif //SWIGPYTHON


%include "casadi/mx/mx_tools.hpp"

#ifdef SWIGPYTHON
namespace CasADi{
  %extend MX{
  
    %python_matrix_helpers

  };
} // namespace CasADi
#endif // SWIGPYTHON

