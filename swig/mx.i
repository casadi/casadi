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
  NATIVERETURN(CasADi::Matrix<double>,m)
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

#ifdef SWIGOCTAVE
%meta_vector(CasADi::MX);
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
  %python_matrix_helpers(CasADi::MX)
  #ifdef SWIGPYTHON
  %pythoncode %{
  __array_priority__ = 1002.0
  
  
  def __array_wrap__(self,out_arr,context=None):
    if context is None:
      return out_arr
    name = context[0].__name__
    args = list(context[1])
    
    if "vectorized" in name:
      name = name[:-len(" (vectorized)")]

    conversion = {"multiply": "mul", "divide": "div", "subtract":"sub","power":"pow"}
    if name in conversion:
      name = conversion[name]
    if len(context[1])==2 and context[1][1] is self:
      name = 'r' + name
      args.reverse()
    if not(hasattr(self,name)):
      name = '__' + name + '__'
    fun=getattr(self, name)
    return fun(*args[1:])
      
  def __array__(self,*args,**kwargs):
    import numpy as n
    if len(args) > 1 and isinstance(args[1],tuple) and isinstance(args[1][0],n.ufunc):
      return n.array([1])
    else:
      raise Exception("MX cannot be converted to an array. MX.__array__ purely exists to allow ufunc/numpy goodies")
  %}
  #endif //SWIGPYTHON
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

