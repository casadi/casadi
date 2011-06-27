%{
#include <sstream>
#include "casadi/stl_vector_tools.hpp"
#include "casadi/printable_object.hpp"
#include "casadi/shared_object.hpp"
#include "casadi/generic_type.hpp"
#include "casadi/options_functionality.hpp"
%}


%include "casadi/printable_object.hpp"
%include "casadi/shared_object.hpp"


%inline %{
template<> swig_type_info** meta< CasADi::GenericType >::name = &SWIGTYPE_p_CasADi__GenericType;
template<> swig_type_info** meta< CasADi::Dictionary >::name = &SWIGTYPE_p_std__mapT_std__string_CasADi__GenericType_t;
%}

%{
template<> char meta< CasADi::GenericType >::expected_message[] = "Expecting number, string, vector(number)";
%}

/// CasADi::GenericType
#ifdef SWIGPYTHON
%inline %{

template <>
int meta< CasADi::GenericType >::as(PyObject * p,CasADi::GenericType &s) {
  NATIVERETURN(CasADi::GenericType, s)
  if (PyBool_Check(p)) {
    std::cout << "I am supposed to be a bool" << std::endl;
    s=CasADi::GenericType((bool) PyInt_AsLong(p));
  } else if (PyInt_Check(p)) {
    s=CasADi::GenericType((int) PyInt_AsLong(p));
  } else if (PyFloat_Check(p)) {
    s=CasADi::GenericType(PyFloat_AsDouble(p));
  } else if (PyString_Check(p)) {
    s=CasADi::GenericType(std::string(PyString_AsString(p)));
  } else if (meta< std::vector<int> >::couldbe(p)) {
    std::vector<int> temp;
    int ret = meta< std::vector<int> >::as(p,temp); 
    if (!ret) return false;
    s = CasADi::GenericType(temp);
  } else if (meta< std::vector<double> >::couldbe(p)) {
    std::vector<double> temp;
    int ret = meta< std::vector<double> >::as(p,temp); 
    if (!ret) return false;
    s = CasADi::GenericType(temp);
  } else {
    return false;
  }
  return true;
}

template <>
bool meta< CasADi::GenericType >::couldbe(PyObject * p) {
  return meta< CasADi::GenericType >::isa(p) || PyBool_Check(p) ||  PyInt_Check(p) || PyFloat_Check(p) || PyString_Check(p) || meta< std::vector<double> >::couldbe(p);
}

template <>
bool meta< CasADi::GenericType >::toPython(CasADi::GenericType &a, PyObject *&p) {
  if (a.isBool()) {
    p=PyBool_FromLong(a.toBool());
  } else if (a.isInt()) {
    p=PyInt_FromLong(a.toInt());
  } else if (a.isDouble()) {
    p=PyFloat_FromDouble(a.toDouble());
  } else if (a.isString()) {
    p=PyString_FromString(a.toString().c_str());
  } else if (a.isIntVector()) {
    p = swig::from(a.toIntVector());
  } else if (a.isDoubleVector()) {
    p = swig::from( a.toDoubleVector());
  } else {
    return false;
  }
  return true;
}
%}
#endif //SWIGPYTHON

/// CasADi::GenericType
#ifdef SWIGOCTAVE
%inline %{

template <>
int meta< CasADi::GenericType >::as(const octave_value& p,CasADi::GenericType &s) {
  NATIVERETURN(CasADi::GenericType, s)
  if (p.is_real_scalar()) {
    s=CasADi::GenericType(p.double_value());
  } else if (meta< std::vector<int> >::couldbe(p)) {
    std::vector<int> temp;
    int ret = meta< std::vector<int> >::as(p,temp); 
    if (!ret) return false;
    s = CasADi::GenericType(temp);
  } else if (meta< std::vector<double> >::couldbe(p)) {
    std::vector<double> temp;
    int ret = meta< std::vector<double> >::as(p,temp); 
    if (!ret) return false;
    s = CasADi::GenericType(temp);
  } else if (p.is_string()) {
    s = CasADi::GenericType(p.string_value());
  } else {
    return false;
  }
  return true;
}

template <>
bool meta< CasADi::GenericType >::couldbe(const octave_value& p) {
  return p.is_real_scalar() || meta< std::vector<int> >::couldbe(p) || meta< std::vector<double> >::couldbe(p) || p.is_string() ;
}

%}
#endif //SWIGOCTAVE

%{
template<> char meta< CasADi::Dictionary >::expected_message[] = "Expecting dictionary of GenericTypes";
%}

/// CasADi::Dictionary
#ifdef SWIGPYTHON
%inline %{

template <>
int meta< CasADi::Dictionary >::as(PyObject * p,CasADi::Dictionary &s) {
  NATIVERETURN(CasADi::Dictionary, s)
  if (!PyDict_Check(p))
    return false;
  PyObject *key, *value;
  Py_ssize_t pos = 0;
  CasADi::GenericType gt;
  while (PyDict_Next(p, &pos, &key, &value)) {
    if (!PyString_Check(key))
      return false;
    bool ret=meta< CasADi::GenericType >::as(value,gt);
    if (!ret)
      return false;
    s[std::string(PyString_AsString(key))] = gt;
  }

  return true;
}

template <>
bool meta< CasADi::Dictionary >::couldbe(PyObject * p) {
  return PyDict_Check(p);
}

template <>
bool meta< CasADi::Dictionary >::toPython(CasADi::Dictionary &a, PyObject *&p) {
  p = PyDict_New();
  //CasADi::Dictionary::const_iterator end = a.end(); 
  CasADi::Dictionary::iterator end = a.end();
  for (CasADi::Dictionary::iterator it = a.begin(); it != end; ++it)
  {
    PyObject * e;
    bool ret=meta< CasADi::GenericType >::toPython(it->second,e);
    if (!ret) {
      Py_DECREF(p);
      return false;
    }
    PyDict_SetItemString(p,(it->first).c_str(),e);
    Py_DECREF(e);
  }
  return true;
}
%}
#endif //SWIGPYTHON


#ifdef SWIGPYTHON

namespace CasADi {
%typemap(out) GenericType {
bool ret=meta<  CasADi::GenericType >::toPython($1,$result);
if (!ret) {
  SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
}

%typemap(out) const Dictionary&  {
bool ret=meta<  CasADi::Dictionary >::toPython(*$1,$result);
if (!ret) {
  SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
}

} // namespace CasADi
#endif // SWIGPYTHON

%my_generic_const_typemap(PRECEDENCE_GENERICTYPE,CasADi::GenericType)
#ifdef SWIGPYTHON
%my_generic_const_typemap(PRECEDENCE_DICTIONARY ,CasADi::Dictionary)
#endif

#ifdef SWIGOCTAVE
%my_generic_const_typemap(PRECEDENCE_DVector,std::vector<double>);
%my_generic_const_typemap(PRECEDENCE_IVector,std::vector<int>);
#endif //SWIGOCTAVE

%include "casadi/generic_type.hpp"
%include "casadi/options_functionality.hpp"

// Exceptions handling
%include "exception.i"
%exception {
try {
  $action
  } catch (const std::exception& e) {
  SWIG_exception(SWIG_RuntimeError, e.what());
  } catch (const char* e) { // depreciated!!
    SWIG_exception(SWIG_RuntimeError, e);
  }
}

namespace CasADi {
  %extend OptionsFunctionality {
    void setOption(const std::string &name, const std::string& val){$self->setOption(name,val);} 
    void setOption(const std::string &name, const std::vector<int>& val){$self->setOption(name,val);} 
    void setOption(const std::string &name, const std::vector<double>& val){$self->setOption(name,val);} 
    void setOption(const std::string &name, double val){$self->setOption(name,val);}
    void setOption(const std::string &name, int val){$self->setOption(name,val);} 
    void setOption(const std::string &name, bool val){$self->setOption(name,val);}  
  }
} // namespace CasADi

