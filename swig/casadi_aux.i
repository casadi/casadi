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
%include "casadi/generic_type.hpp"

#ifdef SWIGPYTHON

%inline %{

bool asPyObject(const CasADi::GenericType &a,PyObject * & p) {
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

bool asPyObject(const CasADi::Dictionary &a,PyObject * & p) {
  p = PyDict_New();
  CasADi::Dictionary::const_iterator end = a.end(); 
  for (CasADi::Dictionary::const_iterator it = a.begin(); it != end; ++it)
  {
    PyObject * e;
    bool ret=asPyObject(it->second,e);
    if (!ret) {
      Py_DECREF(p);
      return false;
    }
    PyDict_SetItemString(p,(it->first).c_str(),e);
    Py_DECREF(e);
  }
  return true;
}

bool asGenericType(PyObject *p, CasADi::GenericType &a) {
  if (PyBool_Check(p)) {
    a=CasADi::GenericType((bool) PyInt_AsLong(p));
  } else if (PyInt_Check(p)) {
    a=CasADi::GenericType((int) PyInt_AsLong(p));
  } else if (PyFloat_Check(p)) {
    a=CasADi::GenericType(PyFloat_AsDouble(p));
  } else if (PyString_Check(p)) {
    a=CasADi::GenericType(std::string(PyString_AsString(p)));
  } else {
    return false;
  }
  return true;
}

bool asDictionary(PyObject *p, CasADi::Dictionary &a) {
  if (!PyDict_Check(p))
    return false;
  PyObject *key, *value;
  Py_ssize_t pos = 0;
  CasADi::GenericType gt;
  while (PyDict_Next(p, &pos, &key, &value)) {
    if (!PyString_Check(key))
      return false;
    bool ret=asGenericType(value,gt);
    if (!ret)
      return false;
    a[std::string(PyString_AsString(key))] = gt;
  }

  return true;
}

%}

namespace CasADi {
%typemap(out) GenericType {
bool ret=asPyObject($1,$result);
if (!ret) {
  SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
}

%typemap(out) const Dictionary&  {
bool ret=asPyObject(*$1,$result);
if (!ret) {
  SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
}

%typemap(in) const Dictionary& (Dictionary d)  {
bool ret=asDictionary($input,d);
$1 = &d;
if (!ret) {
  SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
}

%typemap(in) const GenericType& (GenericType d)  {
bool ret=asGenericType($input,d);
$1 = &d;
if (!ret) {
  SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
}


%typemap(typecheck) const Dictionary& {
if (PyDict_Check($input)) {
  $1=1;
} else {
  $1=0;
}
}


} // namespace CasADi
#endif // SWIGPYTHON

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

