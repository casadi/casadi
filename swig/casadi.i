/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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



%module(package="casadi",directors=1) casadi

#ifdef SWIGJAVA
%begin %{
#define SWIG_JAVA_GUEST_EXCEPTION "org/casadi/CasadiException"
%}
%include <javaguest.swg>
%include <javacleaner.swg>
// Java has single inheritance: the empty mixin bases become interfaces, so each
// class extends its first real base (DM extends GenDM, MX extends GenMX)
%include <swiginterface.i>
%interface_impl(casadi::PrintableCommon);
%interface_impl(casadi::GenericExpressionCommon);
%interface_impl(casadi::GenericMatrixCommon);
%interface_impl(casadi::SparsityInterfaceCommon);
%interface_impl(casadi::MatrixCommon);
// toString() is str(), as Python's __str__ (an interface cannot override Object methods)
%typemap(javacode) casadi::Matrix<double>, casadi::Matrix<casadi::SXElem>, casadi::Sparsity,
  casadi::Function, casadi::Slice, casadi::DaeBuilder, casadi::Importer, casadi::Linsol,
  casadi::NlpBuilder, casadi::Opti, casadi::OptiSol, casadi::XmlFile, casadi::Resource, casadi::GraphBuilder,
  casadi::ModelicaParser %{
  @Override
  public String toString() {
    return str();
  }
%}
// f(a, b, ...) and f(map) as Python's __call__; Kotlin spells it f(a, b)
%extend casadi::Function {
  %proxycode %{
  /** One output is returned as is, several as a List; numbers and double arrays become DM. */
  public Object invoke(Object... args) {
    java.util.List<GenericMatrix> in = new java.util.ArrayList<>();
    for (Object a : args) in.add(swigMatrix(a));
    java.util.List<GenericMatrix> r = call(in);
    return r.size() == 1 ? r.get(0) : r;
  }

  /** Call with inputs by name; returns the outputs by name. */
  public java.util.Map<String, GenericMatrix> invoke(java.util.Map<String, ?> args) {
    java.util.Map<String, GenericMatrix> in = new java.util.LinkedHashMap<>();
    for (java.util.Map.Entry<String, ?> e : args.entrySet()) in.put(e.getKey(), swigMatrix(e.getValue()));
    return call(in);
  }

  private static GenericMatrix swigMatrix(Object x) {
    if (x instanceof GenericMatrix) return (GenericMatrix) x;
    if (x instanceof Number) return new DM(((Number) x).doubleValue());
    if (x instanceof double[]) return new DM((double[]) x);
    if (x instanceof double[][]) return new DM((double[][]) x);
    throw new IllegalArgumentException("Cannot convert " + (x == null ? "null" : x.getClass().getName()) + " to a matrix");
  }
  %}
}
// MX and GenericType have str() from SharedObject, which Java does not inherit
%extend casadi::MX { std::string toString() const { return $self->get_str(); } }
%extend casadi::GenericType { std::string toString() const { return $self->get_str(); } }
#endif // SWIGJAVA

#ifdef SWIGPYTHON
// Collect multiple outputs in a tuple (swig/swig#1933)
#undef SWIG_AppendOutput
#if SWIG_VERSION >= 0x040300
#define SWIG_AppendOutput(result, obj)  SWIG_Python_AppendOutput_Tuple(result, obj, $isvoid)
#else
// $isvoid is SWIG >= 4.3; 1 reproduces the older SWIG_Python_AppendOutput
#define SWIG_AppendOutput(result, obj)  SWIG_Python_AppendOutput_Tuple(result, obj, 1)
#endif

%{
SWIGINTERN PyObject*
SWIG_Python_AppendOutput_Tuple(PyObject* result, PyObject* obj, int is_void) {
  PyObject*   o2;
  PyObject*   o3;
  if (!result) {
    result = obj;
  } else if (result == Py_None && is_void) {
    Py_DECREF(result);
    result = obj;
  } else {
    if (!PyTuple_Check(result)) {
      o2 = result;
      result = PyTuple_New(1);
      PyTuple_SET_ITEM(result, 0, o2);
    }
    o3 = PyTuple_New(1);
    PyTuple_SET_ITEM(o3, 0, obj);
    o2 = result;
    result = PySequence_Concat(o2, o3);
    Py_DECREF(o2);
    Py_DECREF(o3);
  }
  return result;
}
%}
#endif // SWIGPYTHON

#ifdef CASADI_WITH_COPYSIGN_UNDEF
%{
#ifdef copysign
#undef copysign
#endif
%}
#endif // CASADI_WITH_COPYSIGN_UNDEF

// Include all public CasADi C++
%{
#include <casadi/casadi.hpp>
#include <casadi/core/casadi_interrupt.hpp>
%}

// casadi_int type
%include <casadi/core/casadi_types.hpp>

// No-op %stub_* macros without -stubs
#ifndef SWIG_STUBS_ENABLED
%define %stub_method(NAME, RET, ARGS...) %enddef
%define %stub_method0(NAME, RET) %enddef
%define %stub_overload_method(NAME, RET, ARGS...) %enddef
%define %stub_overload_method0(NAME, RET) %enddef
%define %stub_overload_method_selftyped(NAME, RET, SELFTYPE, ARGS...) %enddef
%define %stub_func(NAME, RET, ARGS...) %enddef
%define %stub_overload_func(NAME, RET, ARGS...) %enddef
%define %stub_attr(NAME, TYPE...) %enddef
%define %stub_const(NAME, TYPE...) %enddef
%define %stub_class_begin(NAME) %enddef
%define %stub_alias_in(NAME, TYPES...) %enddef
%define %stub_alias_out(NAME, TYPES...) %enddef
%define %stub_alias_in_key(KEY, TYPES...) %enddef
%define %stub_alias_out_key(KEY, TYPES...) %enddef
%define %ts_alias_in(NAME) %enddef
%define %ts_alias_out(NAME) %enddef
#endif

  /// Data structure in the target language holding data
#ifdef SWIGPYTHON
#define GUESTOBJECT PyObject
#define CASADI_OK(x) casadi_ok(x)
#define CASADI_NOT_NONE(p) p != Py_None
#elif defined(SWIGMATLAB)
#define GUESTOBJECT mxArray
#elif defined(SWIGWASMJS)
#define GUESTOBJECT SWIG_WASMJS_Object
#elif defined(SWIGJULIA)
#define GUESTOBJECT jl_value_t
#define CASADI_NOT_NONE(p) !p || !jl_is_nothing(p)
#elif defined(SWIGJAVA)
#define GUESTOBJECT _jobject
#define CASADI_NOT_NONE(p) p
#else
#define GUESTOBJECT void
#endif
#ifndef SWIGPYTHON
#define CASADI_OK(x) x
#endif

// Define printing routine

#ifdef SWIGPYTHON

#ifdef CASADI_WITH_PYTHON_GIL_RELEASE
%{
  // This .cxx was swig-compiled with WITH_PYTHON_GIL_RELEASE option
  #define CASADI_WITH_PYTHON_GIL_RELEASE
%}
#else //CASADI_WITH_PYTHON_GIL_RELEASE
%{
  // This .cxx was swig-compiled without WITH_PYTHON_GIL_RELEASE option
  #undef CASADI_WITH_PYTHON_GIL_RELEASE
%}
#endif //CASADI_WITH_PYTHON_GIL_RELEASE

%ignore CASADI_SWIG_FLAGS;
%include "swig_config.h"

%{
  namespace casadi {

    // Redirect printout
    static void pythonlogger(const char* s, std::streamsize num, bool error) {
#ifndef CASADI_WITH_PYTHON_GIL_RELEASE
      if (!casadi::InterruptHandler::is_main_thread()) {
        casadi::Logger::writeDefault(s, num, error);
        return;
      }
#endif // CASADI_WITH_PYTHON_GIL_RELEASE
      int n = num;
#ifdef CASADI_WITH_PYTHON_GIL_RELEASE
      SWIG_PYTHON_THREAD_BEGIN_BLOCK;
#endif // CASADI_WITH_PYTHON_GIL_RELEASE
      while (n>0) {
        if (error) {
          PySys_WriteStderr("%.*s", std::min(n, 1000), s);
        } else {
          PySys_WriteStdout("%.*s", std::min(n, 1000), s);
        }
        n -= 1000;
        s += 1000;
      }
#ifdef CASADI_WITH_PYTHON_GIL_RELEASE
      SWIG_PYTHON_THREAD_END_BLOCK;
#endif // CASADI_WITH_PYTHON_GIL_RELEASE
    }

    static bool pythoncheckinterrupted() {
      if (!casadi::InterruptHandler::is_main_thread()) return false;
#ifdef CASADI_WITH_PYTHON_GIL_RELEASE
      SWIG_PYTHON_THREAD_BEGIN_BLOCK;
#endif // CASADI_WITH_PYTHON_GIL_RELEASE
      return PyErr_CheckSignals();
      // SWIG_PYTHON_THREAD_END_BLOCK is not needed, destructor will release GIL
    }

    // Decodes to "" rather than dereferencing NULL when str_py is absent or not a str
    std::string python_string_to_std_string(PyObject *str_py) {
      if (!str_py) return "";
#if SWIG_VERSION < 0x040200
      const char *chars = SWIG_Python_str_AsChar(str_py);
      std::string str(chars ? chars : "");
      SWIG_Python_str_DelForPy3(chars);
#else
      PyObject *bytes = NULL;
      const char *chars = SWIG_PyUnicode_AsUTF8AndSize(str_py, NULL, &bytes);
      std::string str(chars ? chars : "");
      Py_XDECREF(bytes);
#endif
      if (!chars) PyErr_Clear();
      return str;
    }

    void handle_director_exception() {
	    std::string msg = "Exception in SWIG director ";
      // Note: CASADI_WITH_PYTHON_GIL_RELEASE case has SWIG_PYTHON_THREAD_BEGIN_BLOCK in the caller
#ifndef CASADI_WITH_PYTHON_GIL_RELEASE
      SWIG_PYTHON_THREAD_BEGIN_BLOCK;
#endif // CASADI_WITH_PYTHON_GIL_RELEASE
      if (PyErr_ExceptionMatches(PyExc_KeyboardInterrupt)) {
        PyErr_Clear();
#ifndef CASADI_WITH_PYTHON_GIL_RELEASE
        SWIG_PYTHON_THREAD_END_BLOCK;
#endif // CASADI_WITH_PYTHON_GIL_RELEASE
        throw casadi::KeyboardInterruptException();
      }
      PyObject *ptype, *pvalue, *ptraceback;
      PyErr_Fetch(&ptype, &pvalue, &ptraceback);
      PyObject* msg_py = PyObject_Str(pvalue);
      msg = python_string_to_std_string(msg_py);
      Py_XDECREF(msg_py);
      PyErr_Restore(ptype, pvalue, ptraceback);
      PyErr_Print();
#ifndef CASADI_WITH_PYTHON_GIL_RELEASE
      SWIG_PYTHON_THREAD_END_BLOCK;
#endif // CASADI_WITH_PYTHON_GIL_RELEASE
      casadi_error(msg.c_str());
	  }
  }

%}
%init %{
  // Set logger functions
  casadi::Logger::writeFun = casadi::pythonlogger;

  // @jgillis: please document
  casadi::InterruptHandler::checkInterrupted = casadi::pythoncheckinterrupted;

  casadi::InterruptHandler::is_main_thread();

%}
#elif defined(SWIGMATLAB)
%{
  namespace casadi {
    // Redirect printout to mexPrintf
    static void mexlogger(const char* s, std::streamsize num, bool error) {
      if (!casadi::InterruptHandler::is_main_thread()) {
        casadi::Logger::writeDefault(s, num, error);
        return;
      }
      mexPrintf("%.*s", static_cast<int>(num), s);
    }

#ifdef HAVE_OCTAVE
    // Flush the command window buffer (needed in gui mode)
    static void mexflush(bool error) {
    }
    // Never for Octave
    static bool mexcheckinterrupted() {
      return false;
    }
    void mexclearinterrupted() {

    }
#else
    // Undocumented matlab feature
    extern "C" bool utIsInterruptPending(void);
    extern "C" void utSetInterruptPending(bool);

    static bool mexcheckinterrupted() {
      if (!casadi::InterruptHandler::is_main_thread()) return false;
      return utIsInterruptPending();
    }

    void mexclearinterrupted() {
      utSetInterruptPending(false);
    }

    // Flush the command window buffer (needed in gui mode)
    static void mexflush(bool error) {
      if (!casadi::InterruptHandler::is_main_thread()) {
        casadi::Logger::flushDefault(error);
        return;
      }
      if (!mexcheckinterrupted()) {
        if (mexEvalString("drawnow('update');pause(0.0001);")) {
          utSetInterruptPending(true);
        }
      }
    }
#endif

  }
%}
%init %{
  // Get full path
  mxArray *fullpath, *fullpath_cmd = mxCreateString("fullpath");
  mexCallMATLAB(1, &fullpath, 1, &fullpath_cmd, "mfilename");
  mxDestroyArray(fullpath_cmd);
  std::string path = mxArrayToString(fullpath);
  mxDestroyArray(fullpath);

  // Get file separator
  mxArray *filesep;
  mexCallMATLAB(1, &filesep, 0, 0, "filesep");
  std::string sep = mxArrayToString(filesep);
  mxDestroyArray(filesep);

  // Truncate at separator
  path = path.substr(0, path.rfind(sep));

  // Strip a superfluous +casadi (Octave on Windows)
  if (path.rfind(sep)!=std::string::npos && path.substr(path.rfind(sep)+1)=="+casadi")
    path = path.substr(0, path.rfind(sep));

  // Set library path
  casadi::GlobalOptions::setCasadiPath(path);
  casadi::GlobalOptions::setCasadiIncludePath(path+sep+"include");

  // Matlab is index-one based
  casadi::GlobalOptions::start_index = 1;

  // @jgillis: please document
  mxArray *warning_rhs[] = {mxCreateString("error"),

                            mxCreateString("SWIG:OverloadError")};
  mexCallMATLAB(0, 0, 2, warning_rhs, "warning");
  mxDestroyArray(warning_rhs[0]);
  mxDestroyArray(warning_rhs[1]);


  // Set logger functions
  casadi::Logger::writeFun = casadi::mexlogger;
  casadi::Logger::flush = casadi::mexflush;

  // @jgillis: please document
  casadi::InterruptHandler::checkInterrupted = casadi::mexcheckinterrupted;
  casadi::InterruptHandler::clearInterrupted = casadi::mexclearinterrupted;

  casadi::InterruptHandler::is_main_thread();

%}
#endif

#ifdef SWIGJULIA
%{
  namespace casadi {
    // Printout goes to Julia's stdout/stderr
    static jl_value_t* julia_write = 0;
    static jl_value_t* julia_flush = 0;
    static void julialogger(const char* s, std::streamsize num, bool error) {
      // Threads Julia does not know print directly
      if (!jl_get_pgcstack()) return casadi::Logger::writeDefault(s, num, error);
      SwigJlGcUnsafe enter;
      jl_value_t* str = jl_pchar_to_string(s, num);
      JL_GC_PUSH1(&str);
      if (!jl_call2(julia_write, str, error ? jl_true : jl_false)) jl_exception_clear();
      JL_GC_POP();
    }
    static void juliaflush(bool error) {
      if (!jl_get_pgcstack()) return casadi::Logger::flushDefault(error);
      SwigJlGcUnsafe enter;
      if (!jl_call1(julia_flush, error ? jl_true : jl_false)) jl_exception_clear();
    }
  } // namespace casadi
%}
%wrapper %{
// Set from the Julia module's __init__; the module keeps write and flush alive
extern "C" SWIGEXPORT void casadi_julia_logger(jl_value_t* write, jl_value_t* flush) {
  casadi::julia_write = write;
  casadi::julia_flush = flush;
  casadi::Logger::writeFun = casadi::julialogger;
  casadi::Logger::flush = casadi::juliaflush;
}
%}
#endif // SWIGJULIA

// Spurious warnings, e.g. for vertcat(SXVector)
#pragma SWIG nowarn=509,303,302

#define CASADI_EXPORT

// Incude cmath early on, see #622
%begin %{
#include <cmath>
#ifdef _XOPEN_SOURCE
#undef _XOPEN_SOURCE
#endif
#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif
%}

%ignore *::operator->;

%rename(str) get_str;

#ifdef SWIGPYTHON
%pythoncode %{

import contextlib

# copy / deepcopy support is added to the casadi value types via
# `%extend PrintableCommon` (their common base), not by shadowing the global
# `object` -- see further down in this file.

_swig_repr_default = _swig_repr
def _swig_repr(self):
  if hasattr(self,'repr'):
    return self.repr()
  else:
    return _swig_repr_default(self)

def DM_from_array(m, check_only=True):
  import numpy as np
  if isinstance(m, np.ndarray):
    if len(m.shape)>2:
      return False
    try:
      m = m.astype(float,casting="same_kind",copy=False)
    except:
      return False
    if check_only:
      return True
    else:
      shape = m.shape + (1, 1)
      nrow, ncol = shape[0], shape[1]
      return (nrow,ncol,m.flat)
  return False

def IM_from_array(m, check_only=True):
  import numpy as np
  if isinstance(m, np.ndarray):
    if len(m.shape)>2:
      return False
    try:
      m = m.astype(int,casting="same_kind",copy=False)
    except:
      return False
    if check_only:
      return True
    else:
      shape = m.shape + (1, 1)
      nrow, ncol = shape[0], shape[1]
      return (nrow,ncol,m.flat)
  return False

def SX_from_array(m, check_only=True):
  import numpy as np
  if isinstance(m, np.ndarray):
    if len(m.shape)>2:
      return False
    if m.dtype!=object: return None
    shape = m.shape + (1, 1)
    nrow, ncol = shape[0], shape[1]
    return (nrow,ncol,m.flat)
  return False

def DM_from_csc(m, check_only=True):
  if hasattr(m,"tocsc"):
    m = m.tocsc()
  if m.__class__.__name__ in ["csc_matrix", "csc_array"]:
    if len(m.shape)!=2: return False
    if check_only: return True
    return m.shape + (m.indptr.flat,m.indices.flat,m.data.flat)
  return False

%}
#endif // WITH_SWIGPYTHON


// Doc styles: error, overview, group

%feature("autodoc", "1");

%feature("customdoc", "1");

// Docstrings keep the \lang groups of the target language
#if defined(SWIGMATLAB) || defined(SWIGOCTAVE)
%feature("docstring:lang", "matlab");
#elif defined(SWIGJULIA)
%feature("docstring:lang", "julia");
#elif defined(SWIGWASMJS)
%feature("docstring:lang", "javascript");
#elif defined(SWIGJAVA)
%feature("docstring:lang", "java");
#else
%feature("docstring:lang", "python");
#endif

%feature("customdoc:arg:self", "self");

#if defined(SWIGMATLAB) || defined(SWIGOCTAVE)
  #define MNAME "$NAME"
  %feature("customdoc:proto:constructor", "new_obj = $NAME($in)");
  %feature("customdoc:proto:single_out", "$out = $NAME($in)");
  %feature("customdoc:proto:normal", "[$out] = $NAME($in)");
  %feature("customdoc:main", "    $NAME $brief\n\n$overview\n$main");
#else
  #define MNAME "$name"
  %feature("customdoc:proto:constructor", "$name($in)");
  %feature("customdoc:proto:single_out", "$name($in) -> $out");
  %feature("customdoc:proto:normal", "$name($in) -> ($out)");
  %feature("customdoc:main", "  $brief\n\n::\n\n$overview\n$main");
#endif

%feature("customdoc:arg:normal:style_error", "$type");
%feature("customdoc:arg:only:style_error", "$type");
%feature("customdoc:arg:separator:style_error", ",");
%feature("customdoc:proto:void:style_error", MNAME "($in)");
%feature("customdoc:proto:single_out:style_error", MNAME "($in)");
%feature("customdoc:proto:normal:style_error", MNAME "($in)");
%feature("customdoc:proto:constructor:style_error", MNAME "($in)");

%feature("customdoc:arg:normal", "$type $name");
%feature("customdoc:arg:only", "$type $name");
%feature("customdoc:arg:only:out", "$type");
%feature("customdoc:arg:no_name", "out$ip");
%feature("customdoc:arg:separator", ", ");


%feature("customdoc:proto:void", MNAME "($in)");
%feature("customdoc:proto:single_out:style_group", MNAME "($in)");
%feature("customdoc:proto:normal:style_group", MNAME "($in)");
%feature("customdoc:proto:constructor:style_group", MNAME "($in)");

%feature("customdoc:protoline", "    $proto");
%feature("customdoc:protoline:style_overview", "  $proto");
%feature("customdoc:protoline:nobrief:style_overview", "  $proto");

%feature("customdoc:protoline:style_group", "  $proto");

%feature("customdoc:group", "\n.......\n\n::\n\n$group\n$main\n\n.............\n\n");
#ifdef SWIGJULIA
// Julia docstrings are Markdown: prototypes go in code fences
%feature("customdoc:main", "$brief\n\n```\n$overview```\n$main");
%feature("customdoc:group", "\n---\n\n```\n$group```\n$main\n");
#endif

// append works for strings

%naturalvar;

// Make data members read-only
%immutable;

// Make sure that a copy constructor is created
%copyctor;

#if !defined(SWIGXML) && !defined(SWIGJAVA)
%feature("compactdefaultargs","1");
#endif //SWIGXML

// STL
#ifdef SWIGXML
namespace std {
  template<class T> class vector {};
  template<class A, class B> class pair {};
  template<class A, class B> class map {};
}
#else // SWIGXML
#ifdef SWIGPYTHON
// Overrides of stock fragments; fragments are first-registration-wins, so this precedes stl.i

// New reference to p.item() for a size-1 numpy-like object, else 0 (no error left set)
%fragment("casadi_np_item", "header") {
namespace casadi {
  SWIGINTERN PyObject* np_item(PyObject* p) {
    if (PyErr_Occurred()) PyErr_Clear();
    if (!PyObject_HasAttrString(p, "__array__")) return 0;
    PyObject* sz = PyObject_GetAttrString(p, "size");
    if (!sz) { PyErr_Clear(); return 0; }
    long long n = -1;
    int res = SWIG_AsVal(long long)(sz, &n);
    Py_DECREF(sz);
    if (!SWIG_IsOK(res) || n != 1) return 0;
    PyObject* it = PyObject_CallMethod(p, (char*) "item", 0);
    if (!it) PyErr_Clear();
    return it;
  }
} // namespace casadi
}

%fragment(SWIG_Traits_frag(double), "header", fragment=SWIG_AsVal_frag(double),
          fragment=SWIG_From_frag(double), fragment=SWIG_AsVal_frag(long long),
          fragment="StdTraits", fragment="casadi_np_item") {
namespace swig {
  template <> struct traits< double > {
    typedef value_category category;
    static const char* type_name() { return "double"; }
  };
  template <> struct traits_asval< double > {
    static int asval(PyObject* obj, double* val) {
      int res = SWIG_AsVal(double)(obj, val);
      if (SWIG_IsOK(res)) return res;
      PyObject* it = casadi::np_item(obj);
      if (!it) return SWIG_ERROR;
      res = asval(it, val);
      Py_DECREF(it);
      return res;
    }
  };
  template <> struct traits_from< double > {
    static PyObject* from(const double& val) { return SWIG_From(double)(val); }
  };
}
}

%fragment(SWIG_Traits_frag(long long), "header", fragment=SWIG_AsVal_frag(long long),
          fragment=SWIG_From_frag(long long), fragment="StdTraits", fragment="casadi_np_item") {
namespace swig {
  template <> struct traits< long long > {
    typedef value_category category;
    static const char* type_name() { return "long long"; }
  };
  template <> struct traits_asval< long long > {
    static int asval(PyObject* obj, long long* val) {
      int res = SWIG_AsVal(long long)(obj, val);
      if (SWIG_IsOK(res)) return res;
      PyObject* it = casadi::np_item(obj);
      if (!it) return SWIG_ERROR;
      res = asval(it, val);
      Py_DECREF(it);
      return res;
    }
  };
  template <> struct traits_from< long long > {
    static PyObject* from(const long long& val) { return SWIG_From(long long)(val); }
  };
}
}

// std::pair: a 2-tuple only (stock also takes a 2-element sequence)
%fragment("StdPairTraits", "header", fragment="StdTraits") {
namespace swig {
  template <class T, class U> struct traits_asptr<std::pair<T, U> > {
    static int asptr(PyObject* obj, std::pair<T, U>** val) {
      if (!PyTuple_Check(obj) || PyTuple_Size(obj) != 2) return SWIG_ERROR;
      PyObject* first = PyTuple_GetItem(obj, 0);
      PyObject* second = PyTuple_GetItem(obj, 1);
      if (!val) {
        if (SWIG_IsOK(swig::asval(first, (T*)0)) && SWIG_IsOK(swig::asval(second, (U*)0))) return SWIG_OK;
        if (PyErr_Occurred()) PyErr_Clear();
        return SWIG_ERROR;
      }
      std::pair<T, U>* p = new std::pair<T, U>();
      if (!SWIG_IsOK(swig::asval(first, &p->first)) || !SWIG_IsOK(swig::asval(second, &p->second))) {
        delete p;
        if (PyErr_Occurred()) PyErr_Clear();
        return SWIG_ERROR;
      }
      *val = p;
      return SWIG_NEWOBJ;
    }
  };
  template <class T, class U> struct traits_from<std::pair<T, U> > {
    static PyObject* from(const std::pair<T, U>& val) {
      PyObject* first = swig::from(val.first);
      if (!first) return 0;
      PyObject* second = swig::from(val.second);
      if (!second) { Py_DECREF(first); return 0; }
      PyObject* obj = PyTuple_New(2);
      PyTuple_SetItem(obj, 0, first);
      PyTuple_SetItem(obj, 1, second);
      return obj;
    }
  };
}
}
#elif defined(SWIGMATLAB)
// mxCreateString: '' is returned as 0x0
%fragment(SWIG_From_frag(std::string), "header") {
SWIGINTERNINLINE mxArray* SWIG_From_dec(std::string)(const std::string& s) {
  return mxCreateString(s.c_str());
}
}
#endif
%include "stl.i"
#endif // SWIGXML

%include "doc.i"


// Exceptions handling
%include "exception.i"
#ifdef SWIGJULIA
// casadi runs GC-safe: it may wait on threads that call back into Julia
%exception {
  try {
    SwigJlGcSafe gc_safe;
    $action
  } catch(const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}
#elif defined(SWIGJAVA)
// Every wrapper turns C++ exceptions into CasadiException (javaguest.swg)
%exception {
  $action
}
#else
%exception {
  try {
    $action
  } catch(const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}
#endif

// Python expects __int__ to raise on failure
%exception __int__ {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

// numpy printing only catches TypeError/ValueError (#701)
%exception __bool__ {
 try {
    $action
  } catch (const std::exception& e) {
   SWIG_exception(SWIG_TypeError, e.what());
  }
}

#ifdef SWIGPYTHON
%feature("director:except") {
	if ($error != NULL) casadi::handle_director_exception();
}


%{
#define SWIG_FILE_WITH_INIT
%}

#endif // SWIGPYTHON

#ifdef SWIGJAVA
%fragment(SWIG_From_frag(bool));
%fragment(SWIG_From_frag(long long));
%fragment(SWIG_From_frag(double));
%{
  namespace casadi {
    inline jclass java_class(const char* name) {
      JNIEnv* e = SWIG_JavaGuestEnvChecked();
      jclass c = e->FindClass(name);
      SWIG_JavaGuestCheck(e);
      jclass g = static_cast<jclass>(e->NewGlobalRef(c));
      e->DeleteLocalRef(c);
      return g;
    }
    // casadi printout through System.out and System.err
    static void javalogger(const char* s, std::streamsize num, bool error) {
      JNIEnv* e = SWIG_JavaGuestEnv();
      if (!e || e->ExceptionCheck()) return Logger::writeDefault(s, num, error);
      static jclass sys = java_class("java/lang/System");
      static jclass ps = java_class("java/io/PrintStream");
      static jfieldID out = e->GetStaticFieldID(sys, "out", "Ljava/io/PrintStream;");
      static jfieldID err = e->GetStaticFieldID(sys, "err", "Ljava/io/PrintStream;");
      static jmethodID print = e->GetMethodID(ps, "print", "(Ljava/lang/String;)V");
      jobject stream = e->GetStaticObjectField(sys, error ? err : out);
      jobject str = SWIG_JavaGuestNewString(std::string(s, num));
      e->CallVoidMethod(stream, print, str);
      e->DeleteLocalRef(str);
      e->DeleteLocalRef(stream);
      if (e->ExceptionCheck()) e->ExceptionClear();
    }
    static void javaflush(bool error) {
      JNIEnv* e = SWIG_JavaGuestEnv();
      if (!e || e->ExceptionCheck()) return Logger::flushDefault(error);
      static jclass sys = java_class("java/lang/System");
      static jclass ps = java_class("java/io/PrintStream");
      static jfieldID out = e->GetStaticFieldID(sys, "out", "Ljava/io/PrintStream;");
      static jfieldID err = e->GetStaticFieldID(sys, "err", "Ljava/io/PrintStream;");
      static jmethodID flush = e->GetMethodID(ps, "flush", "()V");
      jobject stream = e->GetStaticObjectField(sys, error ? err : out);
      e->CallVoidMethod(stream, flush);
      e->DeleteLocalRef(stream);
      if (e->ExceptionCheck()) e->ExceptionClear();
    }
  } // namespace casadi
%}
#ifndef _WIN32
%{
#include <dlfcn.h>
%}
#endif
%init %{
  casadi::Logger::writeFun = casadi::javalogger;
  casadi::Logger::flush = casadi::javaflush;
#ifndef _WIN32
  // Plugins next to the loaded libcasadi (a data symbol: function addresses may resolve to a PLT stub)
  {
    Dl_info info;
    if (dladdr(reinterpret_cast<void*>(&casadi::Logger::writeFun), &info) && info.dli_fname) {
      std::string p = info.dli_fname;
      casadi::GlobalOptions::setCasadiPath(p.substr(0, p.find_last_of('/')));
    }
  }
#endif
%}
%pragma(java) jniclasscode=%{
  static {
    System.loadLibrary("casadi_java");
  }
%}
// Free functions live in the module class, like java.lang.Math
%pragma(java) moduleclassmodifiers="public final class"
%pragma(java) modulecode=%{
  private ca() {}
%}
// Merged overloads over the matrix types take <T extends GenericMatrix>
%feature("java:typebound", "GenericMatrix");
%typemap(javainterfaces) casadi::Matrix<double>, casadi::Matrix<casadi::SXElem>, casadi::MX "AutoCloseable, GenericMatrix"
#endif // SWIGJAVA

#ifndef SWIGXML

// Can be overloaded by specifying before importing casadi.i
%fragment("casadi_extra_decl", "header") {}

%fragment("casadi_decl", "header", fragment="casadi_extra_decl", fragment="StdTraits") {
  namespace swig {
    // *m is the caller's object filled in place, or redirected to the wrapped object
    template<class T, class C = typename traits<T>::category, class = void> struct traits_asref {
      static bool asref(GUESTOBJECT *p, T** m) {
        T* r = 0;
        int res = asptr(p, &r);
        if (!SWIG_IsOK(res) || !r) return false;
        if (SWIG_IsNewObj(res)) {
          **m = std::move(*r);
          delete r;
        } else {
          *m = r;
        }
        return true;
      }
    };
    template<class T> struct traits_asref<T, value_category, void> {
      static bool asref(GUESTOBJECT *p, T** m) { return SWIG_IsOK(asval(p, *m)); }
    };
    template<class T> struct traits_asref<T, pointer_category, decltype((void)&traits_asptr<T>::in_place)> {
      static bool asref(GUESTOBJECT *p, T** m) { return traits_asptr<T>::in_place(p, m); }
    };

    // Convert p via *m (see traits_asref); m==0: check only
    template<class T> inline bool asref(GUESTOBJECT *p, T** m) {
      return m ? traits_asref<T>::asref(p, m) : check<T>(p);
    }

    // Convert p into *v; v==0: check only
    template<class T> inline bool assign(GUESTOBJECT *p, T* v) {
      T* q = v;
      if (!asref(p, v ? &q : 0)) return false;
      if (q != v) *v = *q;
      return true;
    }
  } // namespace swig

  namespace casadi {
    /* Check if Null or None */
    bool is_null(GUESTOBJECT *p);

#ifdef SWIGMATLAB
    // Get sparsity pattern
    Sparsity get_sparsity(const mxArray* p);

    // Number of nonzeros
    size_t getNNZ(const mxArray* p);
#endif // SWIGMATLAB

    GUESTOBJECT* full(const DM& m, bool simplify=false) {
#ifdef SWIGPYTHON
      PyObject *p = swig::from(m);
      PyObject *method_name = PyUnicode_FromString("toarray");
      PyObject *cr = PyObject_CallMethodObjArgs(p, method_name, (simplify? Py_True: Py_False), 0);
      Py_DECREF(method_name);
      Py_DECREF(p);
      if (cr) return cr;
      return Py_None;
#elif defined(SWIGMATLAB)
      mxArray *p  = mxCreateDoubleMatrix(m.size1(), m.size2(), mxREAL);
      double* d = static_cast<double*>(mxGetData(p));
      casadi_densify(m.ptr(), m.sparsity(), d, false); // Column-major
      return p;
#else
      return 0;
#endif
    }


    // Convert to a sparse matrix
    GUESTOBJECT* sparse(const DM& m) {
#ifdef SWIGPYTHON
      PyObject *p = swig::from(m);
      PyObject *cr = PyObject_CallMethod(p, (char*) "tocsc", 0);
      Py_DECREF(p);
      if (cr) return cr;
      return Py_None;
#elif defined(SWIGMATLAB)
      mxArray *p  = mxCreateSparse(m.size1(), m.size2(), m.nnz(), mxREAL);
      casadi::casadi_copy(m.ptr(), m.nnz(), static_cast<double*>(mxGetData(p)));
      std::copy(m.colind(), m.colind()+m.size2()+1, mxGetJc(p));
      std::copy(m.row(), m.row()+m.nnz(), mxGetIr(p));
      return p;
#else
      return 0;
#endif

    }

    GUESTOBJECT* full_or_sparse(const DM& m, bool simplify=false) {
      if (m.is_dense()) {
        return full(m, simplify);
      } else {
        GUESTOBJECT* p = sparse(m);
        if (is_null(p)) return swig::from(m);
        return p;
      }
    }
#ifdef SWIGPYTHON


    PyObject* get_Python_helper(const std::string& name) {
%#if PY_VERSION_HEX < 0x03070000
      PyObject* module = PyImport_AddModule("casadi");
%#else
      PyObject* c_name = PyUnicode_FromString("casadi");
      PyObject* module = PyImport_GetModule(c_name);
      Py_DECREF(c_name);
%#endif
      if (!module) {
        if (PyErr_Occurred()) {
          PyErr_Clear();
        }
      }
      PyObject* dict = PyModule_GetDict(module);
      return PyDict_GetItemString(dict, (char*) name.c_str());
    }

    template<class T>
    bool casadi_object_from_fun(GUESTOBJECT *p, T** m, const std::string& fun, const std::function<bool(PyObject*, T**)> & conv) {
      PyObject* dm = get_Python_helper(fun);
      if (!dm) return false;
      PyObject *check_only = m? Py_False : Py_True;
      PyObject *cr = PyObject_CallFunctionObjArgs(dm, p, check_only, NULL);
      if (!cr) return false;
      bool ret;
      // None means not handled (#4216)
      if (cr == Py_None) {
        ret = false;
      } else if (PyBool_Check(cr)) {
        ret = PyObject_IsTrue(cr);
      } else {
        ret = conv(cr, m);
      }
      Py_DECREF(cr);
      return ret;
    }

    bool SX_from_array_conv(GUESTOBJECT *p, casadi::SX** m) {
      std::vector<SXElem> data;
      if (!swig::assign(PyTuple_GetItem(p, 2), &data)) return false;
      casadi_int nrow; swig::assign(PyTuple_GetItem(p, 0), &nrow);
      casadi_int ncol; swig::assign(PyTuple_GetItem(p, 1), &ncol);
      if (m) {
        **m = casadi::SX::zeros(nrow, ncol);
        casadi_densify(get_ptr(data), (**m).sparsity().T(), (**m).ptr(), true);
      }
      return true;
    }

    bool IM_from_array_conv(GUESTOBJECT *p, casadi::IM** m) {
      if (!m) return true;
      std::vector<casadi_int> data;
      if (!swig::assign(PyTuple_GetItem(p, 2), &data)) return false;
      casadi_int nrow; swig::assign(PyTuple_GetItem(p, 0), &nrow);
      casadi_int ncol; swig::assign(PyTuple_GetItem(p, 1), &ncol);
      **m = IM::zeros(nrow, ncol);
      casadi_densify(get_ptr(data), (**m).sparsity().T(), (**m).ptr(), true);
      return true;
    }

    bool DM_from_array_conv(GUESTOBJECT *p, casadi::DM** m) {
      if (!m) return true;
      std::vector<double> data;
      if (!swig::assign(PyTuple_GetItem(p, 2), &data)) return false;
      casadi_int nrow; swig::assign(PyTuple_GetItem(p, 0), &nrow);
      casadi_int ncol; swig::assign(PyTuple_GetItem(p, 1), &ncol);
      **m = DM::zeros(nrow, ncol);
      casadi_densify(get_ptr(data), (**m).sparsity().T(), (**m).ptr(), true);
      return true;
    }

    bool DM_from_csc_conv(GUESTOBJECT *p, casadi::DM** m) {
      std::vector<double> data;
      std::vector<casadi_int> colind, row;
      if (!swig::assign(PyTuple_GetItem(p, 4), &data)) return false;
      if (!swig::assign(PyTuple_GetItem(p, 3), &row)) return false;
      if (!swig::assign(PyTuple_GetItem(p, 2), &colind)) return false;
      casadi_int nrow; swig::assign(PyTuple_GetItem(p, 0), &nrow);
      casadi_int ncol; swig::assign(PyTuple_GetItem(p, 1), &ncol);
      **m = casadi::Matrix<double>(casadi::Sparsity(nrow,ncol,colind,row), data, false);
      return true;
    }

    bool SX_from_array(GUESTOBJECT *p, casadi::SX** m) {
      return casadi_object_from_fun<casadi::SX>(p, m, "SX_from_array", SX_from_array_conv);
    }

    bool IM_from_array(GUESTOBJECT *p, casadi::IM** m) {
      return casadi_object_from_fun<casadi::IM>(p, m, "IM_from_array", IM_from_array_conv);
    }

    bool DM_from_array(GUESTOBJECT *p, casadi::DM** m) {
      return casadi_object_from_fun<casadi::DM>(p, m, "DM_from_array", DM_from_array_conv);
    }

    bool DM_from_csc(GUESTOBJECT *p, casadi::DM** m) {
      return casadi_object_from_fun<casadi::DM>(p, m, "DM_from_csc", DM_from_csc_conv);
    }
#endif
  } // namespace casadi
 }

// Conversions: swig::traits specializations for casadi's typemaps

%fragment("casadi_aux", "header", fragment="casadi_decl", fragment="StdTraits") {
  namespace swig {
#ifdef SWIGPYTHON
    // A failed conversion leaves no Python error set
    inline bool casadi_ok(bool ok) {
      if (!ok && PyErr_Occurred()) PyErr_Clear();
      return ok;
    }

#endif // SWIGPYTHON
    // Proxy pointer (SWIG_OLDOBJ) or a new object filled by conv (SWIG_NEWOBJ); m==0: check only
    template<class M> int casadi_asptr(GUESTOBJECT* p, M** m, swig_type_info* desc, bool (*conv)(GUESTOBJECT*, M*)) {
      if (casadi::is_null(p)) return SWIG_ERROR;
      if (desc) {
        void* vp = 0;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, m ? &vp : 0, desc, 0))) {
          if (m) *m = static_cast<M*>(vp);
          return SWIG_OLDOBJ;
        }
      }
      if (!conv) return SWIG_ERROR;
      if (!m) return CASADI_OK(conv(p, 0)) ? SWIG_OK : SWIG_ERROR;
      M* r = new M();
      if (!CASADI_OK(conv(p, r))) {
        delete r;
        return SWIG_ERROR;
      }
      *m = r;
      return SWIG_NEWOBJ;
    }

    // Same, with asref semantics: *m is redirected to the proxy or conv fills **m in place
    template<class M> bool casadi_in_place(GUESTOBJECT* p, M** m, swig_type_info* desc, bool (*conv)(GUESTOBJECT*, M*)) {
      if (casadi::is_null(p)) return false;
      if (desc) {
        void* vp = 0;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, &vp, desc, 0))) {
          *m = static_cast<M*>(vp);
          return true;
        }
      }
      return conv && CASADI_OK(conv(p, *m));
    }
#ifdef SWIGJULIA

    // asval through in_place: container elements and map values convert into their slot
    template<class M> struct casadi_asval {
      static int asval(GUESTOBJECT* p, M* v) {
        if (!v) return traits_asptr<M>::asptr(p, 0);
        M* vp = v;
        if (!traits_asptr<M>::in_place(p, &vp)) return SWIG_ERROR;
        if (vp != v) *v = *vp;
        return SWIG_OK;
      }
    };
#endif // SWIGJULIA
  } // namespace swig

  namespace casadi {
    bool is_null(GUESTOBJECT *p) {
#ifdef SWIGPYTHON
      return p == Py_None;
#elif defined(SWIGMATLAB)
      if (p == 0) return true;
      return false;
#elif defined(SWIGJULIA)
      return !p || jl_is_nothing(p);
#elif defined(SWIGJAVA)
      return !p;
#else
      if (!p) return true;
      emscripten::val v = SWIG_WASMJS_Borrow(p);
      return v.isNull() || v.isUndefined();
#endif
    }
#ifdef SWIGMATLAB

    Sparsity get_sparsity(const mxArray* p) {
      // Get sparsity pattern
      size_t nrow = mxGetM(p);
      size_t ncol = mxGetN(p);

      if (mxIsSparse(p)) {
        // Sparse storage in MATLAB
        mwIndex *Jc = mxGetJc(p);
        mwIndex *Ir = mxGetIr(p);

        // Store in vectors
        std::vector<casadi_int> colind(ncol+1);
        std::copy(Jc, Jc+colind.size(), colind.begin());
        std::vector<casadi_int> row(colind.back());
        std::copy(Ir, Ir+row.size(), row.begin());

        // Create pattern and return
        return Sparsity(nrow, ncol, colind, row);
      } else {
        return Sparsity::dense(nrow, ncol);
      }
    }

    size_t getNNZ(const mxArray* p) {
      // Dimensions
      size_t nrow = mxGetM(p);
      size_t ncol = mxGetN(p);
      if (mxIsSparse(p)) {
        // Sparse storage in MATLAB
        mwIndex *Jc = mxGetJc(p);
        return Jc[ncol];
      } else {
        return nrow*ncol;
      }
    }
#endif // SWIGMATLAB
  } // namespace casadi
}

// Container policy: casadi's acceptance rules on top of SWIG's generic sequence/map/pair traits
#ifdef SWIGPYTHON
%fragment("casadi_map_traits", "header") {} // Python converts std::map in casadi_containers
#else
%fragment("casadi_map_traits", "header", fragment="StdMapTraits") {}
#endif
%fragment("casadi_containers", "header", fragment="casadi_aux", fragment="StdVectorTraits",
          fragment="casadi_map_traits", fragment="StdPairTraits", fragment=SWIG_Traits_frag(bool),
          fragment=SWIG_Traits_frag(long long), fragment=SWIG_Traits_frag(double),
          fragment=SWIG_Traits_frag(std::string)) {
  namespace swig {
#ifdef SWIGPYTHON
    template<class T, class A> struct traits<std::vector<T, A> > {
      typedef pointer_category category;
      static const char* type_name() {
        static std::string n = "std::vector<" + std::string(swig::type_name<T>()) + " >";
        return n.c_str();
      }
    };
    template<class K, class T, class C, class A> struct traits<std::map<K, T, C, A> > {
      typedef pointer_category category;
      static const char* type_name() {
        static std::string n = "std::map<" + std::string(swig::type_name<K>()) + ","
          + std::string(swig::type_name<T>()) + " >";
        return n.c_str();
      }
    };
    template<class T, class U> struct traits<std::pair<T, U> > {
      typedef pointer_category category;
      static const char* type_name() {
        static std::string n = "std::pair<" + std::string(swig::type_name<T>()) + ","
          + std::string(swig::type_name<U>()) + " >";
        return n.c_str();
      }
    };

    // Not a sequence: str/bytes/dict/set, objects with a conversion hook, non-1-D shape
    inline bool casadi_not_a_sequence(PyObject* p) {
      if (PyList_CheckExact(p) || PyTuple_CheckExact(p)) return false;
      if (PyDict_Check(p) || PyBytes_Check(p) || PySet_Check(p) || PyUnicode_Check(p)) return true;
      if (PyObject_HasAttrString(p, "__MX__") || PyObject_HasAttrString(p, "__SX__")
          || PyObject_HasAttrString(p, "__DM__")) return true;
      if (PyErr_Occurred()) PyErr_Clear();
      if (PyObject_HasAttrString(p, "shape")) {
        PyObject* shape = PyObject_GetAttrString(p, "shape");
        if (!shape) { PyErr_Clear(); return true; }
        bool is_1d = PyTuple_Check(shape) && PyTuple_Size(shape)==1;
        Py_DECREF(shape);
        if (!is_1d) return true;
      }
      return false;
    }

    template<class T, class A> struct IteratorProtocol<std::vector<T, A>, T> {
      static void assign(PyObject* p, std::vector<T, A>* seq) {
        if (casadi_not_a_sequence(p)) throw std::invalid_argument("not a sequence");
        Py_ssize_t n = PyObject_Length(p);
        if (n > 0) seq->reserve(n);
        if (n < 0) PyErr_Clear();
        SwigVar_PyObject it = PyObject_GetIter(p);
        if (!it) return;
        T tmp;
        SwigVar_PyObject item = PyIter_Next(it);
        while (item) {
          T* e = &tmp;
          if (!traits_asref<T>::asref(item, &e)) throw std::invalid_argument("bad type");
          seq->push_back(*e);
          item = PyIter_Next(it);
        }
      }
      static bool check(PyObject* p) {
        if (casadi_not_a_sequence(p)) return false;
        SwigVar_PyObject it = PyObject_GetIter(p);
        if (!it) return casadi_ok(false);
        SwigVar_PyObject item = PyIter_Next(it);
        while (item) {
          if (!swig::check<T>(item)) return false;
          item = PyIter_Next(it);
        }
        return casadi_ok(!PyErr_Occurred());
      }
    };

    // A sequence converts into the caller's vector in place
    template<class T, class A> struct traits_asref<std::vector<T, A>, pointer_category, void> {
      static bool asref(PyObject* p, std::vector<T, A>** m) {
        if (p == Py_None || SWIG_Python_GetSwigThis(p)) return false;
        try {
          (*m)->clear();
          IteratorProtocol<std::vector<T, A>, T>::assign(p, *m);
        } catch (std::exception&) {
          return casadi_ok(false);
        }
        return casadi_ok(!PyErr_Occurred());
      }
    };

    // std::vector maps to a Python list
    template<class T, class A> struct traits_from_stdseq<std::vector<T, A>, T> {
      static PyObject* from(const std::vector<T, A>& seq) {
        PyObject* ret = PyList_New(seq.size());
        if (!ret) return 0;
        for (size_t k=0; k<seq.size(); ++k) {
          PyObject* el = swig::from<T>(seq[k]);
          if (!el) {
            Py_DECREF(ret);
            return 0;
          }
          PyList_SetItem(ret, k, el);
        }
        return ret;
      }
    };

    // std::map<std::string, T> maps to a dict with str keys
    template<class T, class C, class A> struct traits_asptr<std::map<std::string, T, C, A> > {
      typedef std::map<std::string, T, C, A> map_type;
      // m==0: check only
      static bool fill(PyObject* p, map_type* m) {
        if (!PyDict_Check(p)) return false;
        if (m) m->clear();
        PyObject *key, *value;
        Py_ssize_t pos = 0;
        while (PyDict_Next(p, &pos, &key, &value)) {
          if (!PyUnicode_Check(key)) return false;
          if (m) {
            T& v = (*m)[casadi::python_string_to_std_string(key)];
            T* vp = &v;
            if (!traits_asref<T>::asref(value, &vp)) return false;
            if (vp != &v) v = *vp;
          } else if (!swig::check<T>(value)) {
            return false;
          }
        }
        return true;
      }
      static int asptr(PyObject* p, map_type** m) {
        if (!m) return fill(p, 0) ? SWIG_OK : SWIG_ERROR;
        map_type* r = new map_type();
        if (!fill(p, r)) {
          delete r;
          return SWIG_ERROR;
        }
        *m = r;
        return SWIG_NEWOBJ;
      }
      static bool in_place(PyObject* p, map_type** m) {
        return fill(p, *m);
      }
    };
    template<class T, class C, class A> struct traits_from<std::map<std::string, T, C, A> > {
      static PyObject* from(const std::map<std::string, T, C, A>& a) {
        PyObject *p = PyDict_New();
        if (!p) return 0;
        for (typename std::map<std::string, T, C, A>::const_iterator it=a.begin(); it!=a.end(); ++it) {
          PyObject * e = swig::from(it->second);
          if (!e) {
            Py_DECREF(p);
            return 0;
          }
          PyDict_SetItemString(p, it->first.c_str(), e);
          Py_DECREF(e);
        }
        return p;
      }
    };
#elif defined(SWIGMATLAB)
    // Containers without a %template'd traits fragment (e.g. std::vector<Dict>)
    template<class T, class A> struct traits<std::vector<T, A> > {
      typedef pointer_category category;
      static const char* type_name() { return "std::vector"; }
    };
    template<class K, class T, class C, class A> struct traits<std::map<K, T, C, A> > {
      typedef pointer_category category;
      static const char* type_name() { return "std::map"; }
    };
    template<class T, class U> struct traits<std::pair<T, U> > {
      typedef pointer_category category;
      static const char* type_name() { return "std::pair"; }
    };

    inline bool casadi_is_vector(mxArray* p) {
      return mxGetNumberOfDimensions(p)==2 && (mxGetM(p)<=1 || mxGetN(p)<=1);
    }

    // Cells only: native arrays are handled by the specializations below
    template<class Seq> inline int casadi_asptr_cell(mxArray* p, Seq** m) {
      if (!p || !mxIsCell(p)) return SWIG_ERROR;
      return traits_asptr_stdseq<Seq>::asptr(p, m);
    }

    // Double vector or cell
    template<> struct traits_asptr<std::vector<double> > {
      static int asptr(mxArray* p, std::vector<double>** m) {
        if (p && mxIsDouble(p) && casadi_is_vector(p)) {
          const double* d = static_cast<const double*>(mxGetData(p));
          if (m) *m = new std::vector<double>(d, d + mxGetM(p)*mxGetN(p));
          return m ? SWIG_NEWOBJ : SWIG_OK;
        }
        return casadi_asptr_cell(p, m);
      }
    };

    // Integral double vector, non-scalar logical vector or cell
    template<> struct traits_asptr<std::vector<casadi_int> > {
      static int asptr(mxArray* p, std::vector<casadi_int>** m) {
        if (p && mxIsDouble(p) && casadi_is_vector(p)) {
          const double* d = static_cast<const double*>(mxGetData(p));
          size_t n = mxGetM(p)*mxGetN(p);
          for (size_t i=0; i<n; ++i) {
            if (d[i]!=static_cast<casadi_int>(d[i])) return casadi_asptr_cell(p, m);
          }
          if (m) *m = new std::vector<casadi_int>(d, d + n);
          return m ? SWIG_NEWOBJ : SWIG_OK;
        }
        if (p && mxIsLogical(p) && !mxIsLogicalScalar(p) && casadi_is_vector(p)) {
          const mxLogical* d = static_cast<const mxLogical*>(mxGetData(p));
          if (m) *m = new std::vector<casadi_int>(d, d + mxGetM(p)*mxGetN(p));
          return m ? SWIG_NEWOBJ : SWIG_OK;
        }
        return casadi_asptr_cell(p, m);
      }
    };

    // Double vector of zeros and ones, non-scalar logical vector or cell
    template<> struct traits_asptr<std::vector<bool> > {
      static int asptr(mxArray* p, std::vector<bool>** m) {
        if (p && mxIsDouble(p) && casadi_is_vector(p)) {
          const double* d = static_cast<const double*>(mxGetData(p));
          size_t n = mxGetM(p)*mxGetN(p);
          for (size_t i=0; i<n; ++i) {
            if (d[i]!=0 && d[i]!=1) return casadi_asptr_cell(p, m);
          }
          if (m) *m = new std::vector<bool>(d, d + n);
          return m ? SWIG_NEWOBJ : SWIG_OK;
        }
        if (p && mxIsLogical(p) && !mxIsLogicalScalar(p) && casadi_is_vector(p)) {
          const mxLogical* d = static_cast<const mxLogical*>(mxGetData(p));
          if (m) *m = new std::vector<bool>(d, d + mxGetM(p)*mxGetN(p));
          return m ? SWIG_NEWOBJ : SWIG_OK;
        }
        return casadi_asptr_cell(p, m);
      }
    };

    // Returned as a char matrix, one string per row
    template<> struct traits_from<std::vector<std::string> > {
      static mxArray* from(const std::vector<std::string>& v) {
        std::vector<const char*> s(v.size());
        for (size_t i=0; i<s.size(); ++i) s[i] = v[i].c_str();
        return mxCreateCharMatrixFromStrings(s.size(), s.empty() ? 0 : &s[0]);
      }
    };

    // Cell of vectors, else any cell with one row per inner vector
    template<class T> struct traits_asptr<std::vector<std::vector<T> > > {
      typedef std::vector<std::vector<T> > VV;
      static int asptr(mxArray* p, VV** m) {
        int res = traits_asptr_stdseq<VV>::asptr(p, m);
        if (SWIG_IsOK(res) || !p || !mxIsCell(p)) return res;
        size_t nrow = mxGetM(p), ncol = mxGetN(p);
        VV* v = m ? new VV(nrow, std::vector<T>(ncol)) : 0;
        for (size_t i=0; i<nrow*ncol; ++i) {
          mxArray* e = mxGetCell(p, i);
          T t;
          if (!e || !SWIG_IsOK(swig::asval<T>(e, v ? &t : 0))) {
            delete v;
            return SWIG_ERROR;
          }
          if (v) (*v)[i % nrow][i / nrow] = t;
        }
        if (!m) return SWIG_OK;
        *m = v;
        return SWIG_NEWOBJ;
      }
    };

    // (casadi_int, casadi_int) is a 1-by-2 double matrix with integral entries
    template<> struct traits_asptr<std::pair<casadi_int, casadi_int> > {
      static int asptr(mxArray* p, std::pair<casadi_int, casadi_int>** m) {
        if (p && mxIsDouble(p) && mxGetNumberOfDimensions(p)==2 && !mxIsSparse(p)
            && mxGetM(p)==1 && mxGetN(p)==2) {
          const double* d = static_cast<const double*>(mxGetData(p));
          casadi_int first = static_cast<casadi_int>(d[0]);
          casadi_int second = static_cast<casadi_int>(d[1]);
          if (d[0]!=first || d[1]!=second) return SWIG_ERROR;
          if (m) *m = new std::pair<casadi_int, casadi_int>(first, second);
          return m ? SWIG_NEWOBJ : SWIG_OK;
        }
        return SWIG_ERROR;
      }
    };
    template<> struct traits_from<std::pair<casadi_int, casadi_int> > {
      static mxArray* from(const std::pair<casadi_int, casadi_int>& v) {
        mxArray* ret = mxCreateDoubleMatrix(1, 2, mxREAL);
        double* d = static_cast<double*>(mxGetData(ret));
        d[0] = v.first;
        d[1] = v.second;
        return ret;
      }
    };
#elif defined(SWIGJAVA)
    // (casadi_int, casadi_int) is a long[2], or a Map.Entry
    template<> struct traits_asptr<std::pair<casadi_int, casadi_int> > {
      static int asptr(jobject p, std::pair<casadi_int, casadi_int>** m) {
        if (!p || !SWIG_JavaGuestIsLongs(p)) return traits_asptr_stdpair<casadi_int, casadi_int>::asptr(p, m);
        std::vector<long long> v = SWIG_JavaGuestLongs(p);
        if (v.size()!=2) return SWIG_ERROR;
        if (m) *m = new std::pair<casadi_int, casadi_int>(v[0], v[1]);
        return m ? SWIG_NEWOBJ : SWIG_OK;
      }
    };
    template<> struct traits_from<std::pair<casadi_int, casadi_int> > {
      static jobject from(const std::pair<casadi_int, casadi_int>& v) {
        std::vector<long long> r(2);
        r[0] = v.first;
        r[1] = v.second;
        return SWIG_JavaGuestNewLongs(r);
      }
    };
#elif defined(SWIGJULIA)
    // Other AbstractVectors (ranges, views, ...): collect, then convert the Array
    template<class T> struct traits_asptr<std::vector<T, std::allocator<T> > > {
      static int asptr(jl_value_t* p, std::vector<T>** m) {
        jl_value_t* av = p && !jl_is_array(p) ? jl_get_global(jl_base_module, jl_symbol("AbstractVector")) : 0;
        if (!av || !jl_isa(p, av)) return traits_asptr_stdseq<std::vector<T> >::asptr(p, m);
        jl_value_t* c = 0;
        JL_GC_PUSH1(&c);
        c = jl_call1(jl_get_function(jl_base_module, "collect"), p);
        int res = c && jl_is_array(c) ? traits_asptr_stdseq<std::vector<T> >::asptr(c, m) : SWIG_ERROR;
        JL_GC_POP();
        return res;
      }
    };
#elif defined(SWIGWASMJS)
    // Numbers also accept a JS boolean
    template<class T> struct casadi_bool_number {
      static bool asref(EM_VAL p, T** m) {
        if (SWIG_IsOK(asval(p, *m))) return true;
        bool b;
        if (!SWIG_IsOK(asval(p, &b))) return false;
        **m = b;
        return true;
      }
      static bool check(EM_VAL p) {
        return SWIG_IsOK(asval(p, (T*)0)) || SWIG_IsOK(asval(p, (bool*)0));
      }
    };
    template<> struct traits_asref<double, value_category, void> : casadi_bool_number<double> {};
    template<> struct traits_check<double, value_category> : casadi_bool_number<double> {};
    template<> struct traits_asref<long long, value_category, void> : casadi_bool_number<long long> {};
    template<> struct traits_check<long long, value_category> : casadi_bool_number<long long> {};

    // A null or undefined dict is an empty dict
    template<class T, class C, class A> struct traits_asref<std::map<std::string, T, C, A>, pointer_category, void> {
      typedef std::map<std::string, T, C, A> map_type;
      static bool asref(EM_VAL p, map_type** m) {
        if (casadi::is_null(p)) return true;
        map_type* r = 0;
        int res = asptr(p, &r);
        if (!SWIG_IsOK(res) || !r) return false;
        if (SWIG_IsNewObj(res)) {
          **m = std::move(*r);
          delete r;
        } else {
          *m = r;
        }
        return true;
      }
    };
    template<class T, class C, class A> struct traits_check<std::map<std::string, T, C, A>, pointer_category> {
      static bool check(EM_VAL p) {
        return casadi::is_null(p) || SWIG_IsOK(asptr(p, (std::map<std::string, T, C, A>**)0));
      }
    };

    // Array elements convert with the casadi rules; in a check, a wrapped element must be of the element type itself
    template<class T, class A> struct JSArrayProtocol<std::vector<T, A>, T> {
      static bool assign(const emscripten::val& arr, std::vector<T, A>* seq) {
        unsigned n = arr["length"].as<unsigned>();
        seq->reserve(n);
        T tmp;
        for (unsigned i = 0; i < n; ++i) {
          emscripten::val el = arr[i];
          T* e = &tmp;
          if (!traits_asref<T>::asref(el.as_handle(), &e)) return false;
          seq->push_back(*e);
        }
        return true;
      }
      static bool exact(EM_VAL p, pointer_category) {
        T* r = 0;
        int res = asptr(p, &r);
        if (SWIG_IsOK(res) && SWIG_IsNewObj(res)) delete r;
        return SWIG_IsOK(res) && !SWIG_IsNewObj(res);
      }
      static bool exact(EM_VAL, value_category) { return false; }
      static bool check(const emscripten::val& arr) {
        unsigned n = arr["length"].as<unsigned>();
        for (unsigned i = 0; i < n; ++i) {
          emscripten::val el = arr[i];
          std::string ty = el.typeOf().as<std::string>();
          bool tagged = (ty == "object" || ty == "function") && !el.isNull()
            && el["_swig_type"].typeOf().as<std::string>() == "string";
          if (tagged) {
            if (!exact(el.as_handle(), typename traits<T>::category())) return false;
          } else {
            T tmp, *e = &tmp;
            if (!swig::asref(el.as_handle(), &e)) return false;
          }
        }
        return true;
      }
    };
#endif
  } // namespace swig
}

// Element traits, one fragment per type, named after the spelling in the %template() lines

// Wrapped classes without conversions
%define %casadi_class_traits(xName, xType...)
%fragment(SWIG_Traits_frag(xType), "header", fragment="casadi_aux") {
  namespace swig {
    template<> struct traits<xType> {
      typedef pointer_category category;
      static const char* type_name() { return xName; }
    };
    template<> struct traits_asptr<xType> {
      static int asptr(GUESTOBJECT* p, xType** m) {
        return casadi_asptr<xType>(p, m, $descriptor(xType*), 0);
      }
#ifdef SWIGJULIA
      static bool in_place(GUESTOBJECT* p, xType** m) {
        return casadi_in_place<xType>(p, m, $descriptor(xType*), 0);
      }
#endif
    };
#ifdef SWIGJULIA
    template<> struct traits_asval<xType> : casadi_asval<xType> {};
#endif
    template<> struct traits_from<xType> {
      static GUESTOBJECT* from(const xType& a) {
        return SWIG_NewPointerObj(new xType(a), $descriptor(xType*), SWIG_POINTER_OWN);
      }
    };
  } // namespace swig
}
%enddef

%casadi_class_traits("Sparsity", casadi::Sparsity)
%casadi_class_traits("Function", casadi::Function)
#ifdef SWIGWASMJS
// Linsol: only as the element of a std::vector in the serializer API
%casadi_class_traits("Linsol", casadi::Linsol)
#endif
#ifdef SWIGJAVA
// An integer i is the slice i:i+1
%fragment(SWIG_Traits_frag(casadi::Slice), "header", fragment="casadi_aux", fragment=SWIG_Traits_frag(long long)) {
  namespace swig {
    template<> struct traits<casadi::Slice> {
      typedef pointer_category category;
      static const char* type_name() { return "Slice"; }
    };
    template<> struct traits_asptr<casadi::Slice> {
      static bool conv(jobject p, casadi::Slice* m) {
        long long i;
        if (!SWIG_IsOK(swig::asval(p, &i))) return false;
        if (m) {
          m->start = i;
          m->stop = m->start+1;
          if (m->stop==0) m->stop = std::numeric_limits<casadi_int>::max();
        }
        return true;
      }
      static int asptr(jobject p, casadi::Slice** m) {
        return casadi_asptr<casadi::Slice>(p, m, $descriptor(casadi::Slice*), conv);
      }
      static bool in_place(jobject p, casadi::Slice** m) {
        return casadi_in_place<casadi::Slice>(p, m, $descriptor(casadi::Slice*), conv);
      }
    };
    template<> struct traits_from<casadi::Slice> {
      static jobject from(const casadi::Slice& a) {
        return SWIG_NewPointerObj(new casadi::Slice(a), $descriptor(casadi::Slice *), SWIG_POINTER_OWN);
      }
    };
  } // namespace swig
}
#elif !defined(SWIGPYTHON)
%casadi_class_traits("Slice", casadi::Slice)
#else
%fragment(SWIG_Traits_frag(casadi::Slice), "header", fragment="casadi_aux") {
  namespace swig {
    template<> struct traits<casadi::Slice> {
      typedef pointer_category category;
      static const char* type_name() { return "Slice"; }
    };
    template<> struct traits_asptr<casadi::Slice> {
      static bool conv(PyObject* p, casadi::Slice* m) {
        // Python casadi_int
        if (PyLong_Check(p)) {
          if (m) {
            m->start = PyLong_AsLong(p);
            m->stop = m->start+1;
            if (m->stop==0) m->stop = std::numeric_limits<casadi_int>::max();
          }
          return true;
        }
        // Python slice - use Limited API compatible approach
        if (PySlice_Check(p)) {
          Py_ssize_t start, stop, step;
%#if PY_VERSION_HEX >= 0x03060100
          int res = PySlice_Unpack(p, &start, &stop, &step);
%#else
          // Python 2.7 and early Python 3.x use _PySlice_Unpack (private API)
          int res = _PySlice_Unpack(p, &start, &stop, &step);
%#endif
          if (res < 0) {
            return false;  // TypeError already set by PySlice_Unpack
          }

          if (m) {
            // Map PySlice_Unpack sentinels (either extreme) to CasADi's limits
            m->start = (start == PY_SSIZE_T_MIN || start == PY_SSIZE_T_MAX)
                ? std::numeric_limits<casadi_int>::min()
                : static_cast<casadi_int>(start);
            m->stop = (stop == PY_SSIZE_T_MAX || stop == PY_SSIZE_T_MIN)
                ? std::numeric_limits<casadi_int>::max()
                : static_cast<casadi_int>(stop);
            if (step != 1) {
              m->step = static_cast<casadi_int>(step);
            }
          }
          return true;
        }
        return false;
      }
      static int asptr(PyObject* p, casadi::Slice** m) {
        return casadi_asptr<casadi::Slice>(p, m, $descriptor(casadi::Slice*), conv);
      }
      static bool in_place(PyObject* p, casadi::Slice** m) {
        return casadi_in_place<casadi::Slice>(p, m, $descriptor(casadi::Slice*), conv);
      }
    };
    template<> struct traits_from<casadi::Slice> {
      static PyObject* from(const casadi::Slice& a) {
        return SWIG_NewPointerObj(new casadi::Slice(a), $descriptor(casadi::Slice *), SWIG_POINTER_OWN);
      }
    };
  } // namespace swig
}
#endif // SWIGPYTHON

%fragment(SWIG_Traits_frag(casadi::Matrix<double>), "header", fragment="casadi_aux",
          fragment=SWIG_Traits_frag(double), fragment="casadi_containers") {
  namespace swig {
    template<> struct traits<casadi::Matrix<double> > {
      typedef pointer_category category;
      static const char* type_name() { return "DM"; }
    };
    template<> struct traits_asptr<casadi::Matrix<double> > {
#ifdef SWIGJULIA
      // Real Julia array of at most two dimensions to dense DM (vector -> column)
      static bool from_julia_array(jl_value_t* p, casadi::DM* m) {
        using namespace casadi;
        jl_value_t* aa = jl_get_global(jl_base_module, jl_symbol("AbstractArray"));
        if (!aa || !jl_isa(p, aa)) return false;
        jl_value_t *et = 0, *sz = 0, *arr = 0, *atype = 0;
        JL_GC_PUSH4(&et, &sz, &arr, &atype);
        et = jl_call1(jl_get_function(jl_base_module, "eltype"), p);
        sz = jl_call1(jl_get_function(jl_base_module, "size"), p);
        size_t nd = sz && jl_is_tuple(sz) ? jl_nfields(sz) : 3;
        bool ok = et && nd <= 2;
        casadi_int dims[2] = {1, 1};
        for (size_t i = 0; ok && i < nd; ++i) dims[i] = jl_unbox_int64(jl_fieldref(sz, i));
        bool real = ok && jl_is_type(et)
          && jl_subtype(et, jl_get_global(jl_base_module, jl_symbol("Real")));
        if (ok && !real) {
          // `[]` has eltype Any: accept it as an empty DM, reject anything else
          ok = dims[0] * dims[1] == 0;
          if (ok && m) *m = DM(0, 0);
        } else if (ok && m) {
          atype = jl_apply_array_type((jl_value_t*)jl_float64_type, nd);
          arr = jl_call2(jl_get_function(jl_base_module, "convert"), atype, p);
          jl_value_t* ptr = arr ? jl_call1(jl_get_function(jl_base_module, "pointer"), arr) : 0;
          ok = ptr != 0;
          if (ok) {
            const double* d = static_cast<const double*>(jl_unbox_voidpointer(ptr));
            *m = DM::zeros(dims[0], dims[1]);
            std::copy(d, d + dims[0] * dims[1], m->ptr());
          }
        }
        JL_GC_POP();
        return ok;
      }
#elif defined(SWIGWASMJS)
      // A JS number or BigInt as a double
      static bool number(const emscripten::val& el, double* d) {
        std::string ty = el.typeOf().as<std::string>();
        if (ty == "number") *d = el.as<double>();
        else if (ty == "bigint") *d = emscripten::val::global("Number")(el).as<double>();
        else return false;
        return true;
      }
#endif
      static bool conv(GUESTOBJECT* p, casadi::DM* m) {
        using namespace casadi;
        // Object is a sparsity pattern
        {
          void *sp = 0;
          if (SWIG_IsOK(SWIG_ConvertPtr(p, &sp, $descriptor(casadi::Sparsity*), 0))) {
            if (m) *m = DM::ones(*static_cast<Sparsity*>(sp));
            return true;
          }
        }

        // Double scalar
        {
          double tmp;
          if (swig::assign(p, m? &tmp: 0)) {
            if (m) *m = tmp;
            return true;
          }
        }

#ifdef SWIGPYTHON
        // Object has __DM__ method; None means not a DM
        if (PyErr_Occurred()) PyErr_Clear(); // Clear pending exception before type check
        if (PyObject_HasAttrString(p,"__DM__")) {
          char name[] = "__DM__";
          PyObject *cr = PyObject_CallMethod(p, name, 0);
          if (!cr) return false;
          if (cr == Py_None) { Py_DECREF(cr); return false; }
          bool result = swig::assign(cr, m);
          Py_DECREF(cr);
          return result;
        }

        DM* mp = m;
        if (DM_from_array(p, m ? &mp : 0)) return true;

        if (DM_from_csc(p, m ? &mp : 0)) return true;

        {
          std::vector<double> t;
          bool res = swig::assign(p, m ? &t : 0);
          if (m) *m = t.empty() ? casadi::Matrix<double>(0,0) : casadi::Matrix<double>(t);
          return res;
        }
#elif defined(SWIGMATLAB)
        // MATLAB double matrix (sparse or dense)
        if (mxIsDouble(p) && mxGetNumberOfDimensions(p)==2) {
          if (m) {
            *m = casadi::DM(get_sparsity(p));
            double* data = static_cast<double*>(mxGetData(p));
            casadi_copy(data, m->nnz(), m->ptr());
          }
          return true;
        }

        // No match
        return false;
#elif defined(SWIGJULIA)
        if (from_julia_array(p, m)) return true;

        // No match
        return false;
#elif defined(SWIGJAVA)
        // double[] or List of numbers: column; double[][] or nested: one entry per row
        {
          std::vector<double> t;
          if (swig::assign(p, m ? &t : 0)) {
            if (m) *m = t.empty() ? casadi::Matrix<double>(0, 0) : casadi::Matrix<double>(t);
            return true;
          }
        }
        {
          std::vector< std::vector<double> > t;
          if (swig::assign(p, m ? &t : 0)) {
            if (m) *m = casadi::Matrix<double>(t);
            return true;
          }
        }
        return false;
#else
        // JS Array of numbers -> column vector; nested [[a,b],[c,d]] -> dense matrix, row-major in JS
        emscripten::val v = SWIG_WASMJS_Borrow(p);
        if (!v.isArray()) return false;
        unsigned n = v["length"].as<unsigned>();
        bool nested = false;
        if (n > 0) {
          emscripten::val first = v[0];
          nested = first.isArray();
        }
        if (nested) {
          unsigned ncols = 0;
          for (unsigned i = 0; i < n; ++i) {
            emscripten::val row = v[i];
            if (!row.isArray()) return false;
            unsigned w = row["length"].as<unsigned>();
            if (i == 0) ncols = w;
            else if (w != ncols) return false;
          }
          std::vector<double> data((size_t)n * (size_t)ncols);
          for (unsigned i = 0; i < n; ++i) {
            emscripten::val row = v[i];
            for (unsigned j = 0; j < ncols; ++j) {
              emscripten::val el = row[j];
              // Column-major: data[j*n + i]
              if (!number(el, &data[(size_t)j * (size_t)n + (size_t)i])) return false;
            }
          }
          if (m) *m = casadi::Matrix<double>::reshape(casadi::Matrix<double>(data), n, ncols);
          return true;
        }
        std::vector<double> data(n);
        for (unsigned i = 0; i < n; ++i) {
          emscripten::val el = v[i];
          if (!number(el, &data[i])) return false;
        }
        if (m) *m = casadi::Matrix<double>(data);
        return true;
#endif
      }
      static int asptr(GUESTOBJECT* p, casadi::DM** m) {
        return casadi_asptr<casadi::DM>(p, m, $descriptor(casadi::Matrix<double>*), conv);
      }
      static bool in_place(GUESTOBJECT* p, casadi::DM** m) {
        return casadi_in_place<casadi::DM>(p, m, $descriptor(casadi::Matrix<double>*), conv);
      }
    };
#ifdef SWIGJULIA
    template<> struct traits_asval<casadi::Matrix<double> > : casadi_asval<casadi::DM> {};
#endif
    template<> struct traits_from<casadi::Matrix<double> > {
      static GUESTOBJECT* from(const casadi::DM& a) {
        return SWIG_NewPointerObj(new casadi::DM(a), $descriptor(casadi::Matrix<double>*), SWIG_POINTER_OWN);
      }
    };
  } // namespace swig
}

%fragment(SWIG_Traits_frag(casadi::Matrix<casadi_int>), "header", fragment="casadi_aux",
          fragment=SWIG_Traits_frag(long long), fragment=SWIG_Traits_frag(casadi::Matrix<double>)) {
  namespace swig {
    template<> struct traits<casadi::Matrix<casadi_int> > {
      typedef pointer_category category;
      static const char* type_name() { return "IM"; }
    };
    template<> struct traits_asptr<casadi::Matrix<casadi_int> > {
      static bool conv(GUESTOBJECT* p, casadi::IM* m) {
        using namespace casadi;
        // Object is a sparsity pattern
        {
          void *sp = 0;
          if (SWIG_IsOK(SWIG_ConvertPtr(p, &sp, $descriptor(casadi::Sparsity*), 0))) {
            if (m) *m = IM::ones(*static_cast<Sparsity*>(sp));
            return true;
          }
        }

        // First convert to integer
        {
          casadi_int tmp;
          if (swig::assign(p, m? &tmp: 0)) {
            if (m) *m = tmp;
            return true;
          }
        }

#ifdef SWIGPYTHON
        // Numpy arrays will be cast to dense Matrix<casadi_int>
        IM* mp = m;
        if (IM_from_array(p, m ? &mp : 0)) return true;

        if (PyErr_Occurred()) PyErr_Clear(); // Clear pending exception before type check
        if (PyObject_HasAttrString(p,"__IM__")) {
          PyObject *cr = PyObject_CallMethod(p, (char*) "__IM__", 0);
          if (!cr) return false;
          bool result = swig::assign(cr, m);
          Py_DECREF(cr);
          return result;
        }

        {
          std::vector <casadi_int> t;
          if (swig::assign(p, m ? &t : 0)) {
            if (m) *m = casadi::Matrix<casadi_int>(t);
            return true;
          }
        }

#elif defined(SWIGMATLAB)
        // In MATLAB, it is common to use floating point values to represent integers
        if (mxIsDouble(p) && mxGetNumberOfDimensions(p)==2) {
          double* data = static_cast<double*>(mxGetData(p));

          // Check if all integers
          bool all_integers=true;
          size_t sz = getNNZ(p);
          for (size_t i=0; i<sz; ++i) {
            if (data[i] != casadi_int(data[i])) {
              all_integers = false;
              break;
            }
          }

          // If successful
          if (all_integers) {
            if (m) {
              *m = casadi::IM(get_sparsity(p));
              for (size_t i=0; i<sz; ++i) {
                (*m)->at(i) = casadi_int(data[i]);
              }
            }
            return true;
          }
        }

#endif
        // Convert from DM
        {
          DM tmp;
          if (swig::assign(p, m? &tmp: 0)) {
            // Check integrality
            for (double d : tmp.nonzeros()) {
              if (d!=casadi_int(d)) return false;
            }
            // Convert
            if (m) *m = casadi::Matrix<double>(tmp);
            return true;
          }
        }

        // No match
        return false;
      }
      static int asptr(GUESTOBJECT* p, casadi::IM** m) {
        return casadi_asptr<casadi::IM>(p, m, 0, conv);
      }
      static bool in_place(GUESTOBJECT* p, casadi::IM** m) {
        return casadi_in_place<casadi::IM>(p, m, 0, conv);
      }
    };
#ifdef SWIGJULIA
    template<> struct traits_asval<casadi::Matrix<casadi_int> > : casadi_asval<casadi::IM> {};
#endif
    // An IM leaves as a DM
    template<> struct traits_from<casadi::Matrix<casadi_int> > {
      static GUESTOBJECT* from(const casadi::IM& a) {
        return swig::from(casadi::DM(a));
      }
    };
  } // namespace swig
}

%fragment(SWIG_Traits_frag(casadi::Matrix<casadi::SXElem>), "header", fragment="casadi_aux",
          fragment=SWIG_Traits_frag(casadi::Matrix<double>)) {
  namespace swig {
    template<> struct traits<casadi::Matrix<casadi::SXElem> > {
      typedef pointer_category category;
      static const char* type_name() { return "SX"; }
    };
    template<> struct traits_asptr<casadi::Matrix<casadi::SXElem> > {
      static bool conv(GUESTOBJECT* p, casadi::SX* m) {
#ifdef SWIGPYTHON
        using namespace casadi;
        // Object has __SX__ method; checked before the numeric heuristics
        if (PyErr_Occurred()) PyErr_Clear(); // Clear pending exception before type check
        if (PyObject_HasAttrString(p,"__SX__")) {
          PyObject *cr = PyObject_CallMethod(p, (char*) "__SX__", 0);
          if (!cr) return false;
          if (cr == Py_None) { Py_DECREF(cr); }
          else {
            // swig::assign COPIES into *m; the hook result dies with cr
            bool flag = swig::assign(cr, m);
            Py_DECREF(cr);
            if (flag) return true;
          }
        }

#endif
        // Try first converting to a temporary DM
        {
          casadi::DM tmp;
          if (swig::assign(p, m? &tmp: 0)) {
            if (m) *m = tmp;
            return true;
          }
        }
#ifdef SWIGPYTHON
        // Clear any error from a failed DM conversion
        if (PyErr_Occurred()) PyErr_Clear();
        // Numpy arrays will be cast to dense SX
        SX* mp = m;
        if (SX_from_array(p, m ? &mp : 0)) return true;

#endif
        // No match
        return false;
      }
      static int asptr(GUESTOBJECT* p, casadi::SX** m) {
        return casadi_asptr<casadi::SX>(p, m, $descriptor(casadi::Matrix<casadi::SXElem>*), conv);
      }
      static bool in_place(GUESTOBJECT* p, casadi::SX** m) {
        return casadi_in_place<casadi::SX>(p, m, $descriptor(casadi::Matrix<casadi::SXElem>*), conv);
      }
    };
#ifdef SWIGJULIA
    template<> struct traits_asval<casadi::Matrix<casadi::SXElem> > : casadi_asval<casadi::SX> {};
#endif
    template<> struct traits_from<casadi::Matrix<casadi::SXElem> > {
      static GUESTOBJECT* from(const casadi::SX& a) {
        return SWIG_NewPointerObj(new casadi::SX(a), $descriptor(casadi::Matrix<casadi::SXElem> *), SWIG_POINTER_OWN);
      }
    };
  } // namespace swig
}

%fragment(SWIG_Traits_frag(casadi::SXElem), "header", fragment="casadi_aux",
          fragment=SWIG_Traits_frag(casadi::Matrix<casadi::SXElem>)) {
  namespace swig {
    template<> struct traits<casadi::SXElem> {
      typedef pointer_category category;
      static const char* type_name() { return "SXElem"; }
    };
    template<> struct traits_asptr<casadi::SXElem> {
      static bool conv(GUESTOBJECT* p, casadi::SXElem* m) {
        // Try first converting to a temporary SX
        casadi::SX tmp, *mt=&tmp;
        if(swig::asref(p, m ? &mt : 0)) {
          if (m && !mt->is_scalar()) return false;
          if (m) *m = mt->scalar();
          return true;
        }
        // No match
        return false;
      }
      static int asptr(GUESTOBJECT* p, casadi::SXElem** m) {
        return casadi_asptr<casadi::SXElem>(p, m, 0, conv);
      }
      static bool in_place(GUESTOBJECT* p, casadi::SXElem** m) {
        return casadi_in_place<casadi::SXElem>(p, m, 0, conv);
      }
    };
#ifdef SWIGJULIA
    template<> struct traits_asval<casadi::SXElem> : casadi_asval<casadi::SXElem> {};
#endif
    // An SXElem leaves as an SX
    template<> struct traits_from<casadi::SXElem> {
      static GUESTOBJECT* from(const casadi::SXElem& a) {
        return swig::from(casadi::SX(a));
      }
    };
  } // namespace swig
}

%fragment(SWIG_Traits_frag(casadi::MX), "header", fragment="casadi_aux",
          fragment=SWIG_Traits_frag(casadi::Matrix<double>)) {
  namespace swig {
    template<> struct traits<casadi::MX> {
      typedef pointer_category category;
      static const char* type_name() { return "MX"; }
    };
    template<> struct traits_asptr<casadi::MX> {
      static bool conv(GUESTOBJECT* p, casadi::MX* m) {
#ifdef SWIGPYTHON
        using namespace casadi;
        // Object has __MX__ method; checked before the numeric heuristics
        if (PyErr_Occurred()) PyErr_Clear(); // Clear pending exception before type check
        if (PyObject_HasAttrString(p,"__MX__")) {
          PyObject *cr = PyObject_CallMethod(p, (char*) "__MX__", 0);
          if (!cr) return false;
          if (cr == Py_None) { Py_DECREF(cr); }
          else {
            // swig::assign copies into *m
            bool flag = swig::assign(cr, m);
            Py_DECREF(cr);
            if (flag) return true;
          }
        }

#endif
        // Try first converting to a temporary DM
        {
          casadi::DM tmp;
          if (swig::assign(p, m ? &tmp : 0)) {
            if (m) *m = tmp;
            return true;
          }
        }
#ifdef SWIGPYTHON
        // Clear any error from a failed DM conversion
        if (PyErr_Occurred()) PyErr_Clear();

#endif
        // No match
        return false;
      }
      static int asptr(GUESTOBJECT* p, casadi::MX** m) {
        return casadi_asptr<casadi::MX>(p, m, $descriptor(casadi::MX*), conv);
      }
      static bool in_place(GUESTOBJECT* p, casadi::MX** m) {
        return casadi_in_place<casadi::MX>(p, m, $descriptor(casadi::MX*), conv);
      }
    };
#ifdef SWIGJULIA
    template<> struct traits_asval<casadi::MX> : casadi_asval<casadi::MX> {};
#endif
    template<> struct traits_from<casadi::MX> {
      static GUESTOBJECT* from(const casadi::MX& a) {
        return SWIG_NewPointerObj(new casadi::MX(a), $descriptor(casadi::MX*), SWIG_POINTER_OWN);
      }
    };
  } // namespace swig
}

// GenericType: declarations here, bodies in casadi_generictype (they need the container traits of GenericType itself)
%fragment(SWIG_Traits_frag(casadi::GenericType), "header", fragment="casadi_aux") {
  namespace swig {
    template<> struct traits<casadi::GenericType> {
      typedef pointer_category category;
      static const char* type_name() { return "GenericType"; }
    };
#ifdef SWIGJULIA
    // Returns as the Julia value, not the proxy
    template<> struct julia_eltype<casadi::GenericType> {
      static jl_value_t* get() { return (jl_value_t*)jl_any_type; }
    };
#endif
    template<> struct traits_asptr<casadi::GenericType> {
      static bool conv(GUESTOBJECT* p, casadi::GenericType* m);
      static int asptr(GUESTOBJECT* p, casadi::GenericType** m) {
#ifdef CASADI_NOT_NONE
        if (CASADI_NOT_NONE(p)) return casadi_asptr<casadi::GenericType>(p, m, $descriptor(casadi::GenericType*), conv);
        if (m) *m = new casadi::GenericType();
        return m ? SWIG_NEWOBJ : SWIG_OK;
#else
        return casadi_asptr<casadi::GenericType>(p, m, $descriptor(casadi::GenericType*), conv);
#endif
      }
      static bool in_place(GUESTOBJECT* p, casadi::GenericType** m) {
#ifdef CASADI_NOT_NONE
        if (CASADI_NOT_NONE(p)) return casadi_in_place<casadi::GenericType>(p, m, $descriptor(casadi::GenericType*), conv);
        **m = casadi::GenericType();
        return true;
#else
        return casadi_in_place<casadi::GenericType>(p, m, $descriptor(casadi::GenericType*), conv);
#endif
      }
    };
#ifdef SWIGJULIA
    template<> struct traits_asval<casadi::GenericType> : casadi_asval<casadi::GenericType> {};
#endif
    template<> struct traits_from<casadi::GenericType> {
      static GUESTOBJECT* from(const casadi::GenericType& a);
    };
  } // namespace swig
}

%fragment("casadi_generictype", "header", fragment=SWIG_Traits_frag(casadi::GenericType),
          fragment="casadi_containers", fragment=SWIG_Traits_frag(bool), fragment=SWIG_Traits_frag(long long),
          fragment=SWIG_Traits_frag(double), fragment=SWIG_Traits_frag(std::string),
          fragment=SWIG_Traits_frag(casadi::Function)) {
  namespace swig {
    // Try one alternative of GenericType
    template<class T> bool casadi_generic(GUESTOBJECT* p, casadi::GenericType* m) {
#ifdef SWIGMATLAB
      if (!m) return swig::check<T>(p);
      T tmp;
      if (!SWIG_IsOK(swig::asval<T>(p, &tmp))) return false;
      *m = casadi::GenericType(tmp);
#else
      if (!m) return CASADI_OK(swig::check<T>(p));
      T tmp, *tmp_ptr = &tmp;
      if (!CASADI_OK(swig::asref(p, &tmp_ptr))) return false;
      *m = casadi::GenericType(*tmp_ptr);
#endif
      return true;
    }

#ifdef SWIGJAVA
    // An empty Java array of objects still has an element type
    inline bool java_empty_generic(jobject p, casadi::GenericType* m) {
      using casadi::GenericType;
      using casadi::java_class;
      JNIEnv* e = SWIG_JavaGuestEnvChecked();
      static jclass s = java_class("[Ljava/lang/String;"), f = java_class("[Lorg/casadi/Function;"),
        dd = java_class("[[D"), ll = java_class("[[J"), ss = java_class("[[Ljava/lang/String;");
      GenericType g;
      if (e->IsInstanceOf(p, s)) g = std::vector<std::string>();
      else if (e->IsInstanceOf(p, f)) g = std::vector<casadi::Function>();
      else if (e->IsInstanceOf(p, dd)) g = std::vector<std::vector<double> >();
      else if (e->IsInstanceOf(p, ll)) g = std::vector<std::vector<casadi_int> >();
      else if (e->IsInstanceOf(p, ss)) g = std::vector<std::vector<std::string> >();
      else return false;
      if (e->GetArrayLength(static_cast<jarray>(p)) != 0) return false;
      if (m) *m = g;
      return true;
    }
#endif
    bool traits_asptr<casadi::GenericType>::conv(GUESTOBJECT* p, casadi::GenericType* m) {
      using casadi::GenericType;
#ifdef SWIGJAVA
      if (java_empty_generic(p, m)) return true;
#endif
      // Try to convert to different types
      return casadi_generic<bool>(p, m)
          || casadi_generic<casadi_int>(p, m)
          || casadi_generic<double>(p, m)
          || casadi_generic<std::string>(p, m)
#ifdef SWIGPYTHON
          || casadi_generic<GenericType::Dict>(p, m) // a dict before any sequence
#endif
          || casadi_generic<std::vector<bool> >(p, m)
          || casadi_generic<std::vector<casadi_int> >(p, m)
          || casadi_generic<std::vector<double> >(p, m)
          || casadi_generic<std::vector<std::string> >(p, m)
          || casadi_generic<std::vector<std::vector<std::string> > >(p, m)
          || casadi_generic<std::vector<std::vector<casadi_int> > >(p, m)
          || casadi_generic<std::vector<std::vector<double> > >(p, m)
          || casadi_generic<casadi::Function>(p, m)
          || casadi_generic<std::vector<casadi::Function> >(p, m)
#ifndef SWIGPYTHON
          || casadi_generic<GenericType::Dict>(p, m)
#endif
          || casadi_generic<std::vector<GenericType::Dict> >(p, m)
          || casadi_generic<std::vector<std::vector<GenericType> > >(p, m)
          || casadi_generic<std::vector<GenericType> >(p, m);
    }

    GUESTOBJECT* traits_from<casadi::GenericType>::from(const casadi::GenericType& a) {
      using namespace casadi;
      switch (a.getType()) {
      case OT_BOOL: return swig::from(a.as_bool());
      case OT_INT: return swig::from(a.as_int());
      case OT_DOUBLE: return swig::from(a.as_double());
      case OT_STRING: return swig::from(a.as_string());
      case OT_INTVECTOR: return swig::from(a.as_int_vector());
      case OT_INTVECTORVECTOR: return swig::from(a.as_int_vector_vector());
      case OT_BOOLVECTOR: return swig::from(a.as_bool_vector());
      case OT_DOUBLEVECTOR: return swig::from(a.as_double_vector());
      case OT_DOUBLEVECTORVECTOR: return swig::from(a.as_double_vector_vector());
      case OT_STRINGVECTOR: return swig::from(a.as_string_vector());
      case OT_STRINGVECTORVECTOR: return swig::from(a.as_string_vector_vector());
      case OT_DICT: return swig::from(a.as_dict());
      case OT_DICTVECTOR: return swig::from(a.as_dict_vector());
      case OT_VECTORVECTOR: return swig::from(a.as_vector_vector());
      case OT_VECTOR: return swig::from(a.as_vector());
      case OT_FUNCTION: return swig::from(a.as_function());
      case OT_FUNCTIONVECTOR: return swig::from(a.as_function_vector());
#ifdef SWIGPYTHON
      case OT_NULL:
      case OT_VOIDPTR:
        // A new reference, like every other branch
        Py_INCREF(Py_None);
        return Py_None;
#elif defined(SWIGMATLAB)
      case OT_VOIDPTR:
        return mxCreateDoubleScalar(0);
#elif defined(SWIGJULIA)
      case OT_NULL:
      case OT_VOIDPTR:
        return jl_nothing;
#endif
      default: return 0;
      }
    }
  } // namespace swig
}

// Collect all fragments
%fragment("casadi_all", "header", fragment="casadi_aux,casadi_extra,casadi_containers,casadi_generictype",
          fragment=SWIG_Traits_frag(casadi::Sparsity), fragment=SWIG_Traits_frag(casadi::Function),
          fragment=SWIG_Traits_frag(casadi::Slice), fragment=SWIG_Traits_frag(casadi::Matrix<double>),
          fragment=SWIG_Traits_frag(casadi::Matrix<casadi_int>), fragment=SWIG_Traits_frag(casadi::Matrix<casadi::SXElem>),
          fragment=SWIG_Traits_frag(casadi::SXElem), fragment=SWIG_Traits_frag(casadi::MX)
#ifdef SWIGWASMJS
          , fragment=SWIG_Traits_frag(casadi::Linsol)
#endif
          ) { }

// Can be overloaded by specifying before importing casadi.i
%fragment("casadi_extra", "header") {}


#endif // SWIGXML

 // Input typemaps; xStubIn/xTsStub are the .pyi/.d.ts annotations
#ifdef SWIGJAVA
 // Director arguments with the JNI descriptor of their Java type
%define %java_type_director(xDesc, xType...)
%typemap(directorin, descriptor=xDesc, noblock=1, fragment="casadi_all") xType, const xType, const xType& { $input = swig::from< xType >($1); }
%enddef
 // Java type of a casadi-typemapped C++ type (parameters, returns, outputs, directors)
%define %java_type(xJava, xDesc, xType...)
%typemap(jstype) xType, const xType, const xType&, xType &OUTPUT, xType &INOUT xJava
%java_type_director(xDesc, xType)
%enddef
 // ... whose parameters also take the Java types in xWiden
%define %java_type_widen(xJava, xDesc, xWiden, xType...)
%java_type(xJava, xDesc, xType)
%typemap(jstype, widen=xWiden) const xType&, xType xJava
%enddef
 // java.util.Map<String, xValue> and java.util.Map.Entry<xKey, xValue>
%define %java_type_map(xValue, xType...)
%java_type("java.util.Map<String," xValue ">", "Ljava/util/Map;", xType)
%enddef
%define %java_type_entry(xKey, xValue, xType...)
%java_type("java.util.Map.Entry<" xKey "," xValue ">", "Ljava/util/Map$Entry;", xType)
%enddef
#endif // SWIGJAVA
%define %casadi_input_typemaps(xName, xStubIn, xTsStub, xPrec, xType...)
#ifdef SWIGWASMJS
 // Wrapped types cross the wasm boundary as EM_VAL
%typemap(ctype) xType, const xType&, xType& "EM_VAL"
#endif
#ifdef SWIGJULIA
%julia_value_type(xType)
#endif
#ifdef SWIGJAVA
 // Values cross as Object and convert in C++; %java_type sets the Java type
%typemap(jni) xType, const xType, const xType&, xType &INOUT "jobject"
%typemap(jtype) xType, const xType, const xType&, xType &INOUT "Object"
%typemap(jstype) xType, const xType, const xType&, xType &OUTPUT, xType &INOUT "Object"
%typemap(javain) xType, const xType&, xType &INOUT "$javainput"
%typemap(javaout) xType, const xType, const xType& { return $imclassname.swigCast($jnicall); }
%typemap(javadirectorin) xType, const xType, const xType& "$imclassname.swigCast($jniinput)"
%typemap(javadirectorout) xType, const xType "$javacall"
#endif
 // Pass input by value, check if matches
%typemap(typecheck, noblock=1, precedence=xPrec, fragment="casadi_all") xType {
  $1 = swig::check< xType >($input);
 }

 // Directorout typemap; as input by value
%typemap(directorout, noblock=1, fragment="casadi_all") xType {
    if (!swig::assign< xType >($input, &$result)) {
      %dirout_fail(SWIG_TypeError,"$type");
    }
 }

 // Pass input by value, convert argument
%typemap(in, doc=xName, pystub_in=xStubIn, tsstub_in=xTsStub, noblock=1, fragment="casadi_all") xType {
  if (!swig::assign< xType >($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input $argnum to type '" xName "'.");
 }

 // Pass input by value, cleanup
%typemap(freearg, noblock=1) xType {}

 // Pass input by reference, check if matches
%typemap(typecheck, noblock=1, precedence=xPrec, fragment="casadi_all") const xType& {
  $1 = swig::check< xType >($input);
 }

 // Pass input by reference, convert argument
%typemap(in, doc=xName, pystub_in=xStubIn, tsstub_in=xTsStub, noblock=1, fragment="casadi_all") const xType & (xType m) {
  $1 = &m;
  if (!swig::asref< xType >($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input $argnum to type '" xName "'.");
 }

 // Pass input by reference, cleanup
%typemap(freearg, noblock=1) const xType & {}

#ifdef SWIGWASMJS
 // Non-const reference: standard SWIG convention
%typemap(in, doc=xName, pystub_in=xStubIn, tsstub_in=xTsStub, noblock=1) xType & {
  if (!SWIG_IsOK(SWIG_ConvertPtr($input, (void**)&$1, $1_descriptor, 0))) {
    SWIG_exception_fail(SWIG_TypeError, "Failed to convert input $argnum to type '" xName "'.");
  }
 }
%typemap(freearg, noblock=1) xType & {}
#endif

%enddef

 // Output typemaps; xStubOut/xTsStub are the .pyi/.d.ts annotations
%define %casadi_output_typemaps(xName, xStubIn, xStubOut, xTsStub, xType...)

 // Return-by-value; Julia returns follow the Lib convention
#ifdef SWIGJAVA
 // null is a Java value
%typemap(out, noblock=1, fragment="casadi_all") xType, const xType { $result = swig::from< xType >($1); }
%typemap(out, noblock=1, fragment="casadi_all") const xType& { $result = swig::from< xType >(*$1); }
#elif !defined(SWIGJULIA)
%typemap(out, doc=xName, pystub_out=xStubOut, tsstub_out=xTsStub, noblock=1, fragment="casadi_all") xType, const xType {
  if(!($result = swig::from< xType >($1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type '" xName "'.");
}

// Return a const-ref behaves like return-by-value
%typemap(out, doc=xName, pystub_out=xStubOut, tsstub_out=xTsStub, noblock=1, fragment="casadi_all") const xType& {
  if(!($result = swig::from< xType >(*$1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type '" xName "'.");
}
#endif

// Inputs marked OUTPUT are also returned by the function, ...
%typemap(argout, noblock=1,fragment="casadi_all") xType &OUTPUT {
  %append_output(swig::from(*$1));
 }

// ... and the corresponding inputs are ignored.
%typemap(in, doc=xName, pystub_out=xStubOut, tsstub_out=xTsStub, noblock=1, numinputs=0) xType &OUTPUT (xType m) {
 $1 = &m;
}

 // Directorin typemap; as output
%typemap(directorin, noblock=1, fragment="casadi_all") xType, const xType {
    if(!($input = swig::from< xType >($1))) %dirout_fail(SWIG_TypeError,"For director inputs, failed to convert input to " xName ".");
 }

%typemap(directorin, noblock=1, fragment="casadi_all") const xType& {
    if(!($input = swig::from< xType >($1))) %dirout_fail(SWIG_TypeError,"For director inputs, failed to convert input to " xName ".");
 }
#ifdef SWIGJAVA
%java_type_director("Ljava/lang/Object;", xType)
#endif

 // Enable dynamic dispatch
%typemap(typecheck, noblock=1, fragment="casadi_all") xType &OUTPUT {
  $1 = swig::check< xType >($input);
 }

// Alternative names
%apply xType &OUTPUT {xType &OUTPUT1};
%apply xType &OUTPUT {xType &OUTPUT2};
%apply xType &OUTPUT {xType &OUTPUT3};
%apply xType &OUTPUT {xType &OUTPUT4};
%apply xType &OUTPUT {xType &OUTPUT5};
%apply xType &OUTPUT {xType &OUTPUT6};

// Inputs marked INOUT are also returned by the function, ...
%typemap(argout,noblock=1,fragment="casadi_all") xType &INOUT {
  %append_output(swig::from(*$1));
 }

// ... but kept as inputs
%typemap(in, doc=xName, pystub_in=xStubIn, pystub_out=xStubOut, tsstub_in=xTsStub, tsstub_out=xTsStub, noblock=1, fragment="casadi_all") xType &INOUT (xType m) {
  $1 = &m;
  if (!swig::asref< xType >($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input to type '" xName "'.");
 }

 // ... also for dynamic dispatch
%typemap(typecheck, noblock=1, fragment="casadi_all") xType& INOUT {
  $1 = swig::check< xType >($input);
 }

// No arguments need to be freed
%typemap(freearg, noblock=1) xType& INOUT {}

// Alternative names
%apply xType &INOUT {xType &INOUT1};
%apply xType &INOUT {xType &INOUT2};
%apply xType &INOUT {xType &INOUT3};
%apply xType &INOUT {xType &INOUT4};
%apply xType &INOUT {xType &INOUT5};
%apply xType &INOUT {xType &INOUT6};
%enddef

 // All-in-one template instantiation + typemaps
%define %casadi_template(xName, xStubIn, xStubOut, xTsStub, xPrec, xType...)
%template() xType;
%casadi_input_typemaps(xName, xStubIn, xTsStub, xPrec, xType)
%casadi_output_typemaps(xName, xStubIn, xStubOut, xTsStub, %arg(xType))
%enddef

// std::vector<E> variant; wasm-js also gets a named %template and %wasm_vec
%define %casadi_template_vec(xName, xStubIn, xStubOut, xTsStub, xPrec, xJSName, xElemType...)
%casadi_template(xName, xStubIn, xStubOut, xTsStub, xPrec, std::vector< xElemType >)
#ifdef SWIGWASMJS
%template(xJSName) std::vector< xElemType >;
%wasm_vec(xElemType, xJSName)
#endif
%enddef

 // Input + output typemaps for an existing type.
%define %casadi_typemaps(xName, xStubIn, xStubOut, xTsStub, xPrec, xType...)
%casadi_input_typemaps(xName, xStubIn, xTsStub, xPrec, xType)
%casadi_output_typemaps(xName, xStubIn, xStubOut, xTsStub, xType)
%enddef

#ifdef SWIGWASMJS
// Director typemaps handing wrapped proxies to JS

%define %casadi_director_class(xType...)
// SWIG_NewPointerObj wraps via clientdata
%typemap(directorin, noblock=1) xType, const xType, const xType& {
  {
    xType *__p = new xType($1);
    $input = SWIG_NewPointerObj(__p, $descriptor(xType *), 0);
  }
}
%typemap(directorout, noblock=1) xType, const xType {
  {
    emscripten::val __r = emscripten::val::take_ownership($input);
    EM_VAL __h = __r.release_ownership();
    void *__p = 0;
    (void)SWIG_WASMJS_ConvertPtr(__h, &__p, 0, 0);
    $result = __p ? *static_cast< xType *>(__p) : xType();
  }
}
%enddef

%define %casadi_director_vec(xElemType...)
%typemap(directorin, noblock=1) std::vector< xElemType >,
                                 const std::vector< xElemType >&,
                                       std::vector< xElemType >& {
  // SWIG_NewPointerObj wraps via clientdata
  {
    emscripten::val __a = emscripten::val::array();
    for (size_t __i = 0; __i < $1.size(); ++__i) {
      xElemType *__ep = new xElemType($1[__i]);
      EM_VAL __eh = SWIG_NewPointerObj(__ep, $descriptor(xElemType *), 0);
      __a.call<void>("push", emscripten::val::take_ownership(__eh));
    }
    $input = __a.release_ownership();
  }
}
%typemap(directorout, noblock=1) std::vector< xElemType > {
  {
    emscripten::val __r = emscripten::val::take_ownership($input);
    std::vector< xElemType > __v;
    size_t __n = __r["length"].as<size_t>();
    for (size_t __i = 0; __i < __n; ++__i) {
      emscripten::val __it = __r[__i];
      EM_VAL __h = __it.release_ownership();
      void *__p = 0;
      (void)SWIG_WASMJS_ConvertPtr(__h, &__p, 0, 0);
      if (__p) __v.push_back(*static_cast< xElemType *>(__p));
    }
    $result = __v;
  }
}
%enddef

// std::vector<std::string>: elements are plain JS strings
%define %casadi_director_vec_str
%typemap(directorin, noblock=1) std::vector< std::string >,
                                 const std::vector< std::string >&,
                                       std::vector< std::string >& {
  {
    emscripten::val __a = emscripten::val::array();
    for (size_t __i = 0; __i < $1.size(); ++__i) {
      __a.call<void>("push", emscripten::val(std::string($1[__i])));
    }
    $input = __a.release_ownership();
  }
}
%typemap(directorout, noblock=1) std::vector< std::string > {
  {
    emscripten::val __r = emscripten::val::take_ownership($input);
    std::vector< std::string > __v;
    size_t __n = __r["length"].as<size_t>();
    for (size_t __i = 0; __i < __n; ++__i) {
      __v.push_back(__r[__i].as<std::string>());
    }
    $result = __v;
  }
}
%enddef

// jsout for std::map<std::string, X>: unbox into a native JS object
%define %casadi_dict_jsout(xValueType...)
%typemap(jsout) std::map< std::string, xValueType >,
                 const std::map< std::string, xValueType >&,
                       std::map< std::string, xValueType >&
"M.__swig_release_handle($call)"
%enddef
#endif // SWIGWASMJS

// Order in typemap matching: Lower value means will be checked first

%define PREC_DICT 21 %enddef
%define PREC_SPARSITY 90 %enddef
%define PREC_IVector 92 %enddef
%define PREC_IVectorVector 92 %enddef
%define PREC_VECTOR 92 %enddef
%define PREC_SLICE 94 %enddef
%define PREC_IM 97 %enddef
%define PREC_DM 100 %enddef
%define PREC_DMVector 101 %enddef
%define PREC_DMVectorVector 101 %enddef
%define PREC_SX 103 %enddef
%define PREC_SXVector 103 %enddef
%define PREC_SXVectorVector 103 %enddef
%define PREC_MX 104 %enddef
%define PREC_MXVector 105 %enddef
%define PREC_MXVectorVector 106 %enddef
%define PREC_STRING 180 %enddef
%define PREC_FUNCTION 200 %enddef
%define PREC_GENERICTYPE 201 %enddef

#ifndef SWIGXML

// std::ostream & is redirected to casadi::uout()
%typemap(in, noblock=1, numinputs=0) std::ostream &stream {
  $1 = &casadi::uout();
}

// Add trailing newline in MATLAB and Octave
#if defined(SWIGMATLAB) || defined(SWIGOCTAVE)
%typemap(argout, noblock=1) std::ostream &stream {
  *$1 << "\n" << std::flush;
}
#endif

#ifdef SWIGJULIA
  #define L_INT "Int"
  #define L_BOOL "Bool"
  #define LPAIR(A,B) "Tuple{" A "," B "}"
  #define L_DOUBLE "Float64"
  #define L_DICT "Dict"
  #define LDICT(ARG) "Dict{String," ARG "}"
  #define LL "Vector{"
  #define LR "}"
  #define L_STR "String"
#else
#define L_INT "int"
#define L_BOOL "bool"
#define LPAIR(A,B) "(" A "," B ")"
#endif

#ifdef SWIGJULIA
#elif defined(SWIGMATLAB) || defined(SWIGOCTAVE)
  #define L_DOUBLE "double"
  #define L_DICT "struct"
  #define LDICT(ARG) L_DICT ":" ARG
  #define LL "{"
  #define LR "}"
  #define L_STR "char"
  #define MATLABSTYLE
#else
  #define LL "["
  #define LR "]"
  #define L_DICT "dict"
  #define L_DOUBLE "float"
  #define LDICT(ARG) L_DICT ":" ARG
  #define L_STR "str"
#endif

#ifdef SWIGPYTHON
%typemap(in, doc="buffer(ro)", pystub_in="memoryview", noblock=1, fragment="casadi_all") (const double * a, casadi_int size) (Py_buffer _global_pybuf_ro) {
  if (PyObject_GetBuffer($input, &_global_pybuf_ro, PyBUF_SIMPLE) != 0) {
    SWIG_exception_fail(SWIG_TypeError, "Must supply a buffer-supporting object.");
  }
  $1 = static_cast<double*>(_global_pybuf_ro.buf);
  $2 = _global_pybuf_ro.len;
 }
%typemap(freearg) (const double * a, casadi_int size) {
  PyBuffer_Release(&_global_pybuf_ro);
 }

%typemap(in, doc="buffer(rw)", pystub_in="memoryview", noblock=1, fragment="casadi_all") (double * a, casadi_int size)  (Py_buffer _global_pybuf_rw) {
  if (PyObject_GetBuffer($input, &_global_pybuf_rw, PyBUF_WRITABLE) != 0) {
    SWIG_exception_fail(SWIG_TypeError, "Must supply a writable buffer-supporting object.");
  }
  $1 = static_cast<double*>(_global_pybuf_rw.buf);
  $2 = _global_pybuf_rw.len;
 }
%typemap(freearg) (double * a, casadi_int size) {
  PyBuffer_Release(&_global_pybuf_rw);
 }

// Directorin typemap; as output
%typemap(directorin, noblock=1, fragment="casadi_all") (const double** arg, const std::vector<casadi_int>& sizes_arg) (PyObject* my_tuple) {
  PyObject * arg_tuple = PyTuple_New($2.size());
  for (casadi_int i=0;i<$2.size();++i) {

    PyObject* buf = $1[i] ? PyMemoryView_FromMemory(reinterpret_cast<char*>(const_cast<double*>($1[i])), $2[i]*sizeof(double), PyBUF_READ) : SWIG_Py_Void();
    PyTuple_SET_ITEM(arg_tuple, i, buf);
  }
  $input = arg_tuple;
}

%typemap(directorin, noblock=1, fragment="casadi_all") (double** res, const std::vector<casadi_int>& sizes_res) {
  PyObject* res_tuple = PyTuple_New($2.size());
  for (casadi_int i=0;i<$2.size();++i) {
    PyObject* buf = $1[i] ? PyMemoryView_FromMemory(reinterpret_cast<char*>(const_cast<double*>($1[i])), $2[i]*sizeof(double), PyBUF_WRITE) : SWIG_Py_Void();
    PyTuple_SET_ITEM(res_tuple, i, buf);
  }
  $input = res_tuple;
}

%typemap(in, doc="void*", pystub_in="Any", noblock=1, fragment="casadi_all") void* raw {
  $1 = PyCapsule_GetPointer($input, NULL);
}

%typemap(out, doc="void*", pystub_out="Any", noblock=1, fragment="casadi_all") void* {
  $result = PyCapsule_New($1, NULL,NULL);
}
#endif

#ifndef SWIGWASMJS
// wasm-js has its own std::string typemaps
%casadi_typemaps(L_STR, "builtins.str", "builtins.str", "string", PREC_STRING, std::string)
#endif
%casadi_template_vec(LL L_STR LR, "Sequence[builtins.str]", "list[builtins.str]", "string[]", PREC_VECTOR, StringVector, std::string)
%casadi_template(LL LL L_STR LR LR, "Sequence[Sequence[builtins.str]]", "list[list[builtins.str]]", "string[][]", PREC_VECTOR, std::vector<std::vector<std::string> >)

// PEP-484 stub aliases and preamble (-stubs)

#ifdef SWIG_STUBS_ENABLED

%insert("stubs_preamble") %{
# Overloads are ordered narrowest-first on purpose (`sin(float) -> float`
# before `sin(_DM) -> DM`), and DM/SX/MX compare elementwise; pyright
# reports both as overlap/override defects when this file itself is opened.
# test/python/pyright_stubs.py checks the stub with this line stripped.
# pyright: reportOverlappingOverload=false, reportIncompatibleMethodOverride=false
# `str` is also a method name on DM/SX/Opti/...; inside those class bodies a
# bare `str` annotation resolves to the method, so type tokens use builtins.str.
import builtins
import numpy as np
from numpy.typing import NDArray
from collections.abc import Callable
from typing import Protocol, runtime_checkable
_T = TypeVar("_T", DM, SX, MX)
# NZproxy is parameterised with Self, so its TypeVar is bounded, not
# value-constrained: a constrained one only accepts its exact constraints.
_TNZ = TypeVar("_TNZ", bound=DM | SX | MX)

@runtime_checkable
class _SupportsDM(Protocol):
    """Any object exposing a ``__DM__()`` method -- the conversion hook
    casadi uses to turn user-defined numeric types into a DM."""
    def __DM__(self) -> DM: ...

@runtime_checkable
class _SupportsTocsc(Protocol):
    """A scipy.sparse matrix or array, which casadi converts to a DM."""
    def tocsc(self) -> Any: ...

@runtime_checkable
class _SupportsSX(Protocol):
    def __SX__(self) -> SX: ...

@runtime_checkable
class _SupportsMX(Protocol):
    def __MX__(self) -> MX: ...

%}

// issue #2959: stub for casadi.ArrayInterface (a %pythoncode class)
#if defined(SWIG_STUBS_ENABLED) && defined(SWIGPYTHON)
%stubcode %{class ArrayInterface(Generic[_T]):
    """numpy-semantics array view over a casadi DM/SX/MX.  Constructing one
    dispatches to the backing-typed subclass ArrayInterfaceDM/SX/MX."""
    __array_priority__: float
    @overload
    def __new__(cls, value: SX, ndim: int | None = None) -> ArrayInterfaceSX: ...
    @overload
    def __new__(cls, value: MX, ndim: int | None = None) -> ArrayInterfaceMX: ...
    @overload
    def __new__(cls, value: object = None, ndim: int | None = None) -> ArrayInterfaceDM: ...
    def __init__(self, value: object, ndim: int | None = None) -> None: ...
    @property
    def shape(self) -> tuple[int, ...]: ...
    @property
    def ndim(self) -> int: ...
    @property
    def size(self) -> int: ...
    def to_casadi(self) -> _T: ...
    @property
    def T(self) -> ArrayInterface[_T]: ...
    def reshape(self, *shape: int) -> ArrayInterface[_T]: ...
    def transpose(self, *axes: int) -> ArrayInterface[_T]: ...
    def squeeze(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def flatten(self) -> ArrayInterface[_T]: ...
    def ravel(self) -> ArrayInterface[_T]: ...
    @property
    def flat(self) -> Any: ...
    def to_DM(self) -> DM: ...
    def sum(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def mean(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def __getitem__(self, idx: Any) -> ArrayInterface[_T]: ...
    def __setitem__(self, idx: Any, value: Any) -> None: ...
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[ArrayInterface[_T]]: ...
    def __array__(self, dtype: Any = None, copy: bool | None = None) -> NDArray[Any]: ...
    def __array_ufunc__(self, ufunc: Any, method: builtins.str, *inputs: Any, **kwargs: Any) -> Any: ...
    def __array_function__(self, func: Any, types: Any, args: Any, kwargs: Any) -> Any: ...
    def __float__(self) -> float: ...
    def __int__(self) -> int: ...
    def __bool__(self) -> bool: ...
    def __add__(self, o: Any) -> ArrayInterface[_T]: ...
    def __radd__(self, o: Any) -> ArrayInterface[_T]: ...
    def __sub__(self, o: Any) -> ArrayInterface[_T]: ...
    def __rsub__(self, o: Any) -> ArrayInterface[_T]: ...
    def __mul__(self, o: Any) -> ArrayInterface[_T]: ...
    def __rmul__(self, o: Any) -> ArrayInterface[_T]: ...
    def __truediv__(self, o: Any) -> ArrayInterface[_T]: ...
    def __rtruediv__(self, o: Any) -> ArrayInterface[_T]: ...
    def __floordiv__(self, o: Any) -> ArrayInterface[_T]: ...
    def __rfloordiv__(self, o: Any) -> ArrayInterface[_T]: ...
    def __mod__(self, o: Any) -> ArrayInterface[_T]: ...
    def __rmod__(self, o: Any) -> ArrayInterface[_T]: ...
    def __pow__(self, o: Any) -> ArrayInterface[_T]: ...
    def __rpow__(self, o: Any) -> ArrayInterface[_T]: ...
    def __matmul__(self, o: Any) -> ArrayInterface[_T]: ...
    def __rmatmul__(self, o: Any) -> ArrayInterface[_T]: ...
    def __lt__(self, o: Any) -> ArrayInterface[_T]: ...
    def __le__(self, o: Any) -> ArrayInterface[_T]: ...
    def __gt__(self, o: Any) -> ArrayInterface[_T]: ...
    def __ge__(self, o: Any) -> ArrayInterface[_T]: ...
    def __eq__(self, o: object) -> ArrayInterface[_T]: ...  # type: ignore[override]
    def __ne__(self, o: object) -> ArrayInterface[_T]: ...  # type: ignore[override]
    __hash__: None  # type: ignore[assignment]
    def __neg__(self) -> ArrayInterface[_T]: ...
    def __pos__(self) -> ArrayInterface[_T]: ...
    def __abs__(self) -> ArrayInterface[_T]: ...
    def prod(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def max(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def min(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def ptp(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def all(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def any(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def var(self, axis: int | tuple[int, ...] | None = None, ddof: int = 0) -> ArrayInterface[_T]: ...
    def std(self, axis: int | tuple[int, ...] | None = None, ddof: int = 0) -> ArrayInterface[_T]: ...
    @classmethod
    def DM(cls, value: Any, shape: Any = ...) -> ArrayInterfaceDM: ...
    @classmethod
    def SX(cls, value: Any, shape: Any = ...) -> ArrayInterfaceSX: ...
    @classmethod
    def MX(cls, value: Any, shape: Any = ...) -> ArrayInterfaceMX: ...

class ArrayInterfaceDM(ArrayInterface[DM]):
    def __DM__(self) -> DM: ...

class ArrayInterfaceSX(ArrayInterface[SX]):
    def __SX__(self) -> SX: ...

class ArrayInterfaceMX(ArrayInterface[MX]):
    def __MX__(self) -> MX: ...

%}
#endif

// Matrix input aliases; SX and MX do not accept each other
%stub_alias_in(DM, bool | int | float | DM | Sparsity | _SupportsDM | _SupportsTocsc | Sequence[bool | int | float] | Sequence[Sequence[bool | int | float]] | NDArray[Any] | %arg(np.floating[Any]) | %arg(np.integer[Any]) | np.bool_)
%stub_alias_in(SX, SX | _SupportsSX | _DM)
%stub_alias_in(MX, MX | _SupportsMX | _DM)
// IM has no Python class: list what swig::asref<IM> accepts
%stub_alias_in(IM, int | Sequence[int] | Sequence[Sequence[int]] | NDArray[Any] | Sparsity | DM)
/* Out: swig::from(IM) hands back a DM (e.g. MX::mapping()). */
%stub_alias_out_key(IM, DM)
%stub_alias_in(SXElem, SXElem | float)

/* Bracketed doc tokens -- not valid SWIG identifiers. */
%stub_alias_in_key([float], Sequence[bool | int | float] | NDArray[Any])
%stub_alias_in_key([int],   Sequence[bool | int] | NDArray[Any])
%stub_alias_in_key([bool],  Sequence[bool | int] | NDArray[Any])

%stub_alias_in(Slice, Slice | int | slice)

// GenericType: precise input union, Any output
// Forward references stay unquoted inside | unions
%stub_alias_in(GenericType, bool | int | float | str | Function | Sequence[_GenericType] | %arg(Mapping[str, _GenericType]) | None)
%stub_alias_out_key(GenericType, Any)

// Matrix index; MX adds MX indices in its own %extend
%stub_alias_in(MIndex, int | slice | Slice | Sequence[bool | int] | NDArray[Any] | Sparsity | DM)

// CasadiMatrix operator overloads, narrowed per class (no MX + SX)
%define %stub_CasadiMatrix_binop(NAME)
%stub_overload_method_selftyped(NAME, DM, DM, other: _DM)
%stub_overload_method_selftyped(NAME, SX, SX, other: _SX)
%stub_overload_method_selftyped(NAME, MX, MX, other: _MX)
%enddef

%define %stub_CasadiMatrix_unop(NAME)
%stub_method0(NAME, Self)
%enddef

#else

%define %stub_CasadiMatrix_binop(NAME) %enddef
%define %stub_CasadiMatrix_unop(NAME)  %enddef

#endif

%casadi_typemaps("Sparsity", "Sparsity", "Sparsity", "Sparsity", PREC_SPARSITY, casadi::Sparsity)
%casadi_template_vec(LL "Sparsity" LR, "Sequence[Sparsity]", "list[Sparsity]", "Sparsity[]", PREC_SPARSITY, SparsityVector, casadi::Sparsity)
%casadi_template(LL LL "Sparsity"  LR  LR, "Sequence[Sequence[Sparsity]]", "list[list[Sparsity]]", "Sparsity[][]", PREC_SPARSITY, std::vector<std::vector< casadi::Sparsity> >)
%casadi_template(LDICT("Sparsity"), "Mapping[builtins.str, Sparsity]", "dict[builtins.str, Sparsity]", "Record<string, Sparsity>", PREC_SPARSITY, std::map<std::string, casadi::Sparsity >)
%casadi_template(LDICT(LL "Sparsity" LR), "Mapping[builtins.str, Sequence[Sparsity]]", "dict[builtins.str, list[Sparsity]]", "Record<string, Sparsity[]>", PREC_SPARSITY, std::map<std::string, std::vector<casadi::Sparsity > >)
%casadi_template(LPAIR(LDICT("Sparsity"),"[" L_STR "]"), "tuple[Mapping[builtins.str, Sparsity], Sequence[builtins.str]]", "tuple[dict[builtins.str, Sparsity], list[builtins.str]]", "[Record<string, Sparsity>, string[]]", PREC_SPARSITY, std::pair<std::map<std::string, casadi::Sparsity >, std::vector<std::string> >)
#ifndef SWIGWASMJS
// wasm-js has its own bool/casadi_int typemaps
%casadi_typemaps(L_BOOL, "bool", "bool", "boolean", SWIG_TYPECHECK_BOOL, bool)
%casadi_template("[" L_BOOL "]", "Sequence[bool | int] | NDArray[Any]", "list[bool]", "boolean[]", SWIG_TYPECHECK_BOOL, std::vector<bool>)
%casadi_template("[[" L_BOOL "]]", "Sequence[Sequence[bool | int]] | Sequence[NDArray[Any]]", "list[list[bool]]", "boolean[][]", SWIG_TYPECHECK_BOOL, std::vector<std::vector<bool> >)
%casadi_typemaps( L_INT , "int", "int", "bigint", SWIG_TYPECHECK_INTEGER, casadi_int)
#endif

#ifdef MATLABSTYLE
#define LABEL "[int,int]"
#elif defined(SWIGJULIA)
#define LABEL LPAIR(L_INT, L_INT)
#else
#define LABEL LPAIR("int","int")
#endif
%casadi_template(LABEL, "tuple[int, int]", "tuple[int, int]", "[bigint, bigint]", SWIG_TYPECHECK_INTEGER, std::pair<casadi_int,casadi_int>)
#undef LABEL
%casadi_template_vec("[" L_INT "]", "Sequence[bool | int] | NDArray[Any]", "list[int]", "bigint[]", PREC_IVector, IntVector, casadi_int)
%casadi_template(LL "[" L_INT "]" LR, "Sequence[Sequence[bool | int] | NDArray[Any]]", "list[list[int]]", "bigint[][]", PREC_IVectorVector, std::vector<std::vector<casadi_int> >)
#ifndef SWIGWASMJS
// wasm_js: double handled by direct typemap (passes through wasm ABI).
%casadi_typemaps(L_DOUBLE, "float", "float", "number", SWIG_TYPECHECK_DOUBLE, double)
#endif
%casadi_template_vec("[" L_DOUBLE "]", "Sequence[bool | int | float] | NDArray[Any]", "list[float]", "number[]", SWIG_TYPECHECK_DOUBLE, DoubleVector, double)
%casadi_template(LL "[" L_DOUBLE "]" LR, "Sequence[Sequence[bool | int | float] | NDArray[Any]]", "list[list[float]]", "number[][]", SWIG_TYPECHECK_DOUBLE, std::vector<std::vector<double> >)
// Most specific first: DM, SX, MX (d.ts overload order)
%casadi_typemaps("DM", "DM", "DM", "DM", PREC_DM, casadi::Matrix<double>)
%casadi_template_vec(LL "DM" LR, "Sequence[DM]", "list[DM]", "DM[]", PREC_DMVector, DMVector, casadi::Matrix<double>)
%casadi_template(LL LL "DM" LR LR, "Sequence[Sequence[DM]]", "list[list[DM]]", "DM[][]", PREC_DMVectorVector, std::vector<std::vector< casadi::Matrix<double> > >)
%casadi_template(LDICT("DM"), "Mapping[builtins.str, DM]", "dict[builtins.str, DM]", "Record<string, DM>", PREC_DM, std::map<std::string, casadi::Matrix<double> >)
%casadi_typemaps("SXElem", "SXElem", "SXElem", "SXElem", PREC_SX, casadi::SXElem)
%casadi_template(LL "SXElem" LR, "Sequence[SXElem]", "list[SXElem]", "SXElem[]", PREC_SXVector, std::vector<casadi::SXElem>)
%casadi_typemaps("SX", "SX", "SX", "SX", PREC_SX, casadi::Matrix<casadi::SXElem>)
%casadi_template_vec(LL "SX" LR, "Sequence[SX]", "list[SX]", "SX[]", PREC_SXVector, SXVector, casadi::Matrix<casadi::SXElem>)
%casadi_template(LL LL "SX" LR LR, "Sequence[Sequence[SX]]", "list[list[SX]]", "SX[][]", PREC_SXVectorVector, std::vector<std::vector< casadi::Matrix<casadi::SXElem> > >)
%casadi_template(LDICT("SX"), "Mapping[builtins.str, SX]", "dict[builtins.str, SX]", "Record<string, SX>", PREC_SX, std::map<std::string, casadi::Matrix<casadi::SXElem> >)
%casadi_typemaps("MX", "MX", "MX", "MX", PREC_MX, casadi::MX)
%casadi_template_vec(LL "MX" LR, "Sequence[MX]", "list[MX]", "MX[]", PREC_MXVector, MXVector, casadi::MX)
%casadi_template(LL LL "MX" LR LR, "Sequence[Sequence[MX]]", "list[list[MX]]", "MX[][]", PREC_MXVectorVector, std::vector<std::vector<casadi::MX> >)
%casadi_template(LDICT("MX"), "Mapping[builtins.str, MX]", "dict[builtins.str, MX]", "Record<string, MX>", PREC_MX, std::map<std::string, casadi::MX>)
%casadi_template(LPAIR("MX","MX"), "tuple[MX, MX]", "tuple[MX, MX]", "[MX, MX]", PREC_MXVector, std::pair<casadi::MX, casadi::MX>)
%casadi_typemaps("IM", "IM", "IM", "IM", PREC_IM, casadi::Matrix<casadi_int>)
%casadi_typemaps("GenericType", "GenericType", "GenericType", "GenericType", PREC_GENERICTYPE, casadi::GenericType)
%casadi_template(LL "GenericType" LR, "Sequence[GenericType]", "list[GenericType]", "GenericType[]", PREC_GENERICTYPE, std::vector<casadi::GenericType>)
%casadi_template(LL LL "GenericType" LR LR, "Sequence[Sequence[GenericType]]", "list[list[GenericType]]", "GenericType[][]", PREC_GENERICTYPE, std::vector<std::vector<casadi::GenericType> >)
%casadi_typemaps("Slice", "Slice", "Slice", "Slice", PREC_SLICE, casadi::Slice)
%casadi_typemaps("Function", "Function", "Function", "Function", PREC_FUNCTION, casadi::Function)
%casadi_template_vec(LL "Function" LR, "Sequence[Function]", "list[Function]", "Function[]", PREC_FUNCTION, FunctionVector, casadi::Function)
%casadi_template(LPAIR("Function","Function"), "tuple[Function, Function]", "tuple[Function, Function]", "[Function, Function]", PREC_FUNCTION, std::pair<casadi::Function, casadi::Function>)
%casadi_template(L_DICT, "Mapping[builtins.str, GenericType]", "dict[builtins.str, GenericType]", "Record<string, GenericType>", PREC_DICT, std::map<std::string, casadi::GenericType>)
%casadi_template(LDICT(LL L_STR LR), "Mapping[builtins.str, Sequence[builtins.str]]", "dict[builtins.str, list[builtins.str]]", "Record<string, string[]>", PREC_DICT, std::map<std::string, std::vector<std::string> >)

#ifdef SWIGJAVA
%java_type("String", "Ljava/lang/String;", std::string)
%java_type("java.util.List<String>", "Ljava/util/List;", std::vector<std::string>)
%java_type("java.util.List<java.util.List<String>>", "Ljava/util/List;", std::vector<std::vector<std::string> >)
%java_type("boolean", "Z", bool)
%java_type("boolean[]", "[Z", std::vector<bool>)
%java_type("boolean[][]", "[[Z", std::vector<std::vector<bool> >)
%java_type("long", "J", casadi_int)
%java_type("long[]", "[J", std::pair<casadi_int,casadi_int>)
%java_type("long[]", "[J", std::vector<casadi_int>)
%java_type("long[][]", "[[J", std::vector<std::vector<casadi_int> >)
%java_type("double", "D", double)
%java_type("double[]", "[D", std::vector<double>)
%java_type("double[][]", "[[D", std::vector<std::vector<double> >)
%java_type("Sparsity", "L$packagepath/Sparsity;", casadi::Sparsity)
%java_type("java.util.List<Sparsity>", "Ljava/util/List;", std::vector<casadi::Sparsity>)
%java_type("java.util.List<java.util.List<Sparsity>>", "Ljava/util/List;", std::vector<std::vector<casadi::Sparsity> >)
%java_type_map("Sparsity", std::map<std::string, casadi::Sparsity>)
%java_type_map("java.util.List<Sparsity>", std::map<std::string, std::vector<casadi::Sparsity> >)
%java_type_entry("java.util.Map<String,Sparsity>", "java.util.List<String>", std::pair<std::map<std::string, casadi::Sparsity>, std::vector<std::string> >)
%java_type_widen("DM", "L$packagepath/DM;", "double", casadi::Matrix<double>)
%java_type("java.util.List<DM>", "Ljava/util/List;", std::vector<casadi::Matrix<double> >)
%java_type("java.util.List<java.util.List<DM>>", "Ljava/util/List;", std::vector<std::vector<casadi::Matrix<double> > >)
%java_type_map("DM", std::map<std::string, casadi::Matrix<double> >)
%java_type("SX", "L$packagepath/SX;", casadi::SXElem)
%java_type("java.util.List<SX>", "Ljava/util/List;", std::vector<casadi::SXElem>)
%java_type_widen("SX", "L$packagepath/SX;", "DM,double", casadi::Matrix<casadi::SXElem>)
%java_type("java.util.List<SX>", "Ljava/util/List;", std::vector<casadi::Matrix<casadi::SXElem> >)
%java_type("java.util.List<java.util.List<SX>>", "Ljava/util/List;", std::vector<std::vector<casadi::Matrix<casadi::SXElem> > >)
%java_type_map("SX", std::map<std::string, casadi::Matrix<casadi::SXElem> >)
%java_type_widen("MX", "L$packagepath/MX;", "DM,double", casadi::MX)
%java_type("java.util.List<MX>", "Ljava/util/List;", std::vector<casadi::MX>)
%java_type("java.util.List<java.util.List<MX>>", "Ljava/util/List;", std::vector<std::vector<casadi::MX> >)
%java_type_map("MX", std::map<std::string, casadi::MX>)
%java_type_entry("MX", "MX", std::pair<casadi::MX, casadi::MX>)
%java_type("DM", "L$packagepath/DM;", casadi::Matrix<casadi_int>)
%java_type("Object", "Ljava/lang/Object;", casadi::GenericType)
%java_type("java.util.List<Object>", "Ljava/util/List;", std::vector<casadi::GenericType>)
%java_type("java.util.List<java.util.List<Object>>", "Ljava/util/List;", std::vector<std::vector<casadi::GenericType> >)
%java_type_widen("Slice", "L$packagepath/Slice;", "long", casadi::Slice)
%java_type("Function", "L$packagepath/Function;", casadi::Function)
%java_type("java.util.List<Function>", "Ljava/util/List;", std::vector<casadi::Function>)
%java_type_entry("Function", "Function", std::pair<casadi::Function, casadi::Function>)
%java_type_map("Object", std::map<std::string, casadi::GenericType>)
// Options are read from any Map with String keys
%typemap(jstype) const std::map<std::string, casadi::GenericType>& "java.util.Map<String,?>"
%java_type_map("java.util.List<String>", std::map<std::string, std::vector<std::string> >)
#endif // SWIGJAVA

#ifdef SWIGJULIA
// Dicts, nested vectors and GenericType return as the Julia value
%define %julia_value_return(xType...)
%typemap(out, ctype="jl_value_t*", jltype="Any", noblock=1, fragment="casadi_all") xType %{ if(!($result = swig::from< xType >($1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output."); %}
%typemap(out, ctype="jl_value_t*", jltype="Any", noblock=1, fragment="casadi_all") const xType& %{ if(!($result = swig::from< xType >(*$1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output."); %}
%typemap(jlout) xType, const xType& "_check($call)"
%enddef
%julia_value_return(casadi::GenericType)
%julia_value_return(std::vector<casadi::GenericType>)
%julia_value_return(std::map<std::string, casadi::Matrix<double> >)
%julia_value_return(std::map<std::string, casadi::Matrix<casadi::SXElem> >)
%julia_value_return(std::map<std::string, casadi::MX>)
%julia_value_return(std::map<std::string, casadi::Sparsity>)
%julia_value_return(std::vector<std::vector<casadi::Matrix<double> > >)
%julia_value_return(std::vector<std::vector<casadi::Matrix<casadi::SXElem> > >)
%julia_value_return(std::vector<std::vector<casadi::MX> >)
%julia_value_return(std::vector<std::vector<casadi::Sparsity> >)
%julia_value_return(std::vector<std::vector<casadi_int> >)
%julia_value_return(std::vector<std::vector<double> >)
%julia_value_return(std::map<std::string, casadi::GenericType>)
%typemap(jlout) casadi::Dict, const casadi::Dict& "_check($call)::AbstractDict"
// std::size_t (hashes) like size_t
%typemap(out, ctype="long long", jltype="Clonglong") std::size_t %{ $result = (long long)$1; %}
%typemap(jlout) std::size_t "Int(_check($call))"
// Matrix<casadi_int> as a Julia Matrix{Int}
%typemap(out, ctype="void *", jltype="Ptr{Cvoid}") casadi::Matrix<casadi_int> %{ $result = (void*) new casadi::Matrix<double>($1); %}
%typemap(jlout) casadi::Matrix<casadi_int> "Base.round.(Int, Base.Matrix{Float64}(DM(_check($call))))"
// an option value is any Julia value GenericType converts from
%typemap(jlparam) const casadi::GenericType& option_value "Any"
// director returns convert as Base.convert does (true for 1, 2 for 2.0), not as the strict input typemaps
%typemap(directorout) bool, casadi_int, double %{ $result = Swig::swig_director_as< $1_ltype >($input); %}
// Julia parameter unions for the matrix types, nested like Python's precedence
%typemap(jlparam) casadi::Matrix<double>, const casadi::Matrix<double>& "Union{DM,Real,AbstractArray{<:Real},Vector{Any}}"
%typemap(jlparam) casadi::Matrix<casadi::SXElem>, const casadi::Matrix<casadi::SXElem>& "Union{SX,DM,Real,AbstractArray{<:Real},Vector{Any}}"
%typemap(jlparam) casadi::MX, const casadi::MX& "Union{MX,DM,Real,AbstractArray{<:Real},Vector{Any}}"
// a Callback is a Function, as in C++
%typemap(jlparam) casadi::Function, const casadi::Function& "Union{Function,CallbackFunction}"
// Containers of matrices: the common element types; other values (Vector{Any}) reach the runtime dispatcher
%define %julia_partial_containers(T, V, D)
%typemap(jlparam, partial="1") std::vector< T >, const std::vector< T >&
  "Union{" V "Vector{DM},Vector{Float64},Vector{Int64}}"
%typemap(jlparam, partial="1") std::map<std::string, T >, const std::map<std::string, T >&
  "Union{" D "Dict{String,DM},Dict{String,Float64},Dict{String,Int64}}"
%enddef
%julia_partial_containers(casadi::Matrix<double>, "", "")
%julia_partial_containers(casadi::Matrix<casadi::SXElem>, "Vector{SX},", "Dict{String,SX},")
%julia_partial_containers(casadi::MX, "Vector{MX},", "Dict{String,MX},")
%typemap(jlparam, partial="1") std::vector<casadi::Sparsity>, const std::vector<casadi::Sparsity>& "Vector{Sparsity}"
#endif // SWIGJULIA

#ifdef SWIGWASMJS
// Director typemaps; after %casadi_typemaps (last registered wins)

%casadi_director_class(casadi::Sparsity)
%casadi_director_class(casadi::Matrix<double>)
%casadi_director_class(casadi::MX)
%casadi_director_class(casadi::Matrix<casadi::SXElem>)
%casadi_director_class(casadi::Matrix<casadi_int>)
%casadi_director_class(casadi::Function)
%casadi_director_class(casadi::Slice)
%casadi_director_class(casadi::SXElem)
// No GenericType: JS sees plain values

%casadi_director_vec(casadi::Sparsity)
%casadi_director_vec(casadi::Matrix<double>)
%casadi_director_vec(casadi::MX)
%casadi_director_vec(casadi::Matrix<casadi::SXElem>)
%casadi_director_vec(casadi::Function)

%casadi_director_vec_str

// jsout for every dict-shaped return type
%casadi_dict_jsout(casadi::GenericType)
%casadi_dict_jsout(casadi::Sparsity)
%casadi_dict_jsout(casadi::Matrix<double>)
%casadi_dict_jsout(casadi::MX)
%casadi_dict_jsout(casadi::Matrix<casadi::SXElem>)
%casadi_dict_jsout(std::vector<std::string>)
%casadi_dict_jsout(std::vector< casadi::Sparsity >)
%casadi_dict_jsout(std::vector< casadi::Matrix<double> >)
%casadi_dict_jsout(std::vector< casadi::MX >)
%casadi_dict_jsout(std::vector< casadi::Matrix<casadi::SXElem> >)
#endif // SWIGWASMJS

#undef L_INT
#undef L_BOOL
#undef LPAIR
#undef L_DOUBLE
#undef L_DICT
#undef LL
#undef LR
#undef L_STR
#undef MATLABSTYLE

// TypeScript aliases for casadi typedefs used in signatures
#ifdef SWIGWASMJS
%insert("stubs") %{
// --- casadi-specific TypeScript type aliases ---
export type Dict         = Record<string, GenericType>;
export type DMDict       = Record<string, DM>;
export type MXDict       = Record<string, MX>;
export type SXDict       = Record<string, SX>;
export type SparsityDict = Record<string, Sparsity>;

// --- Input-position coercion unions (analog of Python pyi's _DM / _SX / _MX).
// xTsStub args in casadi.i use the bare class names (DM / SX / MX / ...);
// the wasm_js.cxx d.ts emitter rewrites bare identifiers to their `_`
// aliases at input positions only, so users can pass JS-natural values
// (numbers, arrays, Sparsity patterns) where DM/SX/MX is expected.  Output
// positions stay as the bare class (DM/SX/MX) so return-type narrowing
// works.
export type _DM = boolean | bigint | number | DM | Sparsity | number[] | number[][];
export type _SX = SX | _DM;
export type _MX = MX | _DM;
export type _IM = IM | bigint | number | number[] | number[][];
export type _SXElem = SXElem | number;
export type _Slice = Slice | number | bigint;
// TS 4.1+ permits inline array/object types in recursive aliases, but a
// bare reference like `_GenericType[]` triggers TS2456 (circular).  Use
// the inline-array form to break the cycle the compiler can see through.
export type _GenericType = boolean | bigint | number | string | Function
                         | readonly _GenericType[]
                         | { readonly [k: string]: _GenericType };
// JS callers regularly pass plain `number` literals where casadi_int
// (-> bigint on the wasm boundary) is expected.  The wasm marshaling
// auto-converts, so widen `bigint` -> `_bigint = bigint | number` on
// INPUT positions only.  Output stays `bigint` for narrow typing of
// size / index / count returns.  Hooked up via %ts_alias_in(bigint)
// below.
export type _bigint = bigint | number;
export type _MIndex = number | bigint | bigint[] | boolean[] | Sparsity | DM;

// Note: DMVector / MXVector / SXVector are NOT type aliases here -- they
// collide with the SWIG-emitted std::vector wrapper classes at the top
// of the d.ts (export class DMVector { ... }).  Functions accepting
// std::vector<DM> use `_DM[]` (input) / `DM[]` (output) via the emitter's
// coercion-union rewrite.  Users who want the wrapper class still get
// `new DMVector(arr)` via the class declaration.
export type DMVectorDict = Record<string, DM[]>;
export type MXVectorDict = Record<string, MX[]>;
export type SXVectorDict = Record<string, SX[]>;
export type IM           = any;  // casadi::Matrix<casadi_int>: no JS-side class today
export type casadi_index = [bigint, bigint];
export type native_DM    = Float64Array;
// Serialization + internal-impl classes -- not exposed to JS.
export type SerializingStream    = any;
export type DeserializingStream  = any;
export type SharedObjectInternal = any;
// Enums that aren't auto-wrapped.
export type ConstraintType = number;
export type VariableType   = number;
export type DomainType     = number;

// ---------------------------------------------------------------------
// JS-side concat aliases.  At runtime, casadi.i's %insert("js")
// block adds variadic + list-form aliases to the concat family:
//   vertcat(a, b, c)  -- variadic
//   vcat([a, b, c])   -- list form
//   (analogous hcat / dcat).
// The auto-generated d.ts only declares the C++-direct names
// (vertcat / horzcat / diagcat).  Re-export the list-form aliases
// so JS-test code calling `M.vcat([...])` typechecks.
export declare const vcat: typeof vertcat;
export declare const hcat: typeof horzcat;
export declare const dcat: typeof diagcat;

// Make Function instances callable at the type level.  At runtime,
// casadi.js wraps Function instances so `f(arg)` evaluates as
// shorthand for `f.call(arg)`.  The auto-generated class body only
// has `.call(...)` so without this merge, `f([3, 7])` would
// typecheck as "instance has no call signature".  This is a
// casadi-specific runtime convention that the wasm_js.cxx emitter
// has no business knowing about.
export interface Function__class {
  (arg: _DM[]): DM[];
  (arg: _SX[]): SX[];
  (arg: _MX[]): MX[];
  (arg: Record<string, _DM>): Record<string, DM>;
  (arg: Record<string, _SX>): Record<string, SX>;
  (arg: Record<string, _MX>): Record<string, MX>;
  // Positional varargs (the runtime accepts `f(x, y, z)` as
  // shorthand for `f([x, y, z])`).  Loose return type since the
  // input flavor is heterogeneous in this form.
  (...args: any[]): any;
}

// Runtime factory shape.
//
// At runtime, `require("./casadi")` returns an async factory function
// (the SWIG `module.exports = async function createcasadi() {...}`
// pattern, see build-wasm/swig/wasm-js/casadi.js).  The factory loads
// the underlying wasm module and returns an object whose properties
// are ALL the exports declared above (`DM`, `SX`, `vertcat`, ...).
//
// Static .ts consumers use the named-export form
// (`import { DM } from "./casadi"`) and never call the factory --
// that path is fine because the d.ts shape AND the JS runtime shape
// happen to agree on the *names* of the exports.
//
// Runtime .js consumers (and TS callers that want to invoke the
// factory) use:
//   const create = require("./casadi");
//   const ca = await create();
//   ca.DM.sym("x");
// To make `ca` typed against the same export shape, the factory's
// return-type is declared self-referentially via
// `typeof import("./casadi")`.  No double bookkeeping needed.
declare const createcasadi: () => Promise<typeof import("./casadi")>;
export default createcasadi;
%}

// d.ts: inputs widen to _<name>, outputs stay bare (GenericType unwraps)
%ts_alias_in(DM)
%ts_alias_in(SX)
%ts_alias_in(MX)
%ts_alias_in(IM)
%ts_alias_in(SXElem)
%ts_alias_in(Slice)
%ts_alias_in(GenericType)
%ts_alias_in(MIndex)
%ts_alias_in(bigint)

%ts_alias_out(GenericType)

// JS-side ergonomics, emitted after the class definitions
%insert("js") %{
  /* Coerce bare numbers/arrays into matrix proxies (jsarg typemaps),
     so `plus(x, 2)` works. */
  const __coerce = (v, Cls) =>
    (v === null || v === undefined ||
     (typeof v === "object" && v._swig_type !== undefined)) ? v : new Cls(v);

  /* Element indexing (wired by %feature("wasmjs:index")), mirroring
     Python: x[k], x[k]=v, x["a:b"], x["0::3"], x["a:b,j"], x.nz[k]
     all route to casadi's 0-based get/set/get_nz/set_nz. */
  const __isIdxKey = (k) => typeof k === "string" && k.length > 0 &&
    (k.charCodeAt(0) === 45 || k.charCodeAt(0) === 58 ||  // '-' or ':'
     (k.charCodeAt(0) >= 48 && k.charCodeAt(0) <= 57));   // digit
  const __idxPart = (s) => {
    if (s.indexOf(":") < 0) return Number(s);
    if (s === ":") return new Slice();
    const p = s.split(":");
    /* Python-style open bounds: "0::3" / ":b" / "1:".  Empty start -> 0;
       empty stop -> the Slice unboundedness sentinel (int64 max, see
       slice.hpp) -- exactly what Python maps a `None` stop to. */
    const start = p[0] === "" ? 0 : Number(p[0]);
    const stop = p[1] === "" ? 9223372036854775807n : Number(p[1]);
    return p.length >= 3 ? new Slice(start, stop, Number(p[2]))
                         : new Slice(start, stop);
  };
  const __idxArgs = (k) => k.split(",").map(__idxPart);
  const __NZ_HANDLER = {
    get(t, k, r) { return __isIdxKey(k) ? t.get_nz(false, Number(k)) : Reflect.get(t, k, r); },
    set(t, k, v, r) { if (__isIdxKey(k)) { t.set_nz(v, false, Number(k)); return true; } return Reflect.set(t, k, v, r); },
  };
  const __INDEX_HANDLER = {
    get(t, k, r) {
      if (k === "__indexed") return true;
      if (k === "nz") return new Proxy(t, __NZ_HANDLER);
      if (__isIdxKey(k)) return t.get(false, ...__idxArgs(k));
      return Reflect.get(t, k, r);
    },
    set(t, k, v, r) {
      if (__isIdxKey(k)) { t.set(v, false, ...__idxArgs(k)); return true; }
      return Reflect.set(t, k, v, r);
    },
  };
  const __mkIndex = (o) =>
    (o === null || typeof o !== "object" || o.__indexed) ? o : new Proxy(o, __INDEX_HANDLER);

  /* x.T as a getter, matching Python. */
  for (const Cls of [DM, MX, SX, Sparsity]) {
    const __orig_T = Cls.prototype.T;
    if (!__orig_T) continue;
    Object.defineProperty(Cls.prototype, "T", {
      get() { return __orig_T.call(this); },
      configurable: true,
    });
  }

  /* Concat family, mirroring Python: vertcat(a, b, c) is variadic,
     vcat([a, b, c]) takes a list (same for horzcat/hcat, diagcat/dcat). */
  const __cat_variadic = (orig) => function (...args) {
    return orig.call(this, args);
  };
  const __cat_list = (orig) => function (arr) {
    return orig.call(this, arr);
  };
  for (const [variadic, list] of [
    ["horzcat", "hcat"],
    ["vertcat", "vcat"],
    ["diagcat", "dcat"],
  ]) {
    const fn = __m[variadic];
    __m[variadic] = __cat_variadic(fn).bind(__m);
    __m[list]     = __cat_list(fn).bind(__m);
  }

  /* Mixed input lists ([DM(0), mx, ...]) must dispatch as the widest
     element type; the generated per-element probes are exact-class, so
     pre-coerce JS-side (matches Python's implicit DM->MX/SX promotion). */
  const __origFnCall = Function.prototype.call;
  Function.prototype.call = function (arg, ...rest) {
    if (Array.isArray(arg)) {
      if (arg.some((e) => e instanceof MX)) arg = arg.map((e) => e instanceof MX ? e : new MX(e));
      else if (arg.some((e) => e instanceof SX)) arg = arg.map((e) => e instanceof SX ? e : new SX(e));
    }
    return __origFnCall.call(this, arg, ...rest);
  };

  // Preserve symbolic lists before permissive scalar-matrix overloads accept arrays.
  const __origExportGraph = __m.export_graph;
  __m.export_graph = function (expression, ...args) {
    if (Array.isArray(expression)) {
      const Cls = expression.every(e => e instanceof SX) ? SXVector
        : expression.every(e => e instanceof MX) ? MXVector : null;
      if (Cls) {
        const vector = __arr_to_vec(expression, Cls);
        try { return __origExportGraph.call(this, vector, ...args); }
        finally { vector.delete(); }
      }
    }
    return __origExportGraph.call(this, expression, ...args);
  };

  /* Function instances are directly callable, mirroring Python:
     f(a, b) positional with n_out-shaped return (0 -> null, 1 -> bare,
     >1 -> list); f({name: v}) dict-in/dict-out.  The dict check must
     exclude casadi types, which are objects too.  Generic machinery
     (prototype transfer, FinalizationRegistry handoff, re-wrap of
     methods returning bare instances) is __callable_class, backend-emitted. */
  Function = __m.Function = __callable_class(Function, __fr_casadi_Function, (f) => function (...args) {
    if (args.length === 1) {
      const a = args[0];
      if (a && typeof a === "object" && !Array.isArray(a)
          && !(a instanceof DM) && !(a instanceof MX) && !(a instanceof SX)
          && !(a instanceof DMVector) && !(a instanceof MXVector) && !(a instanceof SXVector)
          && !(a instanceof Sparsity) && (typeof IM === "undefined" || !(a instanceof IM)))
        return f.call(a);
    }
    const r = f.call(args);
    const n = Number(f.n_out());
    return n === 0 ? null : n === 1 ? r[0] : r;
  });

%}
%insert("js") "wasm-js/load_plugins.js"

#endif

// Matlab is index-1 based
#ifdef SWIGMATLAB
%typemap(in, doc="index", pystub_in="int", noblock=1) casadi_index {
  if (!swig::assign($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input $argnum to type ' index '.");
  if ($1==0) SWIG_exception_fail(SWIG_TypeError,"Index starts at 1, got index '0'.");
  if ($1>=1) $1--;
}
#endif

#endif // SWIGXML

#ifdef SWIGPYTHON
%pythoncode %{
if __name__ != "casadi.casadi":
  raise Exception("""
            CasADi is not running from its package context.

            You probably specified the wrong casadi directory.

            When setting PYTHONPATH or sys.path.append,
            take care not to add a trailing '/casadi'.

        """)

def swigtypeconvertor(*args):
  return swig_typename_convertor_python2cpp(args)

def swig_typename_convertor_python2cpp(a):
  try:
    import numpy as np
  except:
    class NoExist:
      pass
    class Temp(object):
      ndarray = NoExist
    np = Temp()
  if isinstance(a,list):
    if len(a)>0:
      return "[%s]" % "|".join(set([swig_typename_convertor_python2cpp(i) for i in a]))
    else:
      return "[]"
  elif isinstance(a,tuple):
    return "(%s)" % ",".join([swig_typename_convertor_python2cpp(i) for i in a])
  elif isinstance(a,np.ndarray):
    return "np.array(%s)" % ",".join(set([swig_typename_convertor_python2cpp(i) for i in np.array(a).flatten().tolist()]))
  elif isinstance(a,dict):
    if len(a)>0:
      return "|".join(set([swig_typename_convertor_python2cpp(i) for i in a.keys()])) +":"+ "|".join(set([swig_typename_convertor_python2cpp(i) for i in a.values()]))
    else:
      return "dict"
  return type(a).__name__
%}
#endif // SWIGPYTHON

// Init hooks
#ifdef SWIGPYTHON
#ifdef WITH_PYTHON_INTERRUPTS
%{
#include <pythonrun.h>
void SigIntHandler(casadi_int) {
  std::cerr << "Keyboard Interrupt" << std::endl;
  signal(SIGINT, SIG_DFL);
  kill(getpid(), SIGINT);
}
%}

%init %{
PyOS_setsig(SIGINT, SigIntHandler);
%}
#endif // WITH_PYTHON_INTERRUPTS

%pythoncode "numpy_bridge.py"

%pythoncode "casadi_nparray.py"

// No inf/pi stubs: they would clash with numpy's Final declarations

// Stubs for the numpy-flavoured aliases above
%define %stub_unary_math(NAME)
%stub_overload_func(NAME, DM, x: _DM)
%stub_overload_func(NAME, SX, x: _SX)
%stub_overload_func(NAME, MX, x: _MX)
%enddef
%formacro(%stub_unary_math, arcsin, arccos, arctan, arctanh, arcsinh, arccosh)
/* Re-exported from numpy at runtime (`from numpy import pi, inf`). */
#ifdef SWIG_STUBS_ENABLED
%stubcode %{pi: float
inf: float
%}
%stub_overload_func(arctan2, DM, x: _DM, y: _DM)
%stub_overload_func(arctan2, SX, x: _SX, y: _SX)
%stub_overload_func(arctan2, MX, x: _MX, y: _MX)
#endif // SWIG_STUBS_ENABLED
#endif // SWIGPYTHON

#ifdef SWIGWASMJS
// Enum values as 64-bit ints (bigint in JS, like casadi_int)
%feature("wasmjs:enum_int_type","long long");

// Accept a bare JS number/array where SX/MX/DM is expected
%typemap(jstypecheck) casadi::Matrix<casadi::SXElem>, const casadi::Matrix<casadi::SXElem>&, casadi::Matrix<casadi::SXElem>&
  "(M._swig_can_SX(__unwrap($input)) || typeof $input === 'number' || typeof $input === 'bigint' || Array.isArray($input))";
%typemap(jsarg) casadi::Matrix<casadi::SXElem>, const casadi::Matrix<casadi::SXElem>&, casadi::Matrix<casadi::SXElem>&
  "__unwrap(__coerce($input, SX))";
%typemap(jstypecheck) casadi::MX, const casadi::MX&, casadi::MX&
  "(M._swig_can_MX(__unwrap($input)) || typeof $input === 'number' || typeof $input === 'bigint' || Array.isArray($input))";
%typemap(jsarg) casadi::MX, const casadi::MX&, casadi::MX&
  "__unwrap(__coerce($input, MX))";
%typemap(jstypecheck) casadi::Matrix<double>, const casadi::Matrix<double>&, casadi::Matrix<double>&
  "(M._swig_can_DM(__unwrap($input)) || typeof $input === 'number' || typeof $input === 'bigint' || Array.isArray($input))";
%typemap(jsarg) casadi::Matrix<double>, const casadi::Matrix<double>&, casadi::Matrix<double>&
  "__unwrap(__coerce($input, DM))";

// Element indexing via the __mkIndex Proxy
%feature("wasmjs:index") casadi::Matrix<casadi::SXElem>;
%feature("wasmjs:index") casadi::MX;
%feature("wasmjs:index") casadi::Matrix<double>;
#endif // SWIGWASMJS

// Strip leading casadi_ unless followed by ML/int
%rename("%(regex:/casadi_(?!ML|int\\b)(.*)/\\1/)s") "";

%rename(row) get_row;
%rename(colind) get_colind;
%rename(sparsity) get_sparsity;
%rename(nonzeros) get_nonzeros;
%rename(elements) get_elements;

// Explicit conversion to double and casadi_int
#ifdef SWIGPYTHON
%rename(__float__) operator double;
%rename(__int__) operator casadi_int;
#else
%rename(to_double) operator double;
%rename(to_int) operator casadi_int;
#endif
%rename(to_DM) operator Matrix<double>;

#ifdef SWIGJULIA
// Julia concatenates with [a b] / [a; b]; builders stay as _horzcat/_vertcat
%rename(_vertcat) casadi_vertcat;
%rename(_horzcat) casadi_horzcat;
// Python's names; Julia's own spelling is & | ! (and .& .| .!)
%rename(logic_and) casadi_and;
%rename(logic_or) casadi_or;
%rename(logic_not) casadi_not;
%rename(logic_all) casadi_all;
%rename(logic_any) casadi_any;
#endif // SWIGJULIA

#ifdef SWIGPYTHON
%ignore T;

%rename(logic_and) casadi_and;
%rename(logic_or) casadi_or;
%rename(logic_not) casadi_not;
%rename(logic_all) casadi_all;
%rename(logic_any) casadi_any;
%rename(fabs) casadi_abs;

// Concatenations
%rename(_veccat) casadi_veccat;
%rename(_vertcat) casadi_vertcat;
%rename(_horzcat) casadi_horzcat;
%rename(_diagcat) casadi_diagcat;
%pythoncode %{
def veccat(*args):
    try:
        if len(args)==0:
            return DM(0,1)
    except:
        pass
    return _veccat(args)
def vertcat(*args):
    try:
        if len(args)==0:
            return DM(0,1)
    except:
        pass
    return _vertcat(args)
def horzcat(*args):
    try:
        if len(args)==0:
            return DM(1,0)
    except:
        pass
    return _horzcat(args)
def diagcat(*args):
    try:
        if len(args)==0:
            return DM(0,0)
    except:
        pass
    return _diagcat(args)
def vvcat(args):
    try:
        if len(args)==0:
            return DM(0,1)
    except:
        pass
    return _veccat(args)
def vcat(args):
    try:
        if len(args)==0:
            return DM(0,1)
    except:
        pass
    return _vertcat(args)
def hcat(args):
    try:
        if len(args)==0:
            return DM(1,0)
    except:
        pass
    return _horzcat(args)
def dcat(args):
    try:
        if len(args)==0:
            return DM(0,0)
    except:
        pass
    return _diagcat(args)
%}

// Stubs for the *cat wrappers above, narrowest first
%stub_overload_func(veccat,  Sparsity, *args: Sparsity)
%stub_overload_func(veccat,  DM,       *args: _DM)
%stub_overload_func(veccat,  SX,       *args: _SX)
%stub_overload_func(veccat,  MX,       *args: _MX)
%stub_overload_func(vertcat, Sparsity, *args: Sparsity)
%stub_overload_func(vertcat, DM,       *args: _DM)
%stub_overload_func(vertcat, SX,       *args: _SX)
%stub_overload_func(vertcat, MX,       *args: _MX)
%stub_overload_func(horzcat, Sparsity, *args: Sparsity)
%stub_overload_func(horzcat, DM,       *args: _DM)
%stub_overload_func(horzcat, SX,       *args: _SX)
%stub_overload_func(horzcat, MX,       *args: _MX)
%stub_overload_func(diagcat, Sparsity, *args: Sparsity)
%stub_overload_func(diagcat, DM,       *args: _DM)
%stub_overload_func(diagcat, SX,       *args: _SX)
%stub_overload_func(diagcat, MX,       *args: _MX)
%stub_overload_func(vvcat, Sparsity, args: Sequence[Sparsity])
%stub_overload_func(vvcat, DM,       args: Sequence[_DM])
%stub_overload_func(vvcat, SX,       args: Sequence[_SX])
%stub_overload_func(vvcat, MX,       args: Sequence[_MX])
%stub_overload_func(vcat,  Sparsity, args: Sequence[Sparsity])
%stub_overload_func(vcat,  DM,       args: Sequence[_DM])
%stub_overload_func(vcat,  SX,       args: Sequence[_SX])
%stub_overload_func(vcat,  MX,       args: Sequence[_MX])
%stub_overload_func(hcat,  Sparsity, args: Sequence[Sparsity])
%stub_overload_func(hcat,  DM,       args: Sequence[_DM])
%stub_overload_func(hcat,  SX,       args: Sequence[_SX])
%stub_overload_func(hcat,  MX,       args: Sequence[_MX])
%stub_overload_func(dcat,  Sparsity, args: Sequence[Sparsity])
%stub_overload_func(dcat,  DM,       args: Sequence[_DM])
%stub_overload_func(dcat,  SX,       args: Sequence[_SX])
%stub_overload_func(dcat,  MX,       args: Sequence[_MX])

// Non-fatal errors (returning NotImplemented singleton)
%feature("python:maybecall") casadi_plus;
%feature("python:maybecall") casadi_minus;
%feature("python:maybecall") casadi_times;
%feature("python:maybecall") casadi_rdivide;
%feature("python:maybecall") casadi_lt;
%feature("python:maybecall") casadi_le;
%feature("python:maybecall") casadi_eq;
%feature("python:maybecall") casadi_ne;
%feature("python:maybecall") casadi_power;
%feature("python:maybecall") casadi_atan2;
%feature("python:maybecall") casadi_and;
%feature("python:maybecall") casadi_or;
%feature("python:maybecall") casadi_mod;
%feature("python:maybecall") casadi_copysign;
%feature("python:maybecall") casadi_constpow;
#endif // SWIGPYTHON

#ifdef SWIGJULIA
// inv is the matrix inverse in Julia; elementwise inverse is 1 ./ x
%ignore casadi::MX::inv() const;
#endif // SWIGJULIA

#ifdef SWIGJAVA
// Kotlin operator names; Sparsity and Slice keep is_equal / equality members unwrapped
%rename(unaryMinus) operator-() const;
%rename(unaryPlus) operator+() const;
%rename(plus) casadi::Sparsity::operator+;
%rename(times) casadi::Sparsity::operator*;
%rename(minus) casadi::Slice::operator-;
%rename(times) casadi::Slice::operator*;
%ignore operator==;
%ignore operator!=;
%ignore casadi::DaeBuilder::operator();
// Object.finalize() is called by the garbage collector
%rename(finalize_) casadi::Callback::finalize;
// Slice(long, boolean) would clash with the proxy constructor (pointer, ownership)
%ignore casadi::Slice::Slice(casadi_int i, bool ind1);
%extend casadi::Slice {
  Slice(casadi_int i) { return new casadi::Slice(i); }
}
// Interface-internal macros
%ignore CASADI_INT_TYPE_STR; %ignore MNAME; %ignore L_INT; %ignore L_BOOL; %ignore LL; %ignore LR;
%ignore L_DICT; %ignore L_DOUBLE; %ignore L_STR; %ignore LABEL; %ignore SWIG_IND1;
#endif // SWIGJAVA

#ifdef SWIGMATLAB
%rename(uminus) operator-;
// Pattern intersection is elementwise; * and mtimes give the pattern of the matrix product
%rename(times) casadi::Sparsity::operator*;
// inv is the matrix inverse in MATLAB; elementwise inverse is 1./x
%ignore casadi::MX::inv() const;
// Would hide the einstein method of GenericMatrixCommon for MX arguments
%ignore casadi::MX::einstein;
%feature("varargin","1") casadi_vertcat;
%feature("varargin","1") casadi_horzcat;
%feature("varargin","1") casadi_diagcat;
%feature("varargin","1") casadi_veccat;
%feature("optionalunpack","1") size;

// Raise an error if "this" not correct
%typemap(check, noblock=1) SWIGTYPE *self %{
if (!$1) {
  SWIG_Error(SWIG_RuntimeError, "Invalid 'self' object");
  SWIG_fail;
 }
%}

// Workarounds, pending proper fix
%rename(nonzero) __nonzero__;
%rename(hash) __hash__;

%rename(rem) casadi_mod;
#endif // SWIGMATLAB

#ifdef SWIGPYTHON
%ignore casadi_mod;
%rename(__bool__) __nonzero__;

%pythoncode %{
class NZproxy:
  def __init__(self,matrix):
    self.matrix = matrix

  def __getitem__(self,s):
    return self.matrix.get_nz(False, s)

  def __setitem__(self,s,val):
    return self.matrix.set_nz(val, False, s)

  def __len__(self):
    return self.matrix.nnz()

  def __iter__(self):
    for i in range(len(self)):
      yield self[i]

%}

// NZproxy (%pythoncode), generic in the owning matrix type
#ifdef SWIG_STUBS_ENABLED
%stubcode %{class NZproxy(Generic[_TNZ]):
%}
%stub_method(__init__,    None, matrix: _TNZ)
%stub_overload_method_selftyped(__getitem__, MX, NZproxy[MX], s: _MIndex | MX)
%stub_overload_method(__getitem__, _TNZ,  s: _MIndex)
%stub_overload_method_selftyped(__setitem__, None, NZproxy[MX], s: _MIndex | MX, val: float | MX | Sequence[float])
%stub_overload_method(__setitem__, None, s: _MIndex, val: float | _TNZ | Sequence[float])
%stub_method0(__len__,  int)
%stub_method0(__iter__, Iterator[_T])
#endif

%define %python_array_wrappers(arraypriority, arraydtype)
%pythoncode %{

  __array_priority__ = arraypriority

  def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
      return _numpy_ufunc_dispatch(self, ufunc, method, inputs, kwargs)

  def __array_function__(self, func, types, args, kwargs):
      return _numpy_array_function_dispatch(self, func, types, args, kwargs)

  def __array__(self, *args, **kwargs):
      import numpy as _n
      try:
          arr = self.full()
      except Exception:
          if self.is_scalar(True):
              # Box symbolic scalars in an object array (#2743).
              E = _n.empty((), dtype=object)
              E[()] = self
              return E
          raise TypeError(
              "Implicit conversion of symbolic CasADi type to numeric matrix not supported.\n"
              "This may occur when you pass a CasADi object to a numpy function.\n"
              "Use an equivalent CasADi function instead of that numpy function.")
      dtype = kwargs.get("dtype")
      if dtype is not None and dtype is not _n.double:
          return _n.array(arr, dtype=dtype)
      return arr

%}
// __array__ makes DM/SX/MX numpy ArrayLike
%stub_method(__array__, %arg(NDArray[arraydtype]), *args: Any, **kwargs: Any)
#ifdef SWIG_STUBS_ENABLED
%stubcode %{    __array_priority__: float
%}
#endif // SWIG_STUBS_ENABLED
%stub_method(__array_ufunc__, Any, ufunc: Any, method: builtins.str, *inputs: Any, **kwargs: Any)
%stub_method(__array_function__, Any, func: Any, types: Any, args: Any, kwargs: Any)
%enddef
#endif // SWIGPYTHON

#ifdef SWIGJAVA
// Indexing (0-based) as Kotlin's get/set operators: x[i], x[i, j], x[i] = v
%define %java_index_pair(Type, R, C)
    Type get(const R& rr, const C& cc) const {
      Type m;
      $self->get(m, false, rr, cc);
      return m;
    }
    void set(const R& rr, const C& cc, const Type& m) { $self->set(m, false, rr, cc); }
%enddef
#endif // SWIGJAVA

#ifdef SWIGMATLAB
%{
  namespace casadi {
    /// Helper function: Convert ':' to Slice
    inline Slice char2Slice(char ch) {
      casadi_assert_dev(ch==':');
      return Slice();
    }
  } // namespace casadi
%}

#endif // SWIGMATLAB

// MATLAB two-index paren overloads, with ':' passed as char
%define %_paren2(Type, RT, R, CT, C)
    const Type paren(RT rr, CT cc) const { Type m; $self->get(m, true, R, C); return m;}
%enddef
%define %_paren_asgn2(Type, RT, R, CT, C)
    void paren_asgn(const Type& m, RT rr, CT cc) { $self->set(m, true, R, C);}
%enddef
%define %_paren_combos(Macro, Type, I)
    Macro(Type, char, casadi::char2Slice(rr), I, cc)
    Macro(Type, I, rr, char, casadi::char2Slice(cc))
    Macro(Type, I, rr, I, cc)
%enddef

%define %matrix_helpers(Type, xIndex)
#ifdef SWIGPYTHON
%pythoncode %{
    @property
    def shape(self):
        return (self.size1(),self.size2())

    def reshape(self,arg):
        return _casadi.reshape(self,arg)

    @property
    def T(self):
        return _casadi.transpose(self)

    def __getitem__(self, s):
          if isinstance(s, tuple) and len(s)==2:
            if s[1] is None: raise TypeError("Cannot slice with None")
            return self.get(False, s[0], s[1])
          return self.get(False, s)

    def __iter__(self):
      raise TypeError("""CasADi matrices are not iterable by design.
                      Did you mean to iterate over m.nz, with m IM/DM/SX?
                      Did you mean to iterate over horzsplit(m,1)/vertsplit(m,1) with m IM/DM/SX/MX?
                      Did you expect numpy behaviour (looping over rows)? Use ca.array().
                      """)

    def __setitem__(self,s,val):
          if isinstance(s,tuple) and len(s)==2:
            return self.set(val, False, s[0], s[1])
          return self.set(val, False, s)

    @property
    def nz(self):
      return NZproxy(self)

%}
%stub_attr(shape, tuple[int, int])
%stub_attr(T,     Self)
%stub_attr(nz,    NZproxy[Self])
%stub_method(reshape,     Self, arg: tuple[int, int] | int | Sparsity)
// Indexing; MX adds MX indices in its own %extend
%stub_method(__getitem__, Self, s: xIndex | tuple[xIndex, xIndex])
%stub_method(__setitem__, None, s: xIndex | tuple[xIndex, xIndex], val: _DM | SX | MX)
%stub_method0(__iter__, %arg(Iterator[Self]))
#elif defined(SWIGMATLAB)
    // Get a submatrix (index-1)
    const Type paren(char rr) const {
      casadi_assert_dev(rr==':');
      return vec(*$self);
    }
    const Type paren(const Matrix<casadi_int>& rr) const {
      Type m;
      $self->get(m, true, rr);
      return m;
    }
    const Type paren(const Sparsity& sp) const {
      Type m;
      $self->get(m, true, sp);
      return m;
    }
    %_paren2(Type, char, casadi::char2Slice(rr), char, casadi::char2Slice(cc))
    %_paren_combos(%_paren2, Type, const Matrix<casadi_int>&)

    // Set a submatrix (index-1)
    void paren_asgn(const Type& m, char rr) {
      casadi_assert_dev(rr==':');
      $self->set(m, false, casadi::IM(casadi::range($self->numel())));
    }
    void paren_asgn(const Type& m, const Matrix<casadi_int>& rr) { $self->set(m, true, rr);}
    void paren_asgn(const Type& m, const Sparsity& sp) { $self->set(m, true, sp);}
    %_paren_asgn2(Type, char, casadi::char2Slice(rr), char, casadi::char2Slice(cc))
    %_paren_combos(%_paren_asgn2, Type, const Matrix<casadi_int>&)

    // Get nonzeros (index-1)
    const Type brace(char rr) const { Type m; $self->get_nz(m, true, casadi::char2Slice(rr)); return m;}
    const Type brace(const Matrix<casadi_int>& rr) const { Type m; $self->get_nz(m, true, rr); return m;}

    // Set nonzeros (index-1)
    void setbrace(const Type& m, char rr) { $self->set_nz(m, true, casadi::char2Slice(rr));}
    void setbrace(const Type& m, const Matrix<casadi_int>& rr) { $self->set_nz(m, true, rr);}

    // 'end' function (needed for end syntax in MATLAB)
    inline casadi_int end(casadi_int i, casadi_int n) const {
      return n==1 ? $self->numel() : i==1 ? $self->size1() : $self->size2();
    }


    // Needed for brace syntax to access nonzeros
    casadi_int numel(casadi_int k) const {
      return 1;
    }

    // Needed for brace syntax to access nonzeros
    casadi_int numel(char rr) const {
      casadi_assert_dev(rr==':');
      return 1;
    }

    // Needed for brace syntax to access nonzeros
    casadi_int numel(const std::vector<casadi_int> &k) const {
      return 1;
    }

    // Needed because original numel call gets hidden by the above extend overloads
    casadi_int numel() const {
      return $self->numel();
    }


    // Transpose using the A' syntax in addition to A.'
    Type ctranspose() const { return $self->T();}

#elif defined(SWIGJAVA)
    Type get(const casadi::Slice& rr) const {
      Type m;
      $self->get(m, false, rr);
      return m;
    }
    Type get(const casadi::Matrix<casadi_int>& rr) const {
      Type m;
      $self->get(m, false, rr);
      return m;
    }
    Type get(const casadi::Sparsity& sp) const {
      Type m;
      $self->get(m, false, sp);
      return m;
    }
    void set(const casadi::Slice& rr, const Type& m) { $self->set(m, false, rr); }
    void set(const casadi::Matrix<casadi_int>& rr, const Type& m) { $self->set(m, false, rr); }
    void set(const casadi::Sparsity& sp, const Type& m) { $self->set(m, false, sp); }
    %java_index_pair(Type, casadi::Slice, casadi::Slice)
    %java_index_pair(Type, casadi::Slice, casadi::Matrix<casadi_int>)
    %java_index_pair(Type, casadi::Matrix<casadi_int>, casadi::Slice)
    %java_index_pair(Type, casadi::Matrix<casadi_int>, casadi::Matrix<casadi_int>)
#endif
%enddef

#if defined(SWIGMATLAB) || defined(SWIGOCTAVE)
// Abstract mixins avoid the Octave 10+ SwigRef diamond; before their %include
%feature("abstract") casadi::PrintableCommon;
%feature("abstract") casadi::GenericExpressionCommon;
%feature("abstract") casadi::GenericMatrixCommon;
%feature("abstract") casadi::SparsityInterfaceCommon;
%feature("abstract") casadi::MatrixCommon;
// Applied to the template, not to GenSharedObject
%feature("abstract") casadi::GenericShared<casadi::SharedObject, casadi::SharedObjectInternal>;
%feature("abstract") casadi::SharedObject;
%feature("abstract") casadi::IndexAbstraction;
#endif

%include <casadi/core/printable.hpp>

namespace casadi{
%extend PrintableCommon {
#ifdef SWIGPYTHON
  %pythoncode %{
    def __str__(self): return self.str()
    def repr(self): return self.type_name() + '(' + self.str() + ')'
    # Value-copy semantics for all casadi types (DM/SX/MX/Sparsity/Function/
    # ...): copy via the copy constructor.  This base injects __copy__ /
    # __deepcopy__ into every subclass, replacing the old global-`object`
    # shadow hack.
    def __copy__(self): return self.__class__(self)
    def __deepcopy__(self, memo=None): return self.__class__(self)
  %}
  %stub_method0(repr, builtins.str)
  %stub_method0(__copy__, Self)
  %stub_method(__deepcopy__, Self, memo: dict[int, Any] | None = None)
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
  %matlabcode %{
    function s = repr(self)
      s = [self.type_name() '(' self.str() ')'];
    end
  %}
#endif // SWIGMATLAB
}
} // namespace casadi

%include <casadi/core/generic_shared.hpp>

%template(GenSharedObject) casadi::GenericShared<casadi::SharedObject, casadi::SharedObjectInternal>;
%template(GenWeakRef) casadi::GenericWeakRef<casadi::SharedObject, casadi::SharedObjectInternal>;

%include <casadi/core/shared_object.hpp>
%include <casadi/core/casadi_misc.hpp>
%include <casadi/core/casadi_common.hpp>
%include <casadi/core/generic_type.hpp>
%include <casadi/core/calculus.hpp>
%include <casadi/core/sparsity_interface.hpp>
%include <casadi/core/sparsity.hpp>

// Logic for pickling
#ifdef SWIGPYTHON
namespace casadi{
%extend Sparsity {
  %pythoncode %{
    def __setstate__(self, state):
        self.__init__(Sparsity.deserialize(state["serialization"]))

    def __getstate__(self):
        return {"serialization": self.serialize()}
  %}
  %stub_method0(__getstate__, %arg(dict[builtins.str, Any]))
  %stub_method(__setstate__, None, state: Mapping[builtins.str, Any])
  %stub_method(__array__, %arg(NDArray[Any]), *args: Any, **kwargs: Any)
  // Getter-backed attributes without a typemap doc
  %stub_attr(T, Sparsity)
  %stub_attr(shape, tuple[int, int]) // TYPE is variadic; inner comma is fine here
}
%extend Matrix<SXElem> {
  %pythoncode %{
    def __setstate__(self, state):
      ctx = _current_unpickle_context()
      if not ctx:
        raise Exception("Cannot unpickle SX objects without a casadi context. " +
          "Use something like:\n"+
          "with ca.global_unpickle_context(): \n"+
          "  f_ref = pickle.load(open(filename,'rb'))")
      ctx.decode(state)
      self.__init__(ctx.unpack())

    def __getstate__(self):
      ctx = _current_pickle_context()
      if not ctx:
        raise Exception("Cannot pickle SX objects without a casadi context. " +
          "Use something like:\n"+
          "with ca.global_pickle_context(): \n"+
          "  pickle.dump(f,open(filename,'wb'))")
      ctx.pack(self)
      return ctx.encode()
  %}
  %stub_method0(__getstate__, Any)
  %stub_method(__setstate__, None, state: Any)
}
%extend MX {
  %pythoncode %{
    def __setstate__(self, state):
      ctx = _current_unpickle_context()
      if not ctx:
        raise Exception("Cannot unpickle MX objects without a casadi context. " +
          "Use something like:\n"+
          "with ca.global_unpickle_context(): \n"+
          "  f_ref = pickle.load(open(filename,'rb'))")
      ctx.decode(state)
      self.__init__(ctx.unpack())

    def __getstate__(self):
      ctx = _current_pickle_context()
      if not ctx:
        raise Exception("Cannot pickle MX objects without a casadi context. " +
          "Use something like:\n"+
          "with ca.global_pickle_context(): \n"+
          "  pickle.dump(f,open(filename,'wb'))")
      ctx.pack(self)
      return ctx.encode()
  %}
  %stub_method0(__getstate__, Any)
  %stub_method(__setstate__, None, state: Any)
}

} // namespace casadi
#endif // SWIGPYTHON

%include <casadi/core/slice.hpp>


%include <casadi/core/generic_matrix.hpp>

%template(GenDM)        casadi::GenericMatrix<casadi::Matrix<double> >;
%template(GenSX)             casadi::GenericMatrix<casadi::Matrix<casadi::SXElem> >;
%template(GenMX)             casadi::GenericMatrix<casadi::MX>;

%include <casadi/core/generic_expression.hpp>

// Interface-internal flags, not module constants
%ignore IS_GLOBAL; %ignore IS_MEMBER; %ignore IS_SPARSITY; %ignore IS_DMATRIX;
%ignore IS_IMATRIX; %ignore IS_SX; %ignore IS_MX; %ignore IS_DOUBLE; %ignore FLAG;
#define IS_GLOBAL   0x1
#define IS_MEMBER   0x10
#define IS_SPARSITY 0x100
#define IS_DMATRIX  0x1000
#define IS_IMATRIX  0x10000
#define IS_SX       0x100000
#define IS_MX       0x1000000
#define IS_DOUBLE   0x10000000

%define SPARSITY_INTERFACE_FUN_BASE(DECL, FLAG, M)
#if FLAG & IS_MEMBER

DECL M casadi_horzcat(const std::vector< M > &v) { return horzcat(v); }
DECL M casadi_vertcat(const std::vector< M > &v) { return vertcat(v); }
DECL std::vector< M > casadi_horzsplit(const M& v, const std::vector<casadi_int>& offset) { return horzsplit(v, offset); }
DECL std::vector< M > casadi_horzsplit(const M& v, casadi_int incr=1) { return horzsplit(v, incr); }
DECL std::vector< M > casadi_horzsplit_n(const M& v, casadi_int n) { return horzsplit_n(v, n); }
DECL std::vector< M > casadi_vertsplit(const M& v, const std::vector<casadi_int>& offset) { return vertsplit(v, offset); }
DECL std::vector<casadi_int > casadi_offset(const std::vector< M > &v, bool vert=true) { return offset(v, vert); }
DECL std::vector< M > casadi_vertsplit(const M& v, casadi_int incr=1) { return vertsplit(v, incr); }
DECL std::vector< M > casadi_vertsplit_n(const M& v, casadi_int n) { return vertsplit_n(v, n); }
DECL M casadi_blockcat(const M& A, const M& B, const M& C, const M& D) { return vertcat(horzcat(A, B), horzcat(C, D)); }
DECL std::vector< std::vector< M > > casadi_blocksplit(const M& x, const std::vector<casadi_int>& vert_offset, const std::vector<casadi_int>& horz_offset) { return blocksplit(x, vert_offset, horz_offset); }
DECL std::vector< std::vector< M > > casadi_blocksplit(const M& x, casadi_int vert_incr=1, casadi_int horz_incr=1) { return blocksplit(x, vert_incr, horz_incr); }
DECL M casadi_diagcat(const std::vector< M > &A) { return diagcat(A); }
DECL std::vector< M > casadi_diagsplit(const M& x, const std::vector<casadi_int>& output_offset1, const std::vector<casadi_int>& output_offset2) { return diagsplit(x, output_offset1, output_offset2); }
DECL std::vector< M > casadi_diagsplit(const M& x, const std::vector<casadi_int>& output_offset) { return diagsplit(x, output_offset); }
DECL std::vector< M > casadi_diagsplit(const M& x, casadi_int incr=1) { return diagsplit(x, incr); }
DECL std::vector< M > casadi_diagsplit(const M& x, casadi_int incr1, casadi_int incr2) { return diagsplit(x, incr1, incr2); }
DECL M casadi_veccat(const std::vector< M >& x) { return veccat(x); }
DECL M casadi_mtimes(const M& x, const M& y, const std::string& blas = "reference") { return mtimes(x, y, blas); }
DECL M casadi_mtimes(const std::vector< M > &args) { return mtimes(args); }
DECL M casadi_mtimes(const std::vector< M > &args, const std::string& blas) { return mtimes(args, blas); }
DECL M casadi_mac(const M& X, const M& Y, const M& Z, const std::string& blas = "reference") { return mac(X, Y, Z, blas); }
DECL M casadi_transpose(const M& X) { return X.T(); }
DECL M casadi_vec(const M& a) { return vec(a); }
DECL M casadi_reshape(const M& a, casadi_int nrow, casadi_int ncol) { return reshape(a, nrow, ncol); }
DECL M casadi_reshape(const M& a, std::pair<casadi_int, casadi_int> rc) { return reshape(a, rc.first, rc.second); }
DECL M casadi_reshape(const M& a, const Sparsity& sp) { return reshape(a, sp); }
DECL M casadi_sparsity_cast(const M& a, const Sparsity& sp) { return sparsity_cast(a, sp); }
DECL casadi_int casadi_sprank(const M& A) { return sprank(A); }
DECL casadi_int casadi_norm_0_mul(const M& x, const M& y) { return norm_0_mul(x, y); }
DECL M casadi_triu(const M& a, bool includeDiagonal=true) { return triu(a, includeDiagonal); }
DECL M casadi_tril(const M& a, bool includeDiagonal=true) { return tril(a, includeDiagonal); }
DECL M casadi_kron(const M& a, const M& b) { return kron(a, b); }
DECL M casadi_repmat(const M& A, casadi_int n, casadi_int m=1) { return repmat(A, n, m); }
DECL M casadi_repmat(const M& A, const std::pair<casadi_int, casadi_int>& rc) { return repmat(A, rc.first, rc.second); }
DECL M casadi_sum2(const M& x) { return sum2(x); }
DECL M casadi_sum1(const M& x) { return sum1(x); }
#endif
#if defined(SWIGMATLAB) && (FLAG & IS_GLOBAL)
// Package functions for calls without a casadi object first, e.g. offset({...}) or diagcat()
DECL std::vector<casadi_int > casadi_offset(const std::vector< M > &v, bool vert=true) { return offset(v, vert); }
DECL M casadi_blockcat(const M& A, const M& B, const M& C, const M& D) { return vertcat(horzcat(A, B), horzcat(C, D)); }
DECL M casadi_diagcat(const std::vector< M > &A) { return diagcat(A); }
DECL M casadi_veccat(const std::vector< M >& x) { return veccat(x); }
#endif
#if defined(SWIGMATLAB) && (FLAG & IS_MEMBER)
// [] for a dimension: deduce it, as reshape(x, [], 1) in MATLAB
DECL M casadi_reshape(const M& a, const std::vector<casadi_int>& nrow, casadi_int ncol) {
  casadi_assert(nrow.size()<=1, "reshape: expected a scalar or [] for the number of rows");
  return reshape(a, nrow.empty() ? -1 : nrow[0], ncol);
}
DECL M casadi_reshape(const M& a, casadi_int nrow, const std::vector<casadi_int>& ncol) {
  casadi_assert(ncol.size()<=1, "reshape: expected a scalar or [] for the number of columns");
  return reshape(a, nrow, ncol.empty() ? -1 : ncol[0]);
}
#endif
%enddef

%define SPARSITY_INTERFACE_ALL(DECL, FLAG)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_SPARSITY), Sparsity)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_MX), MX)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_DMATRIX), Matrix<double>)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_SX), Matrix<SXElem>)
%enddef

%define SPARSITY_INTERFACE_FUN(DECL, FLAG, M)
SPARSITY_INTERFACE_FUN_BASE(DECL, FLAG, M)
#if FLAG & IS_MEMBER
#ifdef SWIGMATLAB
DECL casadi_int casadi_length(const M &v) {
  if (v.is_empty()) return 0;
  return std::max(v.size1(), v.size2());
}
DECL M casadi_sum(const M& x, casadi_int dim) {
  if (dim==1) return sum1(x);
  if (dim==2) return sum2(x);
  casadi_error(
    "Expected sum(A,1), sum(A,2), sum(A,\"all\") got " + casadi::str(dim) + " instead.");
}
DECL M casadi_sum(const M& x, const std::string& dim) {
  casadi_assert(dim=="all",
    "Expected sum(...,'all'), got '" + dim + "' instead.");
  return sum(x);
}
DECL M casadi_sum(const M& x) {
  if (x.is_vector()) return sum(x);
  return sum1(x);
}
#else
DECL M casadi_sum(const M& x, casadi_int dim) {
  if (dim==0) return sum1(x);
  if (dim==1) return sum2(x);
  casadi_error(
    "Expected sum(A,0), sum(A,1), sum(A,\"all\") got " + casadi::str(dim) + " instead.");
}
DECL M casadi_sum(const M& x) { return sum(x); }
#endif
#endif
%enddef

%define GENERIC_MATRIX_FUN(DECL, FLAG, M)
#if FLAG & IS_MEMBER
DECL M casadi_mpower(const M& x, const M& n) { return mpower(x, n); }
DECL M casadi_mrdivide(const M& x, const M& y) { return mrdivide(x, y); }
DECL M casadi_mldivide(const M& x, const M& y) { return mldivide(x, y); }
DECL std::vector< M > casadi_symvar(const M& x) { return symvar(x); }
DECL M casadi_transform(const M& x, const casadi::Dict& opts = casadi::Dict()) { return transform(x, opts); }
DECL M casadi_transform(const M& x, const std::vector< std::vector< casadi::GenericType > >& passes, const casadi::Dict& opts = casadi::Dict()) { return transform(x, passes, opts); }
DECL std::vector< M > casadi_transform(const std::vector< M >& x, const casadi::Dict& opts = casadi::Dict()) { return transform(x, opts); }
DECL std::vector< M > casadi_transform(const std::vector< M >& x, const std::vector< std::vector< casadi::GenericType > >& passes, const casadi::Dict& opts = casadi::Dict()) { return transform(x, passes, opts); }
DECL M casadi_bilin(const M& A, const M& x, const M& y) { return bilin(A, x, y); }
DECL M casadi_bilin(const M& A, const M& x) { return bilin(A, x); }
DECL M casadi_rank1(const M& A, const M& alpha, const M& x, const M& y) { return rank1(A, alpha, x, y); }
DECL M casadi_sumsqr(const M& X) { return sumsqr(X); }
DECL M casadi_linspace(const M& a, const M& b, casadi_int nsteps) { return linspace(a, b, nsteps); }
DECL M casadi_logsumexp(const M& a) { return logsumexp(a); }
DECL M casadi_logsumexp(const M& a, const M& margin) { return logsumexp(a, margin); }
DECL M casadi_interp1d(const std::vector<double>& x, const M&v, const std::vector<double>& xq, const std::string& mode="linear", bool equidistant=false) { return interp1d(x, v, xq, mode, equidistant); }
DECL M casadi_soc(const M& x, const M& y) { return soc(x, y); }
DECL M casadi_cross(const M& a, const M& b, casadi_int dim = -1) { return cross(a, b, dim); }
DECL M casadi_skew(const M& a) { return skew(a); }
DECL M casadi_inv_skew(const M& a) { return inv_skew(a); }
DECL M casadi_det(const M& A) { return det(A); }
DECL M casadi_det(const M& A, const std::string& lsolver, const casadi::Dict& opts = casadi::Dict()) { return det(A, lsolver, opts); }
DECL M casadi_inv_minor(const M& A) { return inv_minor(A); }
DECL M casadi_inv(const M& A) { return inv(A); }
DECL M casadi_inv(const M& A, const std::string& lsolver, const casadi::Dict& opts = casadi::Dict()) { return inv(A, lsolver, opts); }
DECL M casadi_trace(const M& a) { return trace(a); }
DECL M casadi_tril2symm(const M& a) { return tril2symm(a); }
DECL M casadi_triu2symm(const M& a) { return triu2symm(a); }
DECL M casadi_norm_fro(const M& x) { return norm_fro(x); }
DECL M casadi_norm_2(const M& x) { return norm_2(x); }
DECL M casadi_norm_1(const M& x) { return norm_1(x); }
DECL M casadi_norm_inf(const M& x) { return norm_inf(x); }
DECL M casadi_dot(const M& x, const M& y) { return dot(x, y); }
DECL M casadi_nullspace(const M& A) { return nullspace(A); }
DECL M casadi_polyval(const M& p, const M& x) { return polyval(p, x); }
DECL M casadi_diag(const M& A) { return diag(A); }
DECL M casadi_unite(const M& A, const M& B) { return unite(A, B); }
DECL M casadi_densify(const M& x) { return densify(x); }
DECL M casadi_densify(const M& x, const M& val) { return densify(x, val); }
DECL M casadi_project(const M& A, const Sparsity& sp, bool intersect=false) { return project(A, sp, intersect); }
DECL M casadi_if_else(const M& cond, const M& if_true, const M& if_false, bool short_circuit=false) { return if_else(cond, if_true, if_false, short_circuit); }
DECL M casadi_conditional(const M& ind, const std::vector< M > &x, const M& x_default, bool short_circuit=false) { return conditional(ind, x, x_default, short_circuit); }
DECL bool casadi_depends_on(const M& f, const M& arg) { return depends_on(f, arg); }
DECL bool casadi_contains(const std::vector<M>& v, const M& n) { return contains(v, n); }
DECL bool casadi_contains_all(const std::vector<M>& v, const std::vector<M>& n) { return contains_all(v, n); }
DECL bool casadi_contains_any(const std::vector<M>& v, const std::vector<M>& n) { return contains_any(v, n); }
DECL M casadi_solve(const M& A, const M& b) { return solve(A, b); }
DECL M casadi_solve(const M& A, const M& b, const std::string& lsolver, const casadi::Dict& opts = casadi::Dict()) { return solve(A, b, lsolver, opts); }
DECL M casadi_pinv(const M& A) { return pinv(A); }
DECL M casadi_pinv(const M& A, const std::string& lsolver, const casadi::Dict& opts = casadi::Dict()) { return pinv(A, lsolver, opts); }
DECL M casadi_expm_const(const M& A, const M& t) { return expm_const(A, t); }
DECL M casadi_expm(const M& A) { return expm(A); }
DECL M casadi_jacobian(const M &ex, const M &arg, const Dict& opts=Dict()) { return jacobian(ex, arg, opts); }
DECL M casadi_jtimes(const M& ex, const M& arg, const M& v, bool tr=false, const Dict& opts=Dict()) { return jtimes(ex, arg, v, tr); }
DECL M casadi_linearize(const M& f, const M& x, const M& x0) { return linearize(f, x, x0); }
DECL std::vector<bool> casadi_which_depends(const M& expr, const M& var, casadi_int order=1, bool tr=false) { return which_depends(expr, var, order, tr); }
DECL Sparsity casadi_jacobian_sparsity(const M& f, const M& x) { return jacobian_sparsity(f, x); }
DECL bool casadi_is_linear(const M& expr, const M& var) { return is_linear(expr, var); }
DECL bool casadi_is_quadratic(const M& expr, const M& var) { return is_quadratic(expr, var); }
DECL M casadi_gradient(const M &ex, const M &arg, const Dict& opts=Dict()) { return gradient(ex, arg, opts); }
DECL M casadi_tangent(const M &ex, const M &arg, const Dict& opts=Dict()) { return tangent(ex, arg, opts); }
DECL M casadi_hessian(const M& ex, const M& arg, M& OUTPUT1, const casadi::Dict& opts = casadi::Dict()) { return hessian(ex, arg, OUTPUT1, opts); }
DECL void casadi_quadratic_coeff(const M& ex, const M& arg, M& OUTPUT1, M& OUTPUT2, M& OUTPUT3, bool check=true) { quadratic_coeff(ex, arg, OUTPUT1, OUTPUT2, OUTPUT3, check); }
DECL void casadi_linear_coeff(const M& ex, const M& arg, M& OUTPUT1, M& OUTPUT2, bool check=true) { linear_coeff(ex, arg, OUTPUT1, OUTPUT2, check); }
DECL casadi_int casadi_n_nodes(const M& A) { return n_nodes(A); }
DECL std::string casadi_print_operator(const M& xb, const std::vector<std::string>& args) { return print_operator(xb, args); }
DECL M casadi_repsum(const M& A, casadi_int n, casadi_int m=1) { return repsum(A, n, m); }
DECL M casadi_diff(const M& A, casadi_int n=1, casadi_index axis=-1) { return diff(A, n, axis); }
DECL M casadi_cumsum(const M& A, casadi_index axis=-1) { return cumsum(A, axis); }
DECL M casadi_einstein(const M& A, const M& B, const M& C, const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b, const std::vector<casadi_int>& dim_c, const std::vector<casadi_int>& a, const std::vector<casadi_int>& b, const std::vector<casadi_int>& c) { return einstein(A, B, C, dim_a, dim_b, dim_c, a, b, c); }
DECL M casadi_einstein(const M& A, const M& B, const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b, const std::vector<casadi_int>& dim_c, const std::vector<casadi_int>& a, const std::vector<casadi_int>& b, const std::vector<casadi_int>& c) { return einstein(A, B, dim_a, dim_b, dim_c, a, b, c); }
DECL M casadi_mmin(const M& x) { return mmin(x); }
DECL M casadi_mmax(const M& x) { return mmax(x); }
DECL casadi::DM casadi_evalf(const M& x) { return evalf(x); }
DECL void casadi_separate_linear(const M &expr, const M &sym_lin, const M &sym_const, M& OUTPUT1, M& OUTPUT2, M& OUTPUT3) { separate_linear(expr, sym_lin, sym_const, OUTPUT1, OUTPUT2, OUTPUT3); }
DECL void casadi_separate_linear(const M &expr, const std::vector<M> &sym_lin, const std::vector<M> &sym_const, M& OUTPUT1, M& OUTPUT2, M& OUTPUT3) { separate_linear(expr, sym_lin, sym_const, OUTPUT1, OUTPUT2, OUTPUT3); }
#endif // FLAG & IS_MEMBER

#if FLAG & IS_GLOBAL
DECL std::vector<M> casadi_cse(const std::vector<M>& e) { return cse(e); }
DECL M casadi_cse(const M& e) { return cse(e); }
DECL void casadi_extract_parametric(const M &expr, const M& par, M& OUTPUT1, std::vector<M>& OUTPUT2, std::vector<M>& OUTPUT3, const Dict& opts=Dict()) { extract_parametric(expr, par, OUTPUT1, OUTPUT2, OUTPUT3, opts); }
DECL void casadi_extract_parametric(const M &expr, const std::vector<M>& par, M& OUTPUT1, std::vector<M>& OUTPUT2, std::vector<M>& OUTPUT3, const Dict& opts=Dict()) { extract_parametric(expr, par, OUTPUT1, OUTPUT2, OUTPUT3, opts); }
DECL void casadi_extract_parametric(const std::vector<M> &expr, const M& par, std::vector<M>& OUTPUT1, std::vector<M>& OUTPUT2, std::vector<M>& OUTPUT3, const Dict& opts=Dict()) { extract_parametric(expr, par, OUTPUT1, OUTPUT2, OUTPUT3, opts); }
DECL void casadi_extract_parametric(const std::vector<M> &expr, const std::vector<M>& par, std::vector<M>& OUTPUT1, std::vector<M>& OUTPUT2, std::vector<M>& OUTPUT3, const Dict& opts=Dict()) { extract_parametric(expr, par, OUTPUT1, OUTPUT2, OUTPUT3, opts); }
DECL std::vector<std::vector< M > > casadi_forward(const std::vector< M > &ex, const std::vector< M > &arg, const std::vector<std::vector< M > > &v, const Dict& opts = Dict()) { return forward(ex, arg, v, opts); }
DECL std::vector<std::vector< M > > casadi_reverse(const std::vector< M > &ex, const std::vector< M > &arg, const std::vector<std::vector< M > > &v, const Dict& opts = Dict()) { return reverse(ex, arg, v, opts); }

#if FLAG & (IS_SX | IS_MX)
DECL std::string casadi_export_graph(const std::vector<M>& expressions, const Dict& opts=Dict()) { return export_graph(expressions, opts); }
DECL std::string casadi_export_graph(const M& expression, const Dict& opts=Dict()) { return export_graph(expression, opts); }
DECL void casadi_export_graph(const std::vector<M>& expressions, const std::string& fname, const Dict& opts=Dict()) { export_graph(expressions, fname, opts); }
DECL void casadi_export_graph(const M& expression, const std::string& fname, const Dict& opts=Dict()) { export_graph(expression, fname, opts); }
#endif

DECL M casadi_substitute(const M& ex, const M& v, const M& vdef) { return substitute(ex, v, vdef); }
DECL std::vector< M > casadi_substitute(const std::vector< M >& ex, const std::vector< M >& v, const std::vector< M >& vdef) { return substitute(ex, v, vdef); }
DECL void casadi_substitute_inplace(const std::vector< M >& v, std::vector< M >& INOUT1, std::vector< M >& INOUT2, bool reverse=false) { return substitute_inplace(v, INOUT1, INOUT2, reverse); }
DECL void casadi_extract(const std::vector< M >& ex, std::vector< M >& OUTPUT1, std::vector< M >& OUTPUT2, std::vector< M >& OUTPUT3, const Dict& opts = Dict()) { OUTPUT1 = ex; extract(OUTPUT1, OUTPUT2, OUTPUT3, opts); }
DECL void casadi_shared(const std::vector< M >& ex, std::vector< M >& OUTPUT1, std::vector< M >& OUTPUT2, std::vector< M >& OUTPUT3, const std::string& v_prefix="v_", const std::string& v_suffix="") { OUTPUT1 = ex; shared(OUTPUT1, OUTPUT2, OUTPUT3, v_prefix, v_suffix); }
DECL M casadi_blockcat(const std::vector< std::vector< M > > &v) { return blockcat(v); }
#endif // FLAG & IS_GLOBAL
%enddef

%define GENERIC_MATRIX_ALL(DECL, FLAG)
GENERIC_MATRIX_FUN(DECL, (FLAG | IS_MX), MX)
GENERIC_MATRIX_FUN(DECL, (FLAG | IS_DMATRIX), Matrix<double>)
GENERIC_MATRIX_FUN(DECL, (FLAG | IS_SX), Matrix<SXElem>)
%enddef

%define GENERIC_EXPRESSION_FUN(DECL, FLAG, M)
#if FLAG & IS_MEMBER
DECL M casadi_plus(const M& x, const M& y) { return x+y; }
DECL M casadi_minus(const M& x, const M& y) { return x-y; }
DECL M casadi_times(const M& x, const M& y) { return x*y; }
DECL M casadi_rdivide(const M& x, const M& y) { return x/y; }
DECL M casadi_ldivide(const M& x, const M& y) { return y/x; }
DECL M casadi_lt(const M& x, const M& y) { return x<y; }
DECL M casadi_le(const M& x, const M& y) { return x<=y; }
DECL M casadi_gt(const M& x, const M& y) { return x>y; }
DECL M casadi_ge(const M& x, const M& y) { return x>=y; }
DECL M casadi_eq(const M& x, const M& y) { return x==y; }
DECL M casadi_ne(const M& x, const M& y) { return x!=y; }
DECL M casadi_and(const M& x, const M& y) { return x&&y; }
DECL M casadi_or(const M& x, const M& y) { return x||y; }
DECL M casadi_not(const M& x) { return !x; }
DECL M casadi_abs(const M& x) { return fabs(x); }
DECL M casadi_sqrt(const M& x) { return sqrt(x); }
DECL M casadi_sin(const M& x) { return sin(x); }
DECL M casadi_cos(const M& x) { return cos(x); }
DECL M casadi_tan(const M& x) { return tan(x); }
DECL M casadi_atan(const M& x) { return atan(x); }
DECL M casadi_asin(const M& x) { return asin(x); }
DECL M casadi_acos(const M& x) { return acos(x); }
DECL M casadi_tanh(const M& x) { return tanh(x); }
DECL M casadi_sinh(const M& x) { return sinh(x); }
DECL M casadi_cosh(const M& x) { return cosh(x); }
DECL M casadi_atanh(const M& x) { return atanh(x); }
DECL M casadi_asinh(const M& x) { return asinh(x); }
DECL M casadi_acosh(const M& x) { return acosh(x); }
DECL M casadi_exp(const M& x) { return exp(x); }
DECL M casadi_log(const M& x) { return log(x); }
DECL M casadi_log10(const M& x) { return log10(x); }
DECL M casadi_log1p(const M& x) { return log1p(x); }
DECL M casadi_expm1(const M& x) { return expm1(x); }
DECL M casadi_floor(const M& x) { return floor(x); }
DECL M casadi_ceil(const M& x) { return ceil(x); }
DECL M casadi_erf(const M& x) { return erf(x); }
DECL M casadi_erfinv(const M& x) { using casadi::erfinv; return erfinv(x); }
DECL M casadi_sign(const M& x) { using casadi::sign; return sign(x); }
DECL M casadi_power(const M& x, const M& n) { return pow(x, n); }
DECL M casadi_mod(const M& x, const M& y) { return fmod(x, y); }
DECL M casadi_fmod(const M& x, const M& y) { return fmod(x, y); }
DECL M casadi_remainder(const M& x, const M& y) { return remainder(x, y); }
DECL M casadi_atan2(const M& x, const M& y) { return atan2(x, y); }
DECL M casadi_fmin(const M& x, const M& y) { return fmin(x, y); }
DECL M casadi_fmax(const M& x, const M& y) { return fmax(x, y); }
DECL M casadi_hypot(const M& x, const M& y) { return hypot(x, y); }
DECL M casadi_simplify(const M& x) { using casadi::simplify; return simplify(x); }
DECL bool casadi_is_equal(const M& x, const M& y, casadi_int depth=0) { using casadi::is_equal; return is_equal(x, y, depth); }
DECL M casadi_copysign(const M& x, const M& y) { return copysign(x, y); }
DECL M casadi_constpow(const M& x, const M& y) { using casadi::constpow; return constpow(x, y); }
#endif // FLAG & IS_MEMBER
%enddef

%define GENERIC_EXPRESSION_ALL(DECL, FLAG)
GENERIC_EXPRESSION_FUN(DECL, (FLAG | IS_MX), MX)
GENERIC_EXPRESSION_FUN(DECL, (FLAG | IS_DMATRIX), Matrix<double>)
GENERIC_EXPRESSION_FUN(DECL, (FLAG | IS_SX), Matrix<SXElem>)
GENERIC_EXPRESSION_FUN(DECL, (FLAG | IS_DOUBLE), double)
%enddef

%define MATRIX_FUN(DECL, FLAG, M)
#if FLAG & IS_MEMBER
DECL M casadi_all(const M& x) { return all(x); }
DECL M casadi_any(const M& x) { return any(x); }
DECL M casadi_adj(const M& A) { return adj(A); }
DECL M casadi_minor(const M& x, casadi_int i, casadi_int j) { return minor(x, i, j); }
DECL M casadi_cofactor(const M& x, casadi_int i, casadi_int j) { return cofactor(x, i, j); }
DECL void casadi_qr(const M& A, M& OUTPUT1, M& OUTPUT2) { return qr(A, OUTPUT1, OUTPUT2); }
DECL void casadi_qr_sparse(const M& A, M& OUTPUT1, M& OUTPUT2, M& OUTPUT3, std::vector<casadi_int>& OUTPUT4, std::vector<casadi_int>& OUTPUT5, bool amd=true) { return qr_sparse(A, OUTPUT1, OUTPUT2, OUTPUT3, OUTPUT4, OUTPUT5, amd); }
DECL M casadi_qr_solve(const M& b, const M& v, const M& r, const M& beta, const std::vector<casadi_int>& prinv, const std::vector<casadi_int>& pc, bool tr=false) { return qr_solve(b, v, r, beta, prinv, pc, tr); }
DECL void casadi_ldl(const M& A, M& OUTPUT1, M& OUTPUT2, std::vector<casadi_int>& OUTPUT3, bool amd=true) { return ldl(A, OUTPUT1, OUTPUT2, OUTPUT3, amd); }
DECL M casadi_ldl_solve(const M& b, const M& D, const M& LT, const std::vector<casadi_int>& p) { return ldl_solve(b, D, LT, p); }
DECL M casadi_chol(const M& A) { return chol(A); }
DECL M casadi_norm_inf_mul(const M& x, const M& y) { return norm_inf_mul(x, y); }
DECL M casadi_sparsify(const M& A, double tol=0) { return sparsify(A, tol); }
DECL void casadi_expand(const M& ex, M& OUTPUT1, M& OUTPUT2) { expand(ex, OUTPUT1, OUTPUT2); }
DECL M casadi_pw_const(const M &t, const M& tval, const M& val) { return pw_const(t, tval, val); }
DECL M casadi_pw_lin(const M& t, const M& tval, const M& val) { return pw_lin(t, tval, val); }
DECL M casadi_heaviside(const M& x) { return heaviside(x); }
DECL M casadi_rectangle(const M& x) { return rectangle(x); }
DECL M casadi_triangle(const M& x) { return triangle(x); }
DECL M casadi_ramp(const M& x) { return ramp(x); }
DECL M casadi_gauss_quadrature(const M& f, const M& x, const M& a, const M& b, casadi_int order=5) { return gauss_quadrature(f, x, a, b, order); }
DECL M casadi_gauss_quadrature(const M& f, const M& x, const M& a, const M& b, casadi_int order, const M& w) { return gauss_quadrature(f, x, a, b, order, w); }
DECL M casadi_taylor(const M& ex, const M& x, const M& a=0, casadi_int order=1) { return taylor(ex, x, a, order); }
DECL M casadi_mtaylor(const M& ex, const M& x, const M& a, casadi_int order=1) { return mtaylor(ex, x, a, order); }
DECL M casadi_mtaylor(const M& ex, const M& x, const M& a, casadi_int order, const std::vector<casadi_int>& order_contributions) { return mtaylor(ex, x, a, order, order_contributions); }
DECL M casadi_poly_coeff(const M& ex, const M&x) { return poly_coeff(ex, x); }
DECL M casadi_poly_roots(const M& p) { return poly_roots(p); }
DECL M casadi_eig_symbolic(const M& m) { return eig_symbolic(m); }

#endif
%enddef

%define MATRIX_ALL(DECL, FLAG)
MATRIX_FUN(DECL, (FLAG | IS_DMATRIX), Matrix<double>)
MATRIX_FUN(DECL, (FLAG | IS_SX), Matrix<SXElem>)
%enddef

%define MX_FUN(DECL, FLAG, M)
#if FLAG & IS_MEMBER
DECL M casadi_find(const M& x) { return find(x); }
DECL M casadi_low(const M& v, const M& p, const Dict& options = Dict()) { return low(v, p, options); }
DECL M casadi_inv_node(const M& x) { return inv_node(x); }
#endif // FLAG & IS_MEMBER

#if FLAG & IS_GLOBAL
DECL std::vector< M > casadi_matrix_expand(const std::vector< M >& e, const std::vector< M > &boundary = std::vector< M >(), const Dict& options = Dict()) { return matrix_expand(e, boundary, options); }
DECL M casadi_matrix_expand(const M& e, const std::vector< M > &boundary = std::vector< M >(), const Dict& options = Dict()) { return matrix_expand(e, boundary, options); }
DECL M casadi_graph_substitute(const M& ex, const std::vector< M >& v, const std::vector< M > &vdef) { return graph_substitute(ex, v, vdef); }
DECL std::vector< M > casadi_graph_substitute(const std::vector< M > &ex, const std::vector< M > &v, const std::vector< M > &vdef) { return graph_substitute(ex, v, vdef); }
DECL M casadi_bspline(const M& x, const DM& coeffs, const std::vector< std::vector<double> >& knots, const std::vector<casadi_int>& degree, casadi_int m, const Dict& opts = Dict()) { return bspline(x, coeffs, knots, degree, m, opts); }
DECL M casadi_bspline(const M& x, const M& coeffs, const std::vector< std::vector<double> >& knots, const std::vector<casadi_int>& degree, casadi_int m, const Dict& opts = Dict()) { return bspline(x, coeffs, knots, degree, m, opts); }
DECL MX casadi_bspline(const MX& x, const MX& coeffs, const std::vector<MX>& knots, const std::vector<casadi_int>& degree, casadi_int m, const Dict& opts = Dict()) { return bspline(x, coeffs, knots, degree, m, opts); }
DECL M casadi_convexify(const M& H, const Dict& opts = Dict()) { return convexify(H, opts); }
DECL M casadi_stop_diff(const M& expr, casadi_int order) { return stop_diff(expr, order); }
DECL M casadi_stop_diff(const M& expr, const M& var, casadi_int order) { return stop_diff(expr, var, order); }
DECL std::vector< M > casadi_difference(const std::vector< M >& a, const std::vector< M >& b) { return difference(a, b); }
DECL M casadi_no_hess(const M& expr) { return no_hess(expr); }
DECL M casadi_no_grad(const M& expr) { return no_grad(expr); }
DECL M casadi_kron_contract(const M& m, const M& x, bool inner) { return kron_contract(m, x, inner); }

#endif
%enddef

%define MX_ALL(DECL, FLAG)
MX_FUN(DECL, (FLAG | IS_MX), MX)
%enddef
%include <casadi/core/matrix_fwd.hpp>
%include <casadi/core/matrix_decl.hpp>
%include <casadi/core/dm_fwd.hpp>
%include <casadi/core/sx_fwd.hpp>

// Remove from API
%warnfilter(401) casadi::Matrix<casadi_int>;
%template() casadi::Matrix<casadi_int>;

#ifdef SWIGJAVA
// Kotlin operators (plus, minus, times, div, rem) and named relatives, through
// the module functions, as Python's operators go through _casadi.<op>
%define %java_binop(Self, Op, F)
  public Self Op(Self y) { return ca.F(this, y); }
  public Self Op(double y) { return ca.F(this, y); }
%enddef
%define %java_binop_dm(Op, F)
  %java_binop(DM, Op, F)
  public SX Op(SX y) { return ca.F(this, y); }
  public MX Op(MX y) { return ca.F(this, y); }
%enddef
%define %java_binop_sym(Self, Op, F)
  %java_binop(Self, Op, F)
  public Self Op(DM y) { return ca.F(this, y); }
%enddef
%define %java_binops(M, Self)
  M(Self, plus, plus) M(Self, minus, minus) M(Self, times, times) M(Self, div, rdivide)
  M(Self, rem, mod) M(Self, mtimes, mtimes) M(Self, power, power)
  M(Self, lt, lt) M(Self, le, le) M(Self, gt, gt) M(Self, ge, ge) M(Self, eq, eq) M(Self, ne, ne)
%enddef
%define %java_binop_dm2(Self, Op, F) %java_binop_dm(Op, F) %enddef
%extend casadi::Matrix<double> {
  %proxycode {
  %java_binops(%java_binop_dm2, DM)
  }
}
%extend casadi::Matrix<casadi::SXElem> {
  %proxycode {
  %java_binops(%java_binop_sym, SX)
  }
}
%extend casadi::MX {
  %proxycode {
  %java_binops(%java_binop_sym, MX)
  }
}
#endif // SWIGJAVA

%template(DM) casadi::Matrix<double>;
%extend casadi::Matrix<double> {
   %template(DM) Matrix<SXElem>;
};


namespace casadi{
  %extend Matrix<double> {
#ifndef SWIGWASMJS
    void assign(const casadi::Matrix<double>&rhs) { (*$self)=rhs; }
#endif
    %matrix_helpers(casadi::Matrix<double>, _MIndex)

  }

}

#ifdef SWIGPYTHON
  %feature("nothread") casadi::Matrix<double>::full;
  %feature("nothread") casadi::Matrix<double>::sparse;
#endif

#ifdef SWIGJAVA
namespace casadi{
  %extend Matrix<double> {
    // double[] as a column, double[][] one array per row
    Matrix(const std::vector<double>& x) { return new casadi::Matrix<double>(x); }
    Matrix(const std::vector< std::vector<double> >& x) { return new casadi::Matrix<double>(x); }
    // Dense entries, one array per row
    std::vector< std::vector<double> > full() const {
      std::vector<double> d($self->numel());
      casadi::casadi_densify($self->ptr(), $self->sparsity(), casadi::get_ptr(d), false);
      std::vector< std::vector<double> > r($self->size1(), std::vector<double>($self->size2()));
      for (casadi_int i=0; i<$self->size1(); ++i)
        for (casadi_int j=0; j<$self->size2(); ++j) r[i][j] = d[i + j*$self->size1()];
      return r;
    }
  }
} // namespace casadi
#endif // SWIGJAVA

// Extend DM with SWIG unique features
#if !defined(SWIGWASMJS) && !defined(SWIGJAVA)
namespace casadi{
  %extend Matrix<double> {
    // Convert to a dense matrix
    GUESTOBJECT* full() const {
      return full(*$self);
    }

    // Convert to a sparse matrix
    GUESTOBJECT* sparse() const {
      return sparse(*$self);
    }
  }

} // namespace casadi
#endif


#ifdef SWIGPYTHON
namespace casadi{
%extend Matrix<double> {

%python_array_wrappers(999.0, np.float64)

%pythoncode %{
  def tocsc(self):
    import numpy as np
    import warnings
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      from scipy.sparse import csc_matrix
    return csc_matrix( (self.nonzeros(),self.row(),self.colind()), shape = self.shape, dtype=np.double )
  def tocsc_array(self):
    import numpy as np
    import warnings
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      from scipy.sparse import csc_array
    return csc_array( (self.nonzeros(),self.row(),self.colind()), shape = self.shape, dtype=np.double )
  def toarray(self,simplify=False):
    import numpy as np
    if simplify:
      if self.is_scalar():
        return float(self)
      elif self.is_vector():
        return np.array(self.elements())
    # elements() is column-major; avoids materializing the transpose
    return np.ascontiguousarray(np.array(self.elements()).reshape(self.shape, order='F'))
  def to_numpy(self):
    # Pandas-style densification hook.  matplotlib's
    # cbook._unpack_to_numpy looks for `to_numpy()` before falling back
    # to iteration, so defining it here makes plt.plot(DM) work without
    # forcing densification through NEP-18 shape ops like np.atleast_1d.
    return self.full()
%}
// scipy.sparse result, left untyped
%stub_method0(tocsc,    Any)
%stub_method0(tocsc_array, Any)
%stub_method(toarray,   %arg(NDArray[np.float64] | float), simplify: bool = ...)
%stub_method0(to_numpy, %arg(NDArray[np.float64]))


%pythoncode %{
  def __bool__(self):
    if self.numel()!=1:
      raise Exception("Only a scalar can be cast to a float")
    if self.nnz()==0:
      return False
    return float(self)!=0
%}

}; // extend Matrix<double>


// Logic for pickling

%extend Matrix<double> {

  %pythoncode %{
    def __setstate__(self, state):
        self.__init__(DM.deserialize(state["serialization"]))

    def __getstate__(self):
        return {"serialization": self.serialize()}
  %}
  %stub_method0(__getstate__, %arg(dict[builtins.str, Any]))
  %stub_method(__setstate__, None, state: Mapping[builtins.str, Any])

}


%extend Function {

  %pythoncode %{
    def __setstate__(self, state):
        self.__init__(Function.deserialize(state["serialization"]))

    def __getstate__(self):
        return {"serialization": self.serialize()}
  %}
  %stub_method0(__getstate__, %arg(dict[builtins.str, Any]))
  %stub_method(__setstate__, None, state: Mapping[builtins.str, Any])

}


} // namespace casadi
#endif // SWIGPYTHON


#ifdef SWIGMATLAB
namespace casadi{


%extend Matrix<double> {

  %matlabcode %{
     function s = saveobj(obj)
        try
            s.serialization = obj.serialize();
        catch exception
            warning(['Serializing of CasADi DM failed:' getReport(exception) ]);
            s = struct;
        end
     end
  %}
  %matlabcode_static %{
     function obj = loadobj(s)
        try
          if isstruct(s)
             obj = casadi.DM.deserialize(s.serialization);
          else
             obj = s;
          end
        catch exception
            warning(['Serializing of CasADi DM failed:' getReport(exception) ]);
            s = struct;
        end
     end
  %}
}

%extend Sparsity {
  %matlabcode %{
     function s = saveobj(obj)
        try
            s.serialization = obj.serialize();
        catch exception
            warning(['Serializing of CasADi Sparsity failed:' getReport(exception) ]);
            s = struct;
        end
     end
  %}
  %matlabcode_static %{
     function obj = loadobj(s)
        try
          if isstruct(s)
             obj = casadi.Sparsity.deserialize(s.serialization);
          else
             obj = s;
          end
        catch exception
            warning(['Serializing of CasADi Sparsity failed:' getReport(exception) ]);
            s = struct;
        end
     end
  %}

  // Transpose using the A' syntax in addition to A.'
  Sparsity ctranspose() const { return $self->T();}
}


%extend Function {

  %matlabcode %{
     function s = saveobj(obj)
        try
            s.serialization = obj.serialize();
        catch exception
            warning(['Serializing of CasADi Function failed:' getReport(exception) ]);
            s = struct;
        end
     end
  %}
  %matlabcode_static %{
     function obj = loadobj(s)
        try
          if isstruct(s)
             obj = casadi.Function.deserialize(s.serialization);
          else
             obj = s;
          end
        catch exception
            warning(['Serializing of CasADi Function failed:' getReport(exception) ]);
            s = struct;
        end
     end
  %}

}

} // namespace casadi
#endif // SWIGMATLAB

%include <casadi/core/sx_elem.hpp>

#ifdef SWIGPYTHON
%extend casadi::Sparsity{
    %pythoncode %{
        @property
        def shape(self):
            return (self.size1(),self.size2())

        @property
        def T(self):
            return _casadi.transpose(self)

        def __array__(self,*args,**kwargs):
            return DM.ones(self).full()
    %}
};

%pythoncode %{

try:
  import numpy

  def constpow(x,y):
    pass

  constpow=numpy.frompyfunc(constpow,2,1)
except:
  pass
%}
#endif // SWIGPYTHON

namespace casadi {
%extend Matrix<SXElem>{
    %matrix_helpers(casadi::Matrix<casadi::SXElem>, _MIndex)

  #ifdef SWIGPYTHON
  %python_array_wrappers(1001.0, Any)
  #endif // SWIGPYTHON

};

} // namespace casadi

#ifdef SWIGPYTHON
#include <arrayobject.h>
%template()    std::vector<PyObject*>;
#endif // SWIGPYTHON

%template(SX) casadi::Matrix<casadi::SXElem>;
%extend casadi::Matrix<casadi::SXElem> {
   %template(SX) Matrix<double>;
};

%include <casadi/core/mx.hpp>

%extend casadi::MX{
  %matrix_helpers(casadi::MX, %arg(_MIndex | MX))
  #ifdef SWIGPYTHON
  %python_array_wrappers(1002.0, Any)
  #endif //SWIGPYTHON
};

#ifdef SWIGPYTHON
%pythoncode %{
def attach_return_type(f,t):
  if not(hasattr(f,'func_annotations')):
    f.func_annotations = {}
  if not(isinstance(getattr(f,'func_annotations'),dict)):
    raise Exception("Cannot annotate this python Method to be a sparsitygenerator. Method has func_annotations attribute with unknown type.")
  f.func_annotations["return"] = t
  return f

def pyevaluate(f):
  return attach_return_type(f,None)

def pycallback(f):
  return attach_return_type(f,int)


def pyfunction(inputs,outputs):
  def wrap(f):

    @pyevaluate
    def fcustom(f2):
      res = f([f2.getInput(i) for i in range(f2.n_in())])
      if not isinstance(res,list):
        res = [res]
      for i in range(f2.n_out()):
        f2.setOutput(res[i],i)
    import warnings

    with warnings.catch_warnings():
      warnings.filterwarnings("ignore",category=DeprecationWarning)
      Fun = CustomFunction("CustomFunction",fcustom,inputs,outputs)
      return Fun

  return wrap

def PyFunction(name, obj, inputs, outputs, opts={}):
    @pyevaluate
    def fcustom(f):
      res = [f.getOutput(i) for i in range(f.n_out())]
      obj.evaluate([f.getInput(i) for i in range(f.n_in())],res)
      for i in range(f.n_out()): f.setOutput(res[i], i)

    import warnings

    with warnings.catch_warnings():
      warnings.filterwarnings("ignore",category=DeprecationWarning)
      return CustomFunction("CustomFunction", fcustom,
                            inputs, outputs, opts)

%}
#endif

#ifndef SWIGPYTHON
%ignore FunctionBuffer;
%ignore _function_buffer_eval;
#endif

%include <casadi/core/function.hpp>
#ifdef SWIGPYTHON
namespace casadi{
%extend Function {
  %pythoncode %{
    def __call__(self, *args, **kwargs):
      # Either named inputs or ordered inputs
      if len(args)>0 and len(kwargs)>0:
        raise SyntaxError('Function evaluation requires all arguments to be named or none')
      if len(args)>0:
        # Ordered inputs -> return tuple
        ret = self.call(args)
        if len(ret)==0:
          return None
        elif len(ret)==1:
          return ret[0]
        else:
          return tuple(ret)
      else:
        # Named inputs -> return dictionary
        return self.call(kwargs)

    def buffer(self):
      """
      Create a FunctionBuffer object for evaluating with minimal overhead

      """
      import functools
      fb = FunctionBuffer(self)
      caller = functools.partial(_casadi._function_buffer_eval, fb._self())
      return (fb, caller)
  %}
  // buffer() returns (FunctionBuffer, callable)
  %stub_method0(buffer, %arg(tuple[FunctionBuffer, Any]))
  // __call__ narrows per homogeneous input type; f() and f(**kw) give dicts
  %stub_overload_method(__call__, DM, __arg0: _DM, *args: _DM)
  %stub_overload_method(__call__, SX, __arg0: _SX, *args: _SX)
  %stub_overload_method(__call__, MX, __arg0: _MX, *args: _MX)
  %stub_overload_method(__call__, %arg(dict[builtins.str, DM]), **kwargs: _DM)
  %stub_overload_method(__call__, %arg(dict[builtins.str, SX]), **kwargs: _SX)
  %stub_overload_method(__call__, %arg(dict[builtins.str, MX]), **kwargs: _MX)
  %stub_overload_method(__call__, %arg(dict[builtins.str, DM | SX | MX]), **kwargs: float | DM | SX | MX)
 }

}
#endif // SWIGPYTHON

#ifdef SWIGMATLAB
namespace casadi{
%extend GenericMatrixCommon {
  %matlabcode %{
    function varargout = spy(self,varargin)
      [varargout{1:nargout}] = spy(sparse(casadi.DM(self.sparsity(),1)),varargin{:});
    end
    function varargout = subsref(self,s)
      if strcmp(s(1).type,'()')
        ret = paren(self, s(1).subs{:});
      elseif strcmp(s(1).type,'{}')
        ret = brace(self, s(1).subs{:});
      else
        [varargout{1:nargout}] = builtin('subsref',self,s);
        return
      end
      % Remaining subscripts act on the result, e.g. x(1,:).dim()
      if numel(s)==1
        varargout{1} = ret;
      else
        [varargout{1:nargout}] = subsref(ret, s(2:end));
      end
    end
    function n = numArgumentsFromSubscript(self,s,context)
      % One output for e.g. disp(x.sparsity); statements keep the default
      if ~exist('OCTAVE_VERSION', 'builtin') && context==matlab.mixin.util.IndexingContext.Statement
        n = builtin('numArgumentsFromSubscript',self,s,context);
      else
        n = 1;
      end
    end
    function self = subsasgn(self,s,v)
      if numel(s)==1 && strcmp(s.type,'()')
        paren_asgn(self, v, s.subs{:});
      elseif numel(s)==1 && strcmp(s.type,'{}')
        brace_asgn(self, v, s.subs{:});
      else
        self = builtin('subsasgn',self,s,v);
      end
    end
    function out = norm(self,varargin)
      narginchk(1,2);
      % 2-norm by default
      if nargin==1
        ind = 2;
      else
        ind = varargin{1};
      end
      % Typecheck
      assert((isnumeric(ind) && isscalar(ind)) || ischar(ind))
      % Pick the right norm
      if isnumeric(ind)
        switch ind
          case 1
            out = norm_1(self);
          case 2
            out = norm_2(self);
          case inf
            out = norm_inf(self);
          otherwise
            error(sprintf('Unknown norm argument: %g', ind))
        end
      else
        switch ind
          case 'fro'
            out = norm_fro(self);
          case 'inf'
            out = norm_inf(self);
          otherwise
            error(sprintf('Unknown norm argument: ''%s''', ind))
        end
      end
    end
    function out = min(varargin)
      narginchk(1,2);
      if nargin==1
        out = mmin(varargin{1});
      else
        out = fmin(varargin{1}, varargin{2});
      end
    end
    function out = max(varargin)
      narginchk(1,2);
      if nargin==1
        out = mmax(varargin{1});
      else
        out = fmax(varargin{1}, varargin{2});
      end
    end
    function b = isrow(self)
      b = is_row(self);
    end
    function b = iscolumn(self)
      b = is_column(self);
    end
    function b = isvector(self)
      b = is_vector(self);
    end
    function b = isscalar(self)
      b = is_scalar(self);
    end
  %}
}
%extend DaeBuilder {
  %matlabcode %{
    function varargout = subsref(self,s)
      if strcmp(s(1).type,'()')
        ret = paren(self, s(1).subs{:});
        if numel(s)==1
          varargout{1} = ret;
        else
          [varargout{1:nargout}] = subsref(ret, s(2:end));
        end
      else
        [varargout{1:nargout}] = builtin('subsref',self,s);
      end
    end
    function n = numArgumentsFromSubscript(self,s,context)
      % One output for e.g. fprintf('%g', dae.start_time); statements keep the default
      if ~exist('OCTAVE_VERSION', 'builtin') && context==matlab.mixin.util.IndexingContext.Statement
        n = builtin('numArgumentsFromSubscript',self,s,context);
      else
        n = 1;
      end
    end
  %}
}
%extend Function {
  %matlabcode %{
    function varargout = subsref(self,s)
      if numel(s)==1 && strcmp(s.type,'()')
        [varargout{1:nargout}]= paren(self, s.subs{:});
      else
        [varargout{1:nargout}] = builtin('subsref',self,s);
      end
   end
   function varargout = paren(self, varargin)
      if nargin==1 || (nargin>=2 && ischar(varargin{1}))
        % Named inputs: return struct
        assert(nargout<2, 'Syntax error');
        assert(mod(nargin,2)==1, 'Syntax error');
        arg = struct;
        for i=1:2:nargin-1
          assert(ischar(varargin{i}), 'Syntax error');
          arg.(varargin{i}) = varargin{i+1};
        end
        res = self.call(arg);
        varargout{1} = res;
      else
        % Ordered inputs: return variable number of outputs
        res = self.call(varargin);
        assert(nargout<=numel(res), 'Too many outputs');
        for i=1:max(min(1,numel(res)),nargout)
          varargout{i} = res{i};
        end
      end
    end
  %}
 }
}
#endif // SWIGMATLAB
%include <casadi/core/external.hpp>
%include <casadi/core/integrator.hpp>
%include <casadi/core/conic.hpp>
%include <casadi/core/nlpsol.hpp>
%include <casadi/core/onnx_function.hpp>
%include <casadi/core/graph_model.hpp>
%include <casadi/core/graph_builder.hpp>
%include <casadi/core/rootfinder.hpp>
%include <casadi/core/linsol.hpp>
%include <casadi/core/dple.hpp>
%include <casadi/core/expm.hpp>
%include <casadi/core/interpolant.hpp>
%include <casadi/core/blazing_spline.hpp>
%include <casadi/core/modelica_parser.hpp>

%feature("copyctor", "0") casadi::CodeGenerator;
%include <casadi/core/code_generator.hpp>

#ifdef SWIGMATLAB
// Wrap (static) member functions
%feature("nonstatic");
namespace casadi {
  %extend SparsityInterfaceCommon {
    SPARSITY_INTERFACE_ALL(static inline, IS_MEMBER)
  }
  %extend GenericExpressionCommon {
    GENERIC_EXPRESSION_ALL(static inline, IS_MEMBER)
  }
  %extend GenericMatrixCommon {
    GENERIC_MATRIX_ALL(static inline, IS_MEMBER)
  }
  %extend MatrixCommon {
    MATRIX_ALL(static inline, IS_MEMBER)
  }
  %extend MX {
    MX_ALL(static inline, IS_MEMBER)
    const MX brace(const casadi::MX& rr) const { casadi::MX m; $self->get_nz(m, true, rr); return m;}
    void brace_asgn(const MX& m, const casadi::MX& rr) { $self->set_nz(m, true, rr); }
    const MX paren(const casadi::MX& rr) const {
      casadi::MX m;
      $self->get(m, true, rr);
      return m;
    }
    %_paren_combos(%_paren2, casadi::MX, const casadi::MX&)
    // Needed for brace syntax to access nonzeros
    casadi_int numel(const MX &k) const {
      return 1;
    }
  }
} // namespace casadi
%feature("nonstatic", "");
// Member functions already wrapped
#define FLAG IS_GLOBAL
#else // SWIGMATLAB
// Need to wrap member functions below
#define FLAG (IS_GLOBAL | IS_MEMBER)
#endif // SWIGMATLAB

// Wrap non-member functions, possibly with casadi_ prefix

%inline {
  namespace casadi {
    SPARSITY_INTERFACE_ALL(inline, FLAG)
    GENERIC_EXPRESSION_ALL(inline, FLAG)
    GENERIC_MATRIX_ALL(inline, FLAG)
    MATRIX_ALL(inline, FLAG)
    MX_ALL(inline, FLAG)
  }
}

%inline {
  namespace casadi {
    inline DM casadi_bspline_dual(const std::vector<double>& x,
        const std::vector< std::vector<double> >& knots,
        const std::vector<casadi_int>& degree, const Dict& opts=Dict()) {
      return MX::bspline_dual(x, knots, degree, opts);
    }
  }
}

// Wrap the casadi_ prefixed functions in member functions
#ifdef SWIGPYTHON
%pythoncode %{
def _array_priority_yield(x, y):
    # numpy-style operator deferral: a binary op on a casadi value yields
    # (returns NotImplemented) to a NON-casadi operand that advertises a
    # higher __array_priority__ -- e.g. casadi.ArrayInterface -- so its
    # reflected operator runs and the higher-priority type controls the
    # result.  casadi's own DM/SX/MX promotion is untouched: those are
    # excluded, and their mixing is handled natively by the _casadi.* call.
    py = getattr(y, "__array_priority__", None)
    if py is None or py <= getattr(x, "__array_priority__", 0.0):
        return False
    return not isinstance(y, (DM, SX, MX, Sparsity))
%}
namespace casadi {
  %extend GenericExpressionCommon {
// Upstream SWIG raises TypeError on a failed overload dispatch, the fork NotImplementedError
#ifdef SWIG_STUBS_ENABLED
    %pythoncode %{
      def __hash__(self):
        try:
          return self.element_hash()
        except:
          return SharedObject.__hash__(self)
      def __matmul__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        try:
          return _casadi.mtimes(x, y)
        except NotImplementedError:
          return NotImplemented
      def __rmatmul__(x, y):
        try:
          return _casadi.mtimes(y, x)
        except NotImplementedError:
          return NotImplemented
      def __imatmul__(x, y):
        try:
          return _casadi.mtimes(x, y)
        except NotImplementedError:
          return NotImplemented
    %}
#else
    %pythoncode %{
      def __hash__(self):
        try:
          return self.element_hash()
        except:
          return SharedObject.__hash__(self)
      def __matmul__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        try:
          return _casadi.mtimes(x, y)
        except (NotImplementedError, TypeError):
          return NotImplemented
      def __rmatmul__(x, y):
        try:
          return _casadi.mtimes(y, x)
        except (NotImplementedError, TypeError):
          return NotImplemented
      def __imatmul__(x, y):
        try:
          return _casadi.mtimes(x, y)
        except (NotImplementedError, TypeError):
          return NotImplemented
    %}
#endif
    %stub_method0(__hash__, int)
    %formacro(%stub_CasadiMatrix_binop, __matmul__, __rmatmul__, __imatmul__)
  }
}
namespace casadi {
  %extend GenericExpressionCommon {
    %pythoncode %{
      def __add__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.plus(x, y)
      def __radd__(x, y): return _casadi.plus(y, x)
      def __sub__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.minus(x, y)
      def __rsub__(x, y): return _casadi.minus(y, x)
      def __mul__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.times(x, y)
      def __rmul__(x, y): return _casadi.times(y, x)
      def __truediv__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.rdivide(x, y)
      def __rtruediv__(x, y): return _casadi.rdivide(y, x)
      def __floordiv__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.floor(_casadi.rdivide(x, y))
      def __rfloordiv__(x, y): return _casadi.floor(_casadi.rdivide(y, x))
      def __mod__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return x - y * _casadi.floor(_casadi.rdivide(x, y))
      def __rmod__(x, y):
        return y - x * _casadi.floor(_casadi.rdivide(y, x))
      def __divmod__(x, y):
        q = _casadi.floor(_casadi.rdivide(x, y))
        return (q, x - y * q)
      def __rdivmod__(x, y):
        q = _casadi.floor(_casadi.rdivide(y, x))
        return (q, y - x * q)
      def __lt__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.lt(x, y)
      def __rlt__(x, y): return _casadi.lt(y, x)
      def __le__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.le(x, y)
      def __rle__(x, y): return _casadi.le(y, x)
      def __gt__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.lt(y, x)
      def __rgt__(x, y): return _casadi.lt(x, y)
      def __ge__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.le(y, x)
      def __rge__(x, y): return _casadi.le(x, y)
      def __eq__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        r = _casadi.eq(x, y)
        if r is NotImplemented and isinstance(x, SX) and isinstance(y, MX):
          raise Exception("Cannot compare SX and MX objects for equality")
        return r
      def __ne__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        r = _casadi.ne(x, y)
        if r is NotImplemented and isinstance(x, SX) and isinstance(y, MX):
          raise Exception("Cannot compare SX and MX objects for inequality")
        return r
      def __req__(x, y):
        r = _casadi.eq(y, x)
        if r is NotImplemented and isinstance(x, SX) and isinstance(y, MX):
          raise Exception("Cannot compare SX and MX objects for equality")
        return r
      def __rne__(x, y):
        r = _casadi.ne(y, x)
        if r is NotImplemented and isinstance(x, SX) and isinstance(y, MX):
          raise Exception("Cannot compare SX and MX objects for inequality")
        return r
      def __pow__(x, n, modulo=None):
        if modulo is None and _array_priority_yield(x, n): return NotImplemented
        p = _casadi.power(x, n)
        if modulo is None:
          return p
        return p - modulo * _casadi.floor(_casadi.rdivide(p, modulo))
      def __rpow__(n, x): return _casadi.power(x, n)
      def __round__(x, ndigits=None):
        if ndigits is None:
          return _casadi.if_else(x >= 0, _casadi.floor(x + 0.5),
                                 _casadi.ceil(x - 0.5))
        s = 10.0 ** ndigits
        y = x * s
        return _casadi.if_else(y >= 0, _casadi.floor(y + 0.5),
                               _casadi.ceil(y - 0.5)) / s
      def __trunc__(x):
        return _casadi.if_else(x >= 0, _casadi.floor(x), _casadi.ceil(x))
      def __floor__(x): return _casadi.floor(x)
      def __ceil__(x):  return _casadi.ceil(x)
      def __abs__(x):   return _casadi.fabs(x)
      def __iadd__(x, y):      return _casadi.plus(x, y)
      def __isub__(x, y):      return _casadi.minus(x, y)
      def __imul__(x, y):      return _casadi.times(x, y)
      def __itruediv__(x, y):  return _casadi.rdivide(x, y)
      def __ifloordiv__(x, y): return _casadi.floor(_casadi.rdivide(x, y))
      def __imod__(x, y):      return x - y * _casadi.floor(_casadi.rdivide(x, y))
      def __ipow__(x, n):      return _casadi.power(x, n)
      def __arctan2__(x, y): return _casadi.atan2(x, y)
      def __rarctan2__(y, x): return _casadi.atan2(x, y)
      def fmin(self, other): return _casadi.fmin(self, other)
      def fmax(self, other): return _casadi.fmax(self, other)
      def __fmin__(x, y): return _casadi.fmin(x, y)
      def __rfmin__(y, x): return _casadi.fmin(x, y)
      def __fmax__(x, y): return _casadi.fmax(x, y)
      def __rfmax__(y, x): return _casadi.fmax(x, y)
      def logic_and(self, other): return _casadi.logic_and(self, other)
      def logic_or(self, other): return _casadi.logic_or(self, other)
      def fabs(self): return _casadi.fabs(self)
      def sqrt(self): return _casadi.sqrt(self)
      def sin(self): return _casadi.sin(self)
      def cos(self): return _casadi.cos(self)
      def tan(self): return _casadi.tan(self)
      def arcsin(self): return _casadi.asin(self)
      def arccos(self): return _casadi.acos(self)
      def arctan(self): return _casadi.atan(self)
      def sinh(self): return _casadi.sinh(self)
      def cosh(self): return _casadi.cosh(self)
      def tanh(self): return _casadi.tanh(self)
      def arcsinh(self): return _casadi.asinh(self)
      def arccosh(self): return _casadi.acosh(self)
      def arctanh(self): return _casadi.atanh(self)
      def exp(self): return _casadi.exp(self)
      def log(self): return _casadi.log(self)
      def log10(self): return _casadi.log10(self)
      def log1p(self): return _casadi.log1p(self)
      def expm1(self): return _casadi.expm1(self)
      def floor(self): return _casadi.floor(self)
      def ceil(self): return _casadi.ceil(self)
      def erf(self): return _casadi.erf(self)
      def sign(self): return _casadi.sign(self)
      def fmod(self, other): return _casadi.fmod(self, other)
      def hypot(self, other): return _casadi.hypot(self, other)
      def remainder(self, other): return _casadi.remainder(self, other)
      def __copysign__(x, y): return _casadi.copysign(x, y)
      def __rcopysign__(y, x): return _casadi.copysign(x, y)
      def copysign(self, other): return _casadi.copysign(self, other)
      def rcopysign(self, other): return _casadi.copysign(other, self)
      def __constpow__(x, y): return _casadi.constpow(x, y)
      def __rconstpow__(y, x): return _casadi.constpow(x, y)
      def constpow(self, other): return _casadi.constpow(self, other)
      def rconstpow(self, other): return _casadi.constpow(other, self)
    %}
    /* PEP-484 stubs for the %pythoncode operators above. */
    %formacro(%stub_CasadiMatrix_binop, __add__, __radd__, __sub__, __rsub__, __mul__, __rmul__,
      __truediv__, __rtruediv__, __floordiv__, __rfloordiv__, __mod__, __rmod__)
    %stub_overload_method_selftyped(__pow__, DM, DM, other: _DM, modulo: _DM | None = None)
    %stub_overload_method_selftyped(__pow__, SX, SX, other: _SX, modulo: _SX | None = None)
    %stub_overload_method_selftyped(__pow__, MX, MX, other: _MX, modulo: _MX | None = None)
    %formacro(%stub_CasadiMatrix_binop, __rpow__, __iadd__, __isub__, __imul__, __itruediv__,
      __ifloordiv__, __imod__, __ipow__)
    %formacro(%stub_CasadiMatrix_unop, __abs__, __trunc__, __floor__, __ceil__)
    %formacro(%stub_CasadiMatrix_binop, __lt__, __le__, __gt__, __ge__, __eq__, __ne__, __req__,
      __rne__, __rlt__, __rle__, __rgt__, __rge__)
    %stub_method(__round__, Self, ndigits: int | None = None)
    %stub_overload_method_selftyped(__divmod__,  %arg(tuple[DM, DM]), DM, other: _DM)
    %stub_overload_method_selftyped(__divmod__,  %arg(tuple[SX, SX]), SX, other: _SX)
    %stub_overload_method_selftyped(__divmod__,  %arg(tuple[MX, MX]), MX, other: _MX)
    %stub_overload_method_selftyped(__rdivmod__, %arg(tuple[DM, DM]), DM, other: _DM)
    %stub_overload_method_selftyped(__rdivmod__, %arg(tuple[SX, SX]), SX, other: _SX)
    %stub_overload_method_selftyped(__rdivmod__, %arg(tuple[MX, MX]), MX, other: _MX)
    /* numpy-style method spellings of the elementwise functions. */
    %formacro(%stub_CasadiMatrix_unop, fabs, sqrt, sin, cos, tan, arcsin, arccos, arctan, sinh,
      cosh, tanh, arcsinh, arccosh, arctanh, exp, log, log10, log1p, expm1, floor, ceil, erf, sign)
    %formacro(%stub_CasadiMatrix_binop, fmin, fmax, logic_and, logic_or, fmod, hypot, remainder,
      copysign, rcopysign, constpow, rconstpow, __arctan2__, __rarctan2__, __fmin__, __rfmin__,
      __fmax__, __rfmax__, __copysign__, __rcopysign__, __constpow__, __rconstpow__)
  }

  %extend GenericMatrixCommon {
    %pythoncode %{
      def __mldivide__(x, y): return _casadi.mldivide(x, y)
      def __rmldivide__(y, x): return _casadi.mldivide(x, y)
      def __mrdivide__(x, y): return _casadi.mrdivide(x, y)
      def __rmrdivide__(y, x): return _casadi.mrdivide(x, y)
      def __mpower__(x, y): return _casadi.mpower(x, y)
      def __rmpower__(y, x): return _casadi.mpower(x, y)
    %}
    %formacro(%stub_CasadiMatrix_binop, __mldivide__, __rmldivide__, __mrdivide__, __rmrdivide__,
      __mpower__, __rmpower__)
  }

} // namespace casadi
#endif // SWIGPYTHON

%feature("director") casadi::Callback;

// Callback: hand-written stubs for the director-reshaped methods
#ifdef SWIGPYTHON
%feature("python:stub_skip") casadi::Callback::eval;
%feature("python:stub_skip") casadi::Callback::eval_buffer;
%feature("python:stub_skip") casadi::Callback::get_sparsity_in;
%feature("python:stub_skip") casadi::Callback::get_sparsity_out;
#endif

%include <casadi/core/importer.hpp>
#ifdef SWIGJULIA
// users subtype the abstract Callback, as Python subclasses Callback; the proxy is CallbackFunction
%rename(CallbackFunction) casadi::Callback;
%feature("julia:director", "Callback") casadi::Callback;
// a Callback cannot be copied (as in C++); CallbackFunction(obj) takes a Julia Callback
%ignore casadi::Callback::Callback(const Callback&);
// eval_buffer: arg and res reach Julia as vectors of Vector{Float64} views (nothing for a null)
%{
namespace casadi {
  jl_value_t* julia_buffers(double** p, const std::vector<casadi_int>& sz) {
    jl_value_t* vt = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
    jl_array_t* a = jl_alloc_array_1d(jl_apply_array_type((jl_value_t*)jl_any_type, 1), sz.size());
    JL_GC_PUSH2(&vt, &a);
    for (size_t i = 0; i < sz.size(); ++i) {
      jl_value_t* e = p[i] ? (jl_value_t*)jl_ptr_to_array_1d(vt, p[i], sz[i], 0) : jl_nothing;
      jl_array_ptr_set(a, i, e);
    }
    JL_GC_POP();
    return (jl_value_t*)a;
  }
} // namespace casadi
%}
%typemap(directorin, noblock=1) (const double** arg, const std::vector<casadi_int>& sizes_arg) {
  $input = casadi::julia_buffers(const_cast<double**>($1), $2);
}
%typemap(directorin, noblock=1) (double** res, const std::vector<casadi_int>& sizes_res) {
  $input = casadi::julia_buffers($1, $2);
}
#endif // SWIGJULIA
%include <casadi/core/callback.hpp>

#ifdef SWIGPYTHON
%extend casadi::Callback {
  // eval(arg, /): arg is a sequence of DM
  %stub_method(eval, %arg(Sequence[_DM]), %arg(__arg: Sequence[DM]))
  // eval_buffer(arg, res): tuples of memoryviews
  %stub_method(eval_buffer, int, %arg(__arg: tuple[memoryview, ...], res: tuple[memoryview, ...]))
  %stub_method(get_sparsity_in,  Sparsity, i: int)
  %stub_method(get_sparsity_out, Sparsity, i: int)
}
#endif
%include <casadi/core/global_options.hpp>

%include <casadi/core/casadi_meta.hpp>
#ifdef SWIGPYTHON
%extend casadi::CasadiMeta {
  static const char* swig_flags() { return CASADI_SWIG_FLAGS; }
};
#endif // SWIGPYTHON

%include <casadi/core/integration_tools.hpp>
%include <casadi/core/nlp_tools.hpp>
%include <casadi/core/tools.hpp>
%include <casadi/core/nlp_builder.hpp>
%include <casadi/core/dae_builder.hpp>
%include <casadi/core/xml_file.hpp>
%include <casadi/core/archiver.hpp>
%include <casadi/core/filesystem.hpp>
%include <casadi/core/blas.hpp>
%include <casadi/core/options.hpp>

%feature("copyctor", "0") casadi::SerializerBase;
%feature("copyctor", "0") casadi::DeserializerBase;
%feature("copyctor", "0") casadi::StringSerializer;
%feature("copyctor", "0") casadi::StringDeserializer;
%feature("copyctor", "0") casadi::FileSerializer;
%feature("copyctor", "0") casadi::FileDeserializer;
%nodefaultctor casadi::SerializerBase;
%nodefaultctor casadi::DeserializerBase;

#ifdef SWIGPYTHON
%rename("_pop_type") casadi::DeserializerBase::pop_type;
%rename("%(regex:/(SERIALIZED_\w+)/_\\1/)s", regextarget=1, fullname=1) "casadi::SerializerBase::SERIALIZED_\w+";
#endif // SWIG_PYTHON


#ifdef SWIGMATLAB
%rename("internal_pop_type") casadi::DeserializerBase::pop_type;
%rename("%(regex:/(SERIALIZED_\w+)/internal_\\1/)s", regextarget=1, fullname=1) "casadi::SerializerBase::SERIALIZED_\w+";
#endif // SWIG_PYTHON

%include <casadi/core/serializer.hpp>

#ifdef SWIGPYTHON
%extend casadi::DeserializerBase {
  %stub_method0(unpack, Any)
  %pythoncode %{
    def unpack(self):
      type = SerializerBase.type_to_string(self._pop_type())
      f = getattr(self, "blind_unpack_"+type)
      return f()
  %}
}
%pythoncode %{
veccat.__doc__, vertcat.__doc__, horzcat.__doc__, diagcat.__doc__ = \
    _veccat.__doc__, _vertcat.__doc__, _horzcat.__doc__, _diagcat.__doc__

try:
  import threading
  _thread_local = threading.local()
except:
  threading_available = True
  _thread_local = globals()

def _current_pickle_context():
  return getattr(_thread_local, "casadi_pickle_ctx", None)

def _current_unpickle_context():
  return getattr(_thread_local, "casadi_unpickle_ctx", None)

class global_pickle_context:
    def __enter__(self):
        self.ctx = StringSerializer()
        _thread_local.casadi_pickle_ctx = self.ctx
        return self.ctx

    def __exit__(self, *args):
        _thread_local.casadi_pickle_ctx = None
      
class global_unpickle_context:
    def __enter__(self):
        self.ctx = StringDeserializer("")
        _thread_local.casadi_unpickle_ctx = self.ctx
        return self.ctx

    def __exit__(self, *args):
        _thread_local.casadi_unpickle_ctx = None
%}

// Stubs for pure-pythoncode classes
%stub_class_begin(global_pickle_context)
%stub_method0(__enter__, StringSerializer)
%stub_method(__exit__, None, *args: object)

%stub_class_begin(global_unpickle_context)
%stub_method0(__enter__, StringDeserializer)
%stub_method(__exit__, None, *args: object)


#endif // SWIGPYTHON
#ifdef SWIGMATLAB
%extend casadi::DeserializerBase {
  %matlabcode %{
    function out = unpack(self)
      type = casadi.SerializerBase.type_to_string(self.internal_pop_type);
      out = self.(['blind_unpack_' type]);
    end
  %}
}
#endif // SWIGMATLAB
#ifdef SWIGWASMJS
%extend casadi::DeserializerBase {
  %jscode %{
    DeserializerBase.prototype.unpack = function () {
      const name = SerializerBase.type_to_string(this.pop_type());
      return this["blind_unpack_" + name]();
    };
  %}
}
#endif // SWIGWASMJS

%feature("director") casadi::OptiCallback;
#ifdef SWIGJULIA
// internal: opti.callback takes any Julia callable
%rename(_OptiCallback) casadi::OptiCallback;
%rename(_callback_class) casadi::Opti::callback_class;
%rename(_call) casadi::OptiCallback::call;
#endif

// Return-by-value
#ifdef SWIGJAVA
// A DM: a Java value has one static type
%typemap(jni) casadi::native_DM "jobject"
%typemap(jtype) casadi::native_DM "Object"
%typemap(jstype) casadi::native_DM "DM"
%typemap(javaout) casadi::native_DM { return $imclassname.swigCast($jnicall); }
%typemap(out, noblock=1, fragment="casadi_all") casadi::native_DM { $result = swig::from< casadi::DM >($1); }
#elif defined(SWIGWASMJS)
%typemap(out, doc="DM", tsstub_out="DM", noblock=1, fragment="casadi_all") casadi::native_DM {
  if(!($result = swig::from< casadi::DM >($1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type 'DM'.");
}
#elif defined(SWIGJULIA)
// Native Julia value: Float64, Vector or Matrix
%typemap(jlout) casadi::native_DM "(function (d::DM) n = nonzeros(densify(d)); r, c = size1(d), size2(d); r * c == 1 ? n[1] : (r == 1 || c == 1) ? n : Base.reshape(n, r, c) end)(DM(_check($call)))"
#else
%typemap(out, doc="double", pystub_out="float", noblock=1, fragment="casadi_all") casadi::native_DM {
  if(!($result = full_or_sparse($1, true))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type 'double'.");
}
#endif


%apply casadi_int &OUTPUT { Opti::ConstraintType &OUTPUT };

%typemap(argout, noblock=1,fragment="casadi_all") casadi::Opti::ConstraintType &OUTPUT {
  %append_output(swig::from(static_cast<casadi_int>(*$1)));
}

%typemap(in, doc="Opti.ConstraintType", pystub_in="Any", noblock=1, numinputs=0) casadi::Opti::ConstraintType &OUTPUT (casadi::Opti::ConstraintType m) {
 $1 = &m;
}


#if defined(SWIGPYTHON) || defined(SWIGJULIA)
%define make_property(class, name, type)
#ifdef SWIGPYTHON
  %rename(_ ## name) class ## :: ## name;
  %extend class {
    %pythoncode %{
      @property
      def name(self):
        return self._ ## name()
    %}
    %stub_attr(name, type)
  }
#else
  // Julia: obj.name is the value itself (opti.f, opti.x, sol.debug), not a bound method
  %feature("julia:property") class ## :: ## name;
#endif
%enddef


make_property(casadi::Opti, debug, OptiAdvanced);
make_property(casadi::Opti, advanced, OptiAdvanced);
make_property(casadi::OptiSol, opti, Opti);

%define make_property_opti(name, type)
  make_property(casadi::Opti, name, type);
%enddef

make_property_opti(f, MX)
make_property_opti(g, MX)
make_property_opti(x, MX)
make_property_opti(p, MX)
make_property_opti(lam_g, MX)
make_property_opti(lbg, MX)
make_property_opti(ubg, MX)
make_property_opti(nx, int)
make_property_opti(np, int)
make_property_opti(ng, int)
make_property_opti(x_linear_scale, DM)
make_property_opti(x_linear_scale_offset, DM)
make_property_opti(g_linear_scale, DM)
make_property_opti(f_linear_scale, float)

#ifdef SWIGPYTHON
make_property(casadi::Opti, casadi_solver, Function);
#else
// keep the casadi_ prefix: solver(...) is the setter
%rename(casadi_solver) casadi::OptiAdvanced::casadi_solver;
make_property(casadi::OptiAdvanced, casadi_solver, Function);
#endif
#endif // SWIGPYTHON || SWIGJULIA

// Record the caller's stack trace in the user dict of new variables, parameters and constraints
#ifdef SWIGPYTHON
  %rename(_variable) casadi::Opti::variable;
  %rename(_parameter) casadi::Opti::parameter;
  %rename(_subject_to) casadi::Opti::subject_to;
  %extend casadi::Opti {
    %pythoncode %{
      def parameter(self,*args):
        import sys
        import os
        try:
            frame = sys._getframe(1)
        except:
            frame = {}
        meta = {} if frame is None else {"stacktrace": [{"file":os.path.abspath(frame.f_code.co_filename),"line":frame.f_lineno,"name":frame.f_code.co_name}]}
        ret = self._parameter(*args)
        if len(meta)>0:
            self.update_user_dict(ret, meta)
        return ret

      def variable(self,*args):
        import sys
        import os
        try:
            frame = sys._getframe(1)
        except:
            frame = {}
        meta = {} if frame is None else {"stacktrace": [{"file":os.path.abspath(frame.f_code.co_filename),"line":frame.f_lineno,"name":frame.f_code.co_name}]}
        ret = self._variable(*args)
        if len(meta)>0:
            self.update_user_dict(ret, meta)
        return ret

      def subject_to(self,*args):
        if len(args)==0:
          return self._subject_to()
        import sys
        import os
        stacktrace = []
        for i in range(1,10000):
          try:
            frame = sys._getframe(i)
            stacktrace.append({"file":os.path.abspath(frame.f_code.co_filename),"line":frame.f_lineno,"name":frame.f_code.co_name})
          except Exception as e:
            break
        args = list(args)
        if len(args)==3 and isinstance(args[2],dict):
          args[2] = dict(args[2])
          if "stacktrace" not in args[2]:
            args[2]["stacktrace"] = stacktrace
        elif len(args)==2 and isinstance(args[1],dict):
          args[1] = dict(args[1])
          if "stacktrace" not in args[1]:
            args[1]["stacktrace"] = stacktrace
        elif len(args)==1:
          args = [args[0], {"stacktrace": stacktrace}]
        elif len(args)==2:
          args = [args[0], args[1], {"stacktrace": stacktrace}]
        ret = self._subject_to(*args)
        return ret
    %}
    /* Stubs for the three %pythoncode methods above. */
    %stub_overload_method(variable, MX, n: int = ..., m: int = ..., attribute: builtins.str = ...)
    %stub_overload_method(variable, MX, sp: Sparsity | DM, attribute: builtins.str = ...)
    %stub_overload_method(variable, MX, symbol: MX | SX | DM, attribute: builtins.str = ...)
    %stub_overload_method(parameter, MX, n: int = ..., m: int = ..., attribute: builtins.str = ...)
    %stub_overload_method(parameter, MX, sp: Sparsity | DM, attribute: builtins.str = ...)
    %stub_overload_method(parameter, MX, symbol: MX | SX | DM, attribute: builtins.str = ...)
    // subject_to(g[, linear_scale][, options])
    %stub_overload_method0(subject_to, None)
    %stub_overload_method(subject_to, None, g: MX | SX | DM | float)
    %stub_overload_method(subject_to, None, g: Sequence[MX | SX | DM | float])
    %stub_overload_method(subject_to, None, g: MX | SX | DM | float, options: Mapping[builtins.str, _GenericType])
    %stub_overload_method(subject_to, None, g: Sequence[MX | SX | DM | float], options: Mapping[builtins.str, _GenericType])
    %stub_overload_method(subject_to, None, g: MX | SX | DM | float, linear_scale: _DM)
    %stub_overload_method(subject_to, None, g: Sequence[MX | SX | DM | float], linear_scale: _DM)
    %stub_overload_method(subject_to, None, g: MX | SX | DM | float, linear_scale: _DM, options: Mapping[builtins.str, _GenericType])
    %stub_overload_method(subject_to, None, g: Sequence[MX | SX | DM | float], linear_scale: _DM, options: Mapping[builtins.str, _GenericType])
  }
#elif defined(SWIGMATLAB)
  %rename(internal_variable) casadi::Opti::variable;
  %rename(internal_parameter) casadi::Opti::parameter;
  %rename(internal_subject_to) casadi::Opti::subject_to;
  %extend casadi::Opti {
    %matlabcode %{
      function out = variable(self, varargin)
        st = dbstack('-completenames',1);
        if length(st)>0
          meta = struct('stacktrace', {num2cell(st)});
        else
          meta = struct;
        end
        out = self.internal_variable(varargin{:});
        self.update_user_dict(out, meta);
      end
      function out = parameter(self, varargin)
        st = dbstack('-completenames',1);
        if length(st)>0
          meta = struct('stacktrace', {num2cell(st)});
        else
          meta = struct;
        end
        out = self.internal_parameter(varargin{:});
        self.update_user_dict(out, meta);
      end
      function [] = subject_to(self, varargin)
        if length(varargin)==0
          self.internal_subject_to();
          return;
        end
        st = dbstack('-completenames',1);
        if length(st)>0
          meta = struct('stacktrace', {num2cell(st)});
        else
          meta = struct;
        end
        self.internal_subject_to(varargin{:});
        self.update_user_dict(varargin{1}, meta);
      end
    %}
  }
  // solver() without arguments is the deprecated MATLAB name of casadi_solver()
  %rename(casadi_solver) casadi::OptiAdvanced::casadi_solver;
  %extend casadi::OptiAdvanced {
    %matlabcode %{
      function varargout = solver(self, varargin)
        if isempty(varargin)
          warning('casadi:deprecated', ...
            'OptiAdvanced.solver() is deprecated, use casadi_solver() instead.');
          [varargout{1:nargout}] = self.casadi_solver();
        else
          [varargout{1:nargout}] = solver@casadi.Opti(self, varargin{:});
        end
      end
    %}
  }
#elif defined(SWIGJULIA)
  // variable/parameter/subject_to record the call site; the Julia layer wraps these
  %rename(_variable) casadi::Opti::variable;
  %rename(_parameter) casadi::Opti::parameter;
  %rename(_subject_to) casadi::Opti::subject_to;
  // the Opti proxy keeps the OptiCallback it registers, as Python's self._cb
  %feature("julia:keep") casadi::Opti::callback_class;
  // copy(::Opti) is an independent problem; the copy constructor shares it
  %feature("julia:copy", "copy") casadi::Opti;
#endif

%include <casadi/core/optistack.hpp>


#ifdef SWIGPYTHON
%extend casadi::Opti {
  %pythoncode %{

    @staticmethod
    def _callback(self,fh=None):
      if fh is None:
        self.callback_class();
        return
      class OptiCallbackHelper(OptiCallback):
          def __init__(self, callback):
            OptiCallback.__init__(self)
            self.callback = callback

          def call(self, i):
            self.callback(i)

      self._fh = fh
      self._cb = OptiCallbackHelper(fh);
      self.callback_class(self._cb);


    def callback(self,fh=None):
      self._callback(self,fh)


  %}
  // Opti.callback (%pythoncode); fh gets the iteration index
  %stub_method(callback, None, %arg(fh: Callable[[int], Any] | None = ...))

}
#endif

#ifdef SWIGMATLAB
%extend casadi::Opti {
  %matlabcode %{
    function [] = callback(self, varargin)
      casadi.OptiCallbackHelper.callback_setup(self, varargin{:})
    end
  %}
}

%{
#ifdef HAVE_OCTAVE
  // Mandatory as of Octave 10
  extern "C" const int __octave_mex_soversion__ = 1;
#endif
%}

#endif

%include <casadi/core/resource.hpp>

// Cleanup for dependent modules
%exception {
  $action
}
