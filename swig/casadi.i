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

/* Compat shim for SWIG without the -stubs patch.  The patched fork
 * defines %stub_* primitives (see Lib/python/pyuserdir.swg) and sets
 * SWIG_STUBS_ENABLED when -stubs is passed.  Vanilla SWIG has neither
 * the real macros nor the noop branch, so we declare everything
 * empty here.  Direct `%stubcode %{...%}` use sites must be guarded
 * with `#ifdef SWIG_STUBS_ENABLED` separately -- `%insert("stubs")`
 * isn't a known target in vanilla SWIG. */
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
#endif

  /// Data structure in the target language holding data
#ifdef SWIGPYTHON
#define GUESTOBJECT PyObject
#elif defined(SWIGMATLAB)
#define GUESTOBJECT mxArray
#elif defined(SWIGWASMJS)
#define GUESTOBJECT casadi_emval_struct
#else
#define GUESTOBJECT void
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

    std::string python_string_to_std_string(PyObject *str_py) {
#if SWIG_VERSION < 0x040200
      const char *str_char = SWIG_Python_str_AsChar(str_py);
      std::string str(str_char);
      SWIG_Python_str_DelForPy3(str_char);
#else
      PyObject *bytes = NULL;
      std::string str(SWIG_PyUnicode_AsUTF8AndSize(str_py, NULL, &bytes));
      Py_XDECREF(bytes);
#endif
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
      Py_DECREF(msg_py);
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

  // Octave-on-Windows seems to pick up superfluous +casadi
  // Make sure we exclude it
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

// Turn off the warnings that certain methods are effectively ignored, this seams to be a false warning,
// for example vertcat(SXVector), vertcat(DMVector) and vertcat(MXVector) appears to work fine
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

%begin %{
#define SWIG_PYTHON_OUTPUT_TUPLE
%}

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
  if m.__class__.__name__ == "csc_matrix":
    if len(m.shape)!=2: return False
    if check_only: return True
    return m.shape + (m.indptr.flat,m.indices.flat,m.data.flat)
  return False

%}
#endif // WITH_SWIGPYTHON


//  These are the following styles
// error
// overview
// group

%feature("autodoc", "1");

%feature("customdoc", "1");

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

// append works for strings

%naturalvar;

// Make data members read-only
%immutable;

// Make sure that a copy constructor is created
%copyctor;

#ifndef SWIGXML
%feature("compactdefaultargs","1");
//%feature("compactdefaultargs","0") casadi::taylor; // taylor function has a default argument for which the namespace is not recognised by SWIG
%feature("compactdefaultargs","0") casadi::Function::generateCode; // buggy
#endif //SWIGXML

// STL
#ifdef SWIGXML
namespace std {
  template<class T> class vector {};
  template<class A, class B> class pair {};
  template<class A, class B> class map {};
}
#else // SWIGXML
%include "stl.i"
#endif // SWIGXML

%include "doc.i"


// casadi requires SWIG >= 3.0 (Swig::DirectorException inherits from
// std::exception there, so the std::exception catch already covers it).
// Exceptions handling
%include "exception.i"
%exception {
  try {
    $action
  } catch(const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

// Python sometimes takes an approach to not check, but just try.
// It expects a python error to be thrown.
%exception __int__ {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

// See https://github.com/casadi/casadi/issues/701
// Recent numpys will only catch TypeError or ValueError in printing logic
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
#endif //SWIGPYTHON


#ifdef SWIGPYTHON

%{
#define SWIG_FILE_WITH_INIT
%}

#endif // SWIGPYTHON

%{
#define SWIG_Error_return(code, msg)  { std::cerr << "Error occured in CasADi SWIG interface code:" << std::endl << "  "<< msg << std::endl;SWIG_Error(code, msg); return 0; }
%}

#ifndef SWIGXML

// Can be overloaded by specifying before importing casadi.i
%fragment("casadi_extra_decl", "header") {}

%fragment("casadi_decl", "header",fragment="casadi_extra_decl") {
  namespace casadi {
    /* Check if Null or None */
    bool is_null(GUESTOBJECT *p);

    /* Typemaps from CasADi types to types in the interfaced language:
     *
     * to_ptr: Converts a pointer in interfaced language to C++:
     *   Input: GUESTOBJECT pointer p
     *   Output: Pointer to pointer: At input, pointer to pointer to temporary
     *   The routine will either:
     *     - Do nothing, if 0
     *     - Change the pointer
     *     - Change the temporary object
     *   Returns true upon success, else false
     *
     * from_ptr: Converts result from CasADi to interfaced language
     */

    // Basic types
    bool to_ptr(GUESTOBJECT *p, bool** m);
    GUESTOBJECT* from_ptr(const bool *a);
    bool to_ptr(GUESTOBJECT *p, casadi_int** m);
    GUESTOBJECT* from_ptr(const casadi_int *a);
    bool to_ptr(GUESTOBJECT *p, double** m);
    GUESTOBJECT* from_ptr(const double *a);
    bool to_ptr(GUESTOBJECT *p, std::string** m);
    GUESTOBJECT* from_ptr(const std::string *a);

    // std::vector
#ifdef SWIGMATLAB
    bool to_ptr(GUESTOBJECT *p, std::vector<double> **m);
    GUESTOBJECT* from_ptr(const std::vector<double> *a);
    bool to_ptr(GUESTOBJECT *p, std::vector<casadi_int>** m);
    GUESTOBJECT* from_ptr(const std::vector<casadi_int> *a);
    bool to_ptr(GUESTOBJECT *p, std::vector<bool> **m);
    GUESTOBJECT* from_ptr(const std::vector<bool> *a);
    bool to_ptr(GUESTOBJECT *p, std::vector<std::string>** m);
    GUESTOBJECT* from_ptr(const std::vector<std::string> *a);
#endif // SWIGMATLAB
    template<typename M> bool to_ptr(GUESTOBJECT *p, std::vector<M>** m);
    template<typename M> bool to_ptr(GUESTOBJECT *p, std::vector< std::vector<M> >** m);
    template<typename M> GUESTOBJECT* from_ptr(const std::vector<M> *a);

    // std::pair
#ifdef SWIGMATLAB
    bool to_ptr(GUESTOBJECT *p, std::pair<casadi_int, casadi_int>** m);
    GUESTOBJECT* from_ptr(const std::pair<casadi_int, casadi_int>* a);
#endif // SWIGMATLAB
    template<typename M1, typename M2> bool to_ptr(GUESTOBJECT *p, std::pair<M1, M2>** m);
    template<typename M1, typename M2> GUESTOBJECT* from_ptr(const std::pair<M1, M2>* a);

    // std::map
    template<typename M> bool to_ptr(GUESTOBJECT *p, std::map<std::string, M>** m);
    template<typename M> GUESTOBJECT* from_ptr(const std::map<std::string, M> *a);


    // Slice
    bool to_ptr(GUESTOBJECT *p, casadi::Slice** m);
    GUESTOBJECT* from_ptr(const casadi::Slice *a);

    // Sparsity
    bool to_ptr(GUESTOBJECT *p, casadi::Sparsity** m);
    GUESTOBJECT* from_ptr(const casadi::Sparsity *a);

    // Matrix<>
    bool to_ptr(GUESTOBJECT *p, casadi::DM** m);
    GUESTOBJECT* from_ptr(const casadi::DM *a);
    bool to_ptr(GUESTOBJECT *p, casadi::IM** m);
    GUESTOBJECT* from_ptr(const casadi::IM *a);
    bool to_ptr(GUESTOBJECT *p, casadi::SX** m);
    GUESTOBJECT* from_ptr(const casadi::SX *a);

    // MX
    bool to_ptr(GUESTOBJECT *p, casadi::MX** m);
    GUESTOBJECT* from_ptr(const casadi::MX *a);

    // Function
    bool to_ptr(GUESTOBJECT *p, casadi::Function** m);
    GUESTOBJECT* from_ptr(const casadi::Function *a);

    // SXElem
    bool to_ptr(GUESTOBJECT *p, casadi::SXElem** m);
    GUESTOBJECT* from_ptr(const casadi::SXElem *a);

    // GenericType
    bool to_ptr(GUESTOBJECT *p, casadi::GenericType** m);
    GUESTOBJECT* from_ptr(const casadi::GenericType *a);

    // Same as to_ptr, but with pointer instead of pointer to pointer
    template<typename M> bool to_val(GUESTOBJECT *p, M* m);

#ifdef SWIGWASMJS
    // Generic fallback for wasm-js: classes that aren't explicitly
    // registered (e.g. casadi::SharedObject base) still need to be
    // probe-able via can_convert<>.  Delegates to the runtime
    // SWIG_ConvertPtr (val-based, see Lib/wasm_js/wasm_jsrun.swg) --
    // accepts the input only if it's already an instance carrying a
    // `_ptr` property.  Explicit to_ptr overloads (e.g.
    // `to_ptr(GUESTOBJECT*, casadi::Function**)`) win via standard C++
    // overload resolution; this template is the catch-all.
    template<typename M> bool to_ptr(GUESTOBJECT *p, M** m) {
      if (is_null(p)) return false;
      void* raw = 0;
      if (!SWIG_IsOK(SWIG_ConvertPtr(p, &raw, 0, 0))) return false;
      if (m) *m = static_cast<M*>(raw);
      return true;
    }
#endif

    // Check if conversion is possible
    template<typename M> bool can_convert(GUESTOBJECT *p) { return to_ptr(p, static_cast<M**>(0));}

    // Same as the above, but with reference instead of pointer
    template<typename M> GUESTOBJECT* from_ref(const M& m) { return from_ptr(&m);}

    // Specialization for std::vectors of booleans
    GUESTOBJECT* from_ref(std::vector<bool>::const_reference m) {
      bool tmp = m;
      return from_ptr(&tmp);
    }

    // Same as the above, but with a temporary object
    template<typename M> GUESTOBJECT* from_tmp(M m) { return from_ptr(&m);}
#ifdef SWIGMATLAB
    // Get sparsity pattern
    Sparsity get_sparsity(const mxArray* p);

    // Number of nonzeros
    size_t getNNZ(const mxArray* p);
#endif // SWIGMATLAB

    GUESTOBJECT* full(const DM& m, bool simplify=false) {
#ifdef SWIGPYTHON
      PyObject *p = from_ptr(&m);
      PyObject *method_name = PyString_FromString("toarray");
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
      PyObject *p = from_ptr(&m);
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
        if (is_null(p)) return from_ptr(&m);
        return p;
      }
    }
#ifdef SWIGPYTHON


    PyObject* get_Python_helper(const std::string& name) {
%#if PY_VERSION_HEX < 0x03070000
      PyObject* module = PyImport_AddModule("casadi");
%#else
      PyObject* c_name = PyString_FromString("casadi");
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
      // None signals "not handled by this helper" (issue #4216):
      // without this check, conv() would dereference None and segfault.
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
      if (!to_val(PyTuple_GetItem(p, 2), &data)) return false;
      casadi_int nrow; to_val(PyTuple_GetItem(p, 0), &nrow);
      casadi_int ncol; to_val(PyTuple_GetItem(p, 1), &ncol);
      if (m) {
        **m = casadi::SX::zeros(nrow, ncol);
        casadi_densify(get_ptr(data), (**m).sparsity().T(), (**m).ptr(), true);
      }
      return true;
    }

    bool IM_from_array_conv(GUESTOBJECT *p, casadi::IM** m) {
      if (!m) return true;
      std::vector<casadi_int> data;
      if (!to_val(PyTuple_GetItem(p, 2), &data)) return false;
      casadi_int nrow; to_val(PyTuple_GetItem(p, 0), &nrow);
      casadi_int ncol; to_val(PyTuple_GetItem(p, 1), &ncol);
      **m = IM::zeros(nrow, ncol);
      casadi_densify(get_ptr(data), (**m).sparsity().T(), (**m).ptr(), true);
      return true;
    }

    bool DM_from_array_conv(GUESTOBJECT *p, casadi::DM** m) {
      if (!m) return true;
      std::vector<double> data;
      if (!to_val(PyTuple_GetItem(p, 2), &data)) return false;
      casadi_int nrow; to_val(PyTuple_GetItem(p, 0), &nrow);
      casadi_int ncol; to_val(PyTuple_GetItem(p, 1), &ncol);
      **m = DM::zeros(nrow, ncol);
      casadi_densify(get_ptr(data), (**m).sparsity().T(), (**m).ptr(), true);
      return true;
    }

    bool DM_from_csc_conv(GUESTOBJECT *p, casadi::DM** m) {
      std::vector<double> data;
      std::vector<casadi_int> colind, row;
      if (!to_val(PyTuple_GetItem(p, 4), &data)) return false;
      if (!to_val(PyTuple_GetItem(p, 3), &row)) return false;
      if (!to_val(PyTuple_GetItem(p, 2), &colind)) return false;
      casadi_int nrow; to_val(PyTuple_GetItem(p, 0), &nrow);
      casadi_int ncol; to_val(PyTuple_GetItem(p, 1), &ncol);
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

    bool is_scalar_np_array(GUESTOBJECT *p) {
      if (PyErr_Occurred()) PyErr_Clear(); // Clear pending exception before type check
      if (PyObject_HasAttrString(p, "__array__")) {
        PyObject *cr = PyObject_GetAttrString(p, (char*) "size");
        if (cr) {
          casadi_int size;
          casadi_int res = to_val(cr, &size);
          Py_DECREF(cr);
          if (!res) return false;
          return size==1;
        } else {
          PyErr_Clear();
          return false;
        }
      }
      return false;
   }

#endif


  } // namespace CasADi
 }

%fragment("casadi_aux", "header", fragment="casadi_decl") {
  namespace casadi {
    template<typename M> bool to_val(GUESTOBJECT *p, M* m) {
      // Copy the pointer
      M *m2 = m;
      bool ret = to_ptr(p, m ? &m2 : 0);
      // If pointer changed, copy the object
      if (m!=m2) *m=*m2;
      return ret;
    }

    // Same as to_ptr, but with GenericType
    template<typename M> bool to_generic(GUESTOBJECT *p, GenericType** m) {
      if (m) {
        // Temporary
        M tmp, *tmp_ptr=&tmp;
        bool ret = to_ptr(p, &tmp_ptr);
        if (!ret) return ret;
        **m = GenericType(*tmp_ptr);
        return ret;
      } else {
        return to_ptr(p, static_cast<M**>(0));
      }
    }

    // Check if casadi_int
    template<typename T> struct is_int {
      static inline bool check() {return false;}
    };

    template<> struct is_int<casadi_int> {
      static inline bool check() {return true;}
    };

    bool is_null(GUESTOBJECT *p) {
#ifdef SWIGPYTHON
      if (p == Py_None) return true;
#endif
#ifdef SWIGMATLAB
      if (p == 0) return true;
#endif
#ifdef SWIGWASMJS
      // EM_VAL == 0 is the uninitialized handle.  Otherwise inspect the
      // val with borrow semantics (take + release ownership) so the
      // caller's refcount stays intact.
      if (!p) return true;
      emscripten::val v = emscripten::val::take_ownership(p);
      bool nul = v.isNull() || v.isUndefined();
      v.release_ownership();
      return nul;
#endif
      return false;
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


%fragment("casadi_bool", "header", fragment="casadi_aux", fragment=SWIG_AsVal_frag(bool)) {
  namespace casadi {
    bool to_ptr(GUESTOBJECT *p, bool** m) {
      // Treat Null
      if (is_null(p)) return false;

      // Standard typemaps
      if (SWIG_IsOK(SWIG_AsVal(bool)(p, m ? *m : 0))) return true;

#ifdef SWIGMATLAB
      if (mxIsLogicalScalar(p)) {
        if (m) **m = mxIsLogicalScalarTrue(p);
        return true;
      }
#endif

      // No match
      return false;
    }

    GUESTOBJECT * from_ptr(const bool *a) {
#ifdef SWIGPYTHON
      return PyBool_FromLong(*a);
#elif defined(SWIGMATLAB)
      return mxCreateLogicalScalar(*a);
#elif defined(SWIGWASMJS)
      // Pack as EM_VAL of a JS boolean.  Used by director paths and by
      // GenericType variant-dispatch when the variant is OT_BOOL.
      emscripten::val v(static_cast<bool>(*a));
      return reinterpret_cast<GUESTOBJECT*>(v.release_ownership());
#else
      return 0;
#endif
    }
  } // namespace casadi
 }

%fragment("casadi_int", "header", fragment="casadi_aux", fragment=SWIG_AsVal_frag(int), fragment=SWIG_AsVal_frag(long), fragment=SWIG_AsVal_frag(long long)) {
  namespace casadi {
    bool to_ptr(GUESTOBJECT *p, casadi_int** m) {
      // Treat Null
      if (is_null(p)) return false;

      // long long
      {
        long long tmp;
        if (SWIG_IsOK(SWIG_AsVal(long long)(p, &tmp))) {
          if (m) **m = static_cast<casadi_int>(tmp);
          return true;
        }
      }

#ifdef SWIGPYTHON
      if (is_scalar_np_array(p)) {
        PyObject *cr = PyObject_CallMethod(p, (char*) "item", 0);
        if (cr) {
          casadi_int res = to_ptr(cr, m);
          Py_DECREF(cr);
          if (!res) return false;
          return true;
        } else {
          PyErr_Clear();
          return false;
        }
      }
#endif // SWIGPYTHON

      bool tmp;
      if (to_val(p, m? &tmp : 0)) {
        if (m) **m = tmp;
        return true;
      }

      // No match
      return false;
    }

    GUESTOBJECT * from_ptr(const casadi_int *a) {
#ifdef SWIGPYTHON
      return PyLong_FromLongLong(*a);
#elif defined(SWIGMATLAB)
      return mxCreateDoubleScalar(static_cast<double>(*a));
#elif defined(SWIGWASMJS)
      // Pack as EM_VAL of a JS Number when in safe-integer range, BigInt
      // otherwise.  Mirrors how the rest of the wasm-js boundary surfaces
      // casadi_int -- JS callers tolerate either form.
      emscripten::val v = (*a >  9007199254740992LL || *a < -9007199254740992LL)
          ? emscripten::val(static_cast<long long>(*a))   // -> BigInt
          : emscripten::val(static_cast<double>(*a));     // -> Number
      return reinterpret_cast<GUESTOBJECT*>(v.release_ownership());
#else
      return 0;
#endif
    }
  } // namespace casadi
 }

%fragment("casadi_double", "header", fragment="casadi_aux", fragment=SWIG_AsVal_frag(double)) {
  namespace casadi {
    bool to_ptr(GUESTOBJECT *p, double** m) {
      // Treat Null
      if (is_null(p)) return false;

      // Standard typemaps
      if (SWIG_IsOK(SWIG_AsVal(double)(p, m ? *m : 0))) return true;

#ifdef SWIGPYTHON
      if (is_scalar_np_array(p)) {
        PyObject *cr = PyObject_CallMethod(p, (char*) "item", 0);
        if (cr) {
          casadi_int res = to_ptr(cr, m);
          Py_DECREF(cr);
          if (!res) return false;
          return true;
        } else {
          PyErr_Clear();
          return false;
        }
      }
#endif // SWIGPYTHON

      casadi_int tmp;
      if (to_val(p, m? &tmp: 0)) {
        if (m) **m = tmp;
        return true;
      }

      // No match
      return false;
    }

    GUESTOBJECT * from_ptr(const double *a) {
#ifdef SWIGPYTHON
      return PyFloat_FromDouble(*a);
#elif defined(SWIGMATLAB)
      return mxCreateDoubleScalar(*a);
#elif defined(SWIGWASMJS)
      emscripten::val v(*a);
      return reinterpret_cast<GUESTOBJECT*>(v.release_ownership());
#else
      return 0;
#endif
    }
  } // namespace casadi
 }


%fragment("casadi_vector", "header", fragment="casadi_aux") {
  namespace casadi {

#ifdef SWIGMATLAB

    // Cell array
    template<typename M> bool to_ptr_cell(GUESTOBJECT *p, std::vector<M>** m) {
      // Cell arrays (only row vectors)
      if (mxGetClassID(p)==mxCELL_CLASS) {
        casadi_int nrow = mxGetM(p), ncol = mxGetN(p);
        if (nrow==1 || (nrow==0 && ncol==0) || ncol==1) {
          casadi_int n = (nrow==0 || ncol==0) ? 0 : std::max(nrow, ncol);
          // Allocate elements
          if (m) {
            (**m).clear();
            (**m).reserve(n);
          }

          // Temporary
          M tmp;

          // Loop over elements
          for (casadi_int i=0; i<n; ++i) {
            // Get element
            mxArray* pe = mxGetCell(p, i);
            if (pe==0) return false;

            // Convert element
            M *m_i = m ? &tmp : 0;
            if (!to_ptr(pe, m_i ? &m_i : 0)) {
              return false;
            }
            if (m) (**m).push_back(*m_i);
          }
          return true;
        }
      }
      return false;
    }

    // MATLAB row/column vector maps to std::vector<double>
    bool to_ptr(GUESTOBJECT *p, std::vector<double> **m) {
      // Treat Null
      if (is_null(p)) return false;

      if (mxIsDouble(p) && mxGetNumberOfDimensions(p)==2
          && (mxGetM(p)<=1 || mxGetN(p)<=1)) {
        if (m) {
          double* data = static_cast<double*>(mxGetData(p));
          casadi_int n = mxGetM(p)*mxGetN(p);
          (**m).resize(n);
          std::copy(data, data+n, (**m).begin());
        }
        return true;
      }

      // Cell array
      if (to_ptr_cell(p, m)) return true;

      // No match
      return false;
    }

    bool to_ptr(GUESTOBJECT *p, std::vector<casadi_int>** m) {
      if (mxIsDouble(p) && mxGetNumberOfDimensions(p)==2
          && (mxGetM(p)<=1 || mxGetN(p)<=1)) {
        double* data = static_cast<double*>(mxGetData(p));
        casadi_int n = mxGetM(p)*mxGetN(p);

        // Check if all integers
        bool all_integers=true;
        for (casadi_int i=0; all_integers && i<n; ++i) {
          if (data[i]!=static_cast<casadi_int>(data[i])) {
            all_integers = false;
            break;
          }
        }

        // Successful conversion
        if (all_integers) {
          if (m) {
            (**m).resize(n);
            std::copy(data, data+n, (**m).begin());
          }
          return true;
        }
      }

      if (mxIsLogical(p) && !mxIsLogicalScalar(p) &&mxGetNumberOfDimensions(p)==2
          && (mxGetM(p)<=1 || mxGetN(p)<=1) ) {
        casadi_int n = mxGetM(p)*mxGetN(p);
        mxLogical* data = static_cast<mxLogical*>(mxGetData(p));
        if (m) {
          (**m).resize(n);
          std::copy(data, data+n, (**m).begin());
        }
        return true;
      }

      // Cell array
      if (to_ptr_cell(p, m)) return true;

      return false;
    }

    bool to_ptr(GUESTOBJECT *p, std::vector<bool>** m) {
      if (mxIsDouble(p) && mxGetNumberOfDimensions(p)==2
          && (mxGetM(p)<=1 || mxGetN(p)<=1)) {
        double* data = static_cast<double*>(mxGetData(p));
        casadi_int n = mxGetM(p)*mxGetN(p);

        // Check if all integers
        bool all_0_or_1 = true;
        for (casadi_int i=0; all_0_or_1 && i<n; ++i) {
          double d = data[i];
          all_0_or_1 = all_0_or_1 && (d==1 || d==0);
        }

        // Successful conversion
        if (all_0_or_1) {
          if (m) {
            (**m).resize(n);
            std::copy(data, data+n, (**m).begin());
          }
          return true;
        }
      }

      if (mxIsLogical(p) && !mxIsLogicalScalar(p) &&mxGetNumberOfDimensions(p)==2
          && (mxGetM(p)<=1 || mxGetN(p)<=1) ) {
        casadi_int n = mxGetM(p)*mxGetN(p);
        mxLogical* data = static_cast<mxLogical*>(mxGetData(p));
        if (m) {
          (**m).resize(n);
          std::copy(data, data+n, (**m).begin());
        }
        return true;
      }

      // Cell array
      if (to_ptr_cell(p, m)) return true;

      return false;
    }

    // MATLAB n-by-m char array mapped to vector of length m
    bool to_ptr(GUESTOBJECT *p, std::vector<std::string>** m) {
      // MATLAB string class (R2016b+): convert via cellstr -> cell of char,
      // then reuse the cell path. (Note: char() on a multi-element string
      // array yields a 3-D char in modern MATLAB, so cellstr is the right
      // bridge here.)
      if (mxIsClass(p, "string")) {
        mxArray *cell_arr = 0;
        mxArray *prhs[1] = { p };
        mxArray *err = mexCallMATLABWithTrap(1, &cell_arr, 1, prhs, "cellstr");
        if (err) {
          mxDestroyArray(err);
          if (cell_arr) mxDestroyArray(cell_arr);
          return false;
        }
        if (!cell_arr) return false;
        bool ok = to_ptr(cell_arr, m);
        mxDestroyArray(cell_arr);
        return ok;
      }

      if (mxIsChar(p)) {
	if (m) {
          // Get data
	  size_t nrow = mxGetM(p);
	  size_t ncol = mxGetN(p);
          mxChar *data = mxGetChars(p);

          // Allocate space for output
          (**m).resize(nrow);
          std::vector<std::string> &m_ref = **m;

          // For all strings
          for (size_t j=0; j!=nrow; ++j) {
            // Get length without trailing spaces
            size_t len = ncol;
            while (len!=0 && data[j + nrow*(len-1)]==' ') --len;

            // Check if null-terminated
            for (size_t i=0; i!=len; ++i) {
              if (data[j + nrow*i]=='\0') {
                len = i;
                break;
              }
            }

            // Create a string of the desired length
            m_ref[j] = std::string(len, ' ');

            // Get string content
            for (size_t i=0; i!=len; ++i) {
              m_ref[j][i] = data[j + nrow*i];
            }
          }
        }
	return true;
      }

      // Cell array
      if (to_ptr_cell(p, m)) return true;

      // No match
      return false;
    }
#endif // SWIGMATLAB

    template<typename M> bool to_ptr(GUESTOBJECT *p, std::vector<M>** m) {
      // Treat Null
      if (is_null(p)) return false;
#ifdef SWIGWASMJS
      {
        emscripten::val v = emscripten::val::take_ownership(
            reinterpret_cast< ::emscripten::EM_VAL >(p));
        bool is_arr = v.isArray();
        bool is_obj = (v.typeOf().as<std::string>() == "object") && !is_arr;
        if (is_obj && v.hasOwnProperty("_ptr")) {
          // Only accept SWIG-tagged proxies whose _swig_type marks
          // them as std::vector<...>.  Without this check, ANY proxy
          // with `_ptr` (Function, MX, Callback, ...) gets reinterpret
          // cast'd as a vector<M>* -- which is what made Callback
          // proxies passed as nlpsol options' iteration_callback get
          // misidentified as OT_BOOLVECTOR.
          emscripten::val swig_t = v["_swig_type"];
          bool is_vec_proxy =
              (swig_t.typeOf().as<std::string>() == "string") &&
              (swig_t.as<std::string>().rfind("_p_std__vectorT_", 0) == 0);
          if (is_vec_proxy) {
            uintptr_t raw = v["_ptr"].as<uintptr_t>();
            v.release_ownership();
            if (m) *m = reinterpret_cast<std::vector<M>*>(raw);
            return true;
          }
          v.release_ownership();
          return false;
        }
        if (is_arr) {
          unsigned n = v["length"].as<unsigned>();
          if (m) (*m)->resize(n);
          bool ok = true;
          for (unsigned i = 0; i < n && ok; ++i) {
            emscripten::val el = v[i];
            ::emscripten::EM_VAL eh = el.as_handle();
            if (!m) {
              // Probe: tagged proxy must validate via SWIG_ConvertPtr
              // with the M descriptor (cast-chain check, no fuzzy).
              // Untagged values (number, array) fall through to fuzzy
              // to_ptr<M> below.
              emscripten::val elv = emscripten::val::take_ownership(eh);
              emscripten::val elt = elv["_swig_type"];
              bool tagged = (elt.typeOf().as<std::string>() == "string");
              elv.release_ownership();
              if (tagged) {
                void* raw = 0;
                M* probe_typed = 0;
                M sentinel;
                probe_typed = &sentinel;
                M* before = probe_typed;
                bool accept = ::casadi::to_ptr(reinterpret_cast<GUESTOBJECT*>(eh), &probe_typed);
                if (!accept || probe_typed == before) ok = false;
              } else {
                M tmp;
                M* tmp_p = &tmp;
                if (!::casadi::to_ptr(reinterpret_cast<GUESTOBJECT*>(eh), &tmp_p)) ok = false;
              }
            } else {
              M tmp;
              M* tmp_p = &tmp;
              if (!::casadi::to_ptr(reinterpret_cast<GUESTOBJECT*>(eh), &tmp_p)) {
                ok = false;
              } else {
                (**m)[i] = *tmp_p;
              }
            }
          }
          v.release_ownership();
          return ok;
        }
        v.release_ownership();
        return false;
      }
#endif
#ifdef SWIGPYTHON

      // Some built-in types are iterable
      if (PyDict_Check(p) || PyString_Check(p) || PySet_Check(p) || PyUnicode_Check(p)) return false;

      // An object exposing a casadi conversion hook is a SINGLE matrix, not
      // a sequence of them -- don't iterate it as a vector.  (A 1-D numeric-
      // like iterable such as the experimental NumpyArray would otherwise be
      // walked element-by-element with float()/int(), which throws on a
      // symbolic element and can crash the interpreter mid-iteration.)
      if (PyObject_HasAttrString(p, "__MX__") || PyObject_HasAttrString(p, "__SX__")
          || PyObject_HasAttrString(p, "__DM__")) return false;

      // Make sure shape is 1D, if defined.
      if (PyErr_Occurred()) PyErr_Clear(); // Clear pending exception before type check
      if (PyObject_HasAttrString(p, "shape")) {
        PyObject * shape = PyObject_GetAttrString(p, "shape");
        if(!PyTuple_Check(shape) || PyTuple_Size(shape)!=1) {
          Py_DECREF(shape);
          return false;
        }
      }

      // Iterator to the sequence
      PyObject *it = PyObject_GetIter(p);
      if (!it) {
        PyErr_Clear();
        return false;
      }

      // Allocate elements
      if (m) (**m).clear();

      // Temporary
      M tmp;

      PyObject *pe;
      // Iterate over sequence
      while ((pe=PyIter_Next(it))) {
        // Convert element
        M *m_i = m ? &tmp : 0;
        if (!to_ptr(pe, m_i ? &m_i : 0)) {
          // Failure
          Py_DECREF(pe);
          Py_DECREF(it);
          return false;
        }
        if (m) (**m).push_back(*m_i);
        Py_DECREF(pe);
      }
      Py_DECREF(it);
      return true;
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
      // Cell array
      if (to_ptr_cell(p, m)) return true;
#endif // SWIGMATLAB
      // No match
      return false;
    }

#ifdef SWIGMATLAB
    GUESTOBJECT* from_ptr(const std::vector<double> *a) {
      mxArray* ret = mxCreateDoubleMatrix(1, a->size(), mxREAL);
      std::copy(a->begin(), a->end(), static_cast<double*>(mxGetData(ret)));
      return ret;
    }
    GUESTOBJECT* from_ptr(const std::vector<casadi_int> *a) {
      mxArray* ret = mxCreateDoubleMatrix(1, a->size(), mxREAL);
      std::copy(a->begin(), a->end(), static_cast<double*>(mxGetData(ret)));
      return ret;
    }
    GUESTOBJECT* from_ptr(const std::vector<bool> *a) {
      mxArray* ret = mxCreateLogicalMatrix(1, a->size());
      std::copy(a->begin(), a->end(), static_cast<bool*>(mxGetData(ret)));
      return ret;
    }
    GUESTOBJECT* from_ptr(const std::vector<std::string> *a) {
      // Collect arguments as char arrays
      std::vector<const char*> str(a->size());
      for (size_t i=0; i<str.size(); ++i) str[i] = (*a)[i].c_str();

      // std::vector<string> maps to MATLAB char array with multiple columns
      return mxCreateCharMatrixFromStrings(str.size(), str.empty() ? 0 : &str[0]);
    }
#endif // SWIGMATLAB
// wasm-js intentionally does NOT define `from_ptr(const std::vector<std::string>*)`:
// the regular return path falls through to the generic
// `from_ptr<std::vector<M>>` template below, which produces a `{_ptr}`
// carrier wrapped on the JS side as the `StringVector` proxy class
// (`__vec_to_arr` drains it).  Adding a wasm-js named overload that
// returned a JS array would break that wrapping (the JS side reads
// `_ptr` from the return and would find none on a plain array).
//
// Director paths instead opt into JS-array form via the
// `%casadi_director_vec_str` typemap (see casadi.i below) -- registered
// LATER than the generic directorin so it wins via last-registered.

#ifdef SWIGWASMJS
    // Per-M JS-class-name trait so from_ptr<std::vector<M>> can route
    // through M.__wrap_<NAME>(ptr) and return a fully-wrapped XVector
    // proxy.  We can't easily go through swig_type_info::clientdata
    // here because std::vector<M>* isn't always registered in
    // swig_type_initial[] (SWIG only emits pointer-type entries for
    // types referenced via $descriptor in a wrapper signature, and
    // most vectors flow as values).  Specialized below by the
    // %wasm_vec(T, NAME) macro for each %casadi_template_vec entry.
    // Unspecialized M's return null -> SWIG_NewPointerObj fallback to
    // bare {_ptr} carrier (the legacy behavior).
    template<typename M> struct wasmjs_vector_jsname {
      static const char* get() { return 0; }
    };
#endif

    template<typename M> GUESTOBJECT* from_ptr(const std::vector<M> *a) {
#ifdef SWIGWASMJS
      std::vector<M>* fresh = new std::vector<M>(*a);
      const char* nm = wasmjs_vector_jsname<M>::get();
      if (nm) return reinterpret_cast<GUESTOBJECT*>(SWIG_WASMJS_WrapByName(fresh, nm));
      return SWIG_NewPointerObj(fresh, 0, SWIG_POINTER_OWN);
#endif
#ifdef SWIGPYTHON
      // std::vector maps to Python list
      PyObject* ret = PyList_New(a->size());
      if (!ret) return 0;
      for (casadi_int k=0; k<a->size(); ++k) {
        PyObject* el = from_ref(a->at(k));
        if (!el) {
          Py_DECREF(ret);
          return 0;
        }
        PyList_SetItem(ret, k, el);
      }
      return ret;
#elif defined(SWIGMATLAB)
      // std::vector maps to MATLAB cell array
      mxArray* ret = mxCreateCellMatrix(1, a->size());
      if (!ret) return 0;
      for (casadi_int k=0; k<a->size(); ++k) {
        mxArray* el = from_ref(a->at(k));
        if (!el) return 0;
        mxSetCell(ret, k, el);
      }
      return ret;
#else
      return 0;
#endif
    }
  } // namespace casadi
}


%fragment("casadi_vectorvector", "header", fragment="casadi_aux") {
  namespace casadi {

#ifdef SWIGMATLAB

    // Cell array
    template<typename M> bool to_ptr_cell2(GUESTOBJECT *p, std::vector< std::vector<M> >** m) {
      // Cell arrays (only row vectors)
      if (mxGetClassID(p)==mxCELL_CLASS) {
        casadi_int nrow = mxGetM(p), ncol = mxGetN(p);
        if (true) {
          // Allocate elements
          if (m) {
            (**m).clear();
            (**m).resize(nrow, std::vector<M>(ncol));
          }

          // Temporary
          M tmp;

          // Loop over elements
          for (casadi_int i=0; i<nrow*ncol; ++i) {
            // Get element
            mxArray* pe = mxGetCell(p, i);
            if (pe==0) return false;

            // Convert element
            M *m_i = m ? &tmp : 0;
            if (!to_ptr(pe, m_i ? &m_i : 0)) {
              return false;
            }

            if (m) (**m)[i % nrow][i / nrow] = tmp;
          }
          return true;
        }
      }
      return false;
    }

#endif // SWIGMATLAB

    template<typename M> bool to_ptr(GUESTOBJECT *p, std::vector< std::vector<M> >** m) {
      // Treat Null
      if (is_null(p)) return false;

      // Pass on to to_ptr(GUESTOBJECT *p, std::vector<M>** m)
      if (to_ptr< std::vector<M> >(p, m)) return true;

#ifdef SWIGMATLAB
      // Cell array
      if (to_ptr_cell2(p, m)) return true;
#endif // SWIGMATLAB
      return false;
    }

  } // namespace casadi
}

%fragment("casadi_function", "header", fragment="casadi_aux") {
  namespace casadi {
    bool to_ptr(GUESTOBJECT *p, Function** m) {
      // Treat Null
      if (is_null(p)) return false;

      // Function already?
      if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                    $descriptor(casadi::Function*), 0))) {
        return true;
      }

      // No match
      return false;
    }

    GUESTOBJECT* from_ptr(const Function *a) {
      return SWIG_NewPointerObj(new Function(*a), $descriptor(casadi::Function *), SWIG_POINTER_OWN);
    }
  } // namespace casadi
}

%fragment("casadi_generictype", "header", fragment="casadi_aux") {
  namespace casadi {
    bool to_ptr(GUESTOBJECT *p, GenericType** m) {
#ifdef SWIGPYTHON
      if (p==Py_None) {
        if (m) **m=GenericType();
        return true;
      }
#endif // SWIGPYTHON

      // Treat Null
      if (is_null(p)) return false;

      // GenericType already?
      if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                    $descriptor(casadi::GenericType*), 0))) {
        return true;
      }

      // Try to convert to different types
      if (to_generic<bool>(p, m)
          || to_generic<casadi_int>(p, m)
          || to_generic<double>(p, m)
          || to_generic<std::string>(p, m)
          || to_generic<std::vector<bool> >(p, m)
          || to_generic<std::vector<casadi_int> >(p, m)
          || to_generic<std::vector<double> >(p, m)
          || to_generic<std::vector<std::string> >(p, m)
          || to_generic<std::vector<std::vector<std::string> > >(p, m)
          || to_generic<std::vector<std::vector<casadi_int> > >(p, m)
          || to_generic<std::vector<std::vector<double> > >(p, m)
          || to_generic<casadi::Function>(p, m)
          || to_generic<std::vector<casadi::Function> >(p, m)
          || to_generic<casadi::GenericType::Dict>(p, m)
          || to_generic<std::vector<casadi::GenericType::Dict> >(p, m)
          || to_generic<std::vector<std::vector<casadi::GenericType> > >(p, m)
          || to_generic<std::vector<casadi::GenericType> >(p, m)) {
        return true;
      }

      // Check if it can be converted to boolean (last as e.g. can be converted to boolean)
      if (to_generic<bool>(p, m)) return true;

      // No match
      return false;
    }

    GUESTOBJECT * from_ptr(const GenericType *a) {
      switch (a->getType()) {
      case OT_BOOL: return from_tmp(a->as_bool());
      case OT_INT: return from_tmp(a->as_int());
      case OT_DOUBLE: return from_tmp(a->as_double());
      case OT_STRING: return from_tmp(a->as_string());
      case OT_INTVECTOR: return from_tmp(a->as_int_vector());
      case OT_INTVECTORVECTOR: return from_tmp(a->as_int_vector_vector());
      case OT_BOOLVECTOR: return from_tmp(a->as_bool_vector());
      case OT_DOUBLEVECTOR: return from_tmp(a->as_double_vector());
      case OT_DOUBLEVECTORVECTOR: return from_tmp(a->as_double_vector_vector());
      case OT_STRINGVECTOR: return from_tmp(a->as_string_vector());
      case OT_STRINGVECTORVECTOR: return from_tmp(a->as_string_vector_vector());
      case OT_DICT: return from_tmp(a->as_dict());
      case OT_DICTVECTOR: return from_tmp(a->as_dict_vector());
      case OT_VECTORVECTOR: return from_tmp(a->as_vector_vector());
      case OT_VECTOR: return from_tmp(a->as_vector());
      case OT_FUNCTION: return from_tmp(a->as_function());
      case OT_FUNCTIONVECTOR: return from_tmp(a->as_function_vector());
#ifdef SWIGPYTHON
      case OT_NULL:
      case OT_VOIDPTR:
        return Py_None;
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
      case OT_VOIDPTR:
        return mxCreateDoubleScalar(0);
#endif
      default: return 0;
      }
    }
  } // namespace casadi
}

%fragment("casadi_string", "header", fragment="casadi_aux") {
  namespace casadi {
    bool to_ptr(GUESTOBJECT *p, std::string** m) {
      // Treat Null
      if (is_null(p)) return false;

      // String already?
      if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                    $descriptor(std::string*), 0))) {
        return true;
      }

#ifdef SWIGWASMJS
      // Phase 3.4: accept JS string primitive.
      {
        emscripten::val v = emscripten::val::take_ownership(
            reinterpret_cast< ::emscripten::EM_VAL >(p));
        bool is_str = (v.typeOf().as<std::string>() == "string");
        if (is_str) {
          if (m) **m = v.as<std::string>();
          v.release_ownership();
          return true;
        }
        v.release_ownership();
      }
#endif

#ifdef SWIGPYTHON
      if (PyString_Check(p) || PyUnicode_Check(p)) {
        if (m) {
          (*m)->clear();
          (*m)->append(python_string_to_std_string(p));
        }
        return true;
      }
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
      // MATLAB string class (R2016b+): the MEX wrapper is a 1x1 opaque
      // container, so we cannot read element count via mxGetNumberOfElements.
      // Bridge via cellstr -> cell of char, then require exactly one element
      // for the scalar std::string slot.
      if (mxIsClass(p, "string")) {
        mxArray *cell_arr = 0;
        mxArray *prhs[1] = { p };
        mxArray *err = mexCallMATLABWithTrap(1, &cell_arr, 1, prhs, "cellstr");
        if (err) {
          mxDestroyArray(err);
          if (cell_arr) mxDestroyArray(cell_arr);
          return false;
        }
        if (!cell_arr) return false;
        if (mxGetClassID(cell_arr) != mxCELL_CLASS
            || mxGetNumberOfElements(cell_arr) != 1) {
          mxDestroyArray(cell_arr);
          return false;
        }
        mxArray *char_arr = mxGetCell(cell_arr, 0);
        bool ok = char_arr && to_ptr(char_arr, m);
        mxDestroyArray(cell_arr);
        return ok;
      }
      if (mxIsChar(p) && mxGetM(p)<=1) {
        if (m) {
          if (mxGetM(p)==0) return true;
          size_t len=mxGetN(p);
          std::vector<char> s(len+1);
          if (mxGetString(p, &s[0], (len+1)*sizeof(char))) {
            casadi_warning("mxGetString returned NULL");
            return false;
          }
          // Matlab silent failure; see #4034
          if (s[0]=='\0' && len>0) {
            casadi_warning("mxGetString failure, see https://github.com/casadi/casadi/issues/4034");
            return false;
          }
          **m = std::string(&s[0], len);
        }
        return true;
      }
#endif // SWIGMATLAB

      // No match
      return false;
    }

    GUESTOBJECT* from_ptr(const std::string *a) {
#ifdef SWIGPYTHON
      return PyString_FromString(a->c_str());
#elif defined(SWIGMATLAB)
      return mxCreateString(a->c_str());
#elif defined(SWIGWASMJS)
      emscripten::val v(*a);  // emscripten::val has a std::string ctor (UTF-8)
      return reinterpret_cast<GUESTOBJECT*>(v.release_ownership());
#else
      return 0;
#endif
    }
  } // namespace casadi
}

%fragment("casadi_slice", "header", fragment="casadi_aux") {
  namespace casadi {
    bool to_ptr(GUESTOBJECT *p, Slice** m) {
      // Treat Null
      if (is_null(p)) return false;

      // Slice already?
      if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                    $descriptor(casadi::Slice*), 0))) {
        return true;
      }

#ifdef SWIGPYTHON

      // Python casadi_int
      if (PyInt_Check(p)) {
        if (m) {
          (**m).start = PyInt_AsLong(p);
          (**m).stop = (**m).start+1;
          if ((**m).stop==0) (**m).stop = std::numeric_limits<casadi_int>::max();
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
          // Map sentinel values from PySlice_Unpack to CasADi's limits
          // PySlice_Unpack returns PY_SSIZE_T_MIN or PY_SSIZE_T_MAX as sentinels
          // depending on step direction, so check both extremes
          (**m).start = (start == PY_SSIZE_T_MIN || start == PY_SSIZE_T_MAX)
              ? std::numeric_limits<casadi_int>::min()
              : static_cast<casadi_int>(start);
          (**m).stop = (stop == PY_SSIZE_T_MAX || stop == PY_SSIZE_T_MIN)
              ? std::numeric_limits<casadi_int>::max()
              : static_cast<casadi_int>(stop);
          if (step != 1) {
            (**m).step = static_cast<casadi_int>(step);
          }
        }
        return true;
      }
#endif // SWIGPYTHON

      // No match
      return false;
    }

    GUESTOBJECT* from_ptr(const Slice *a) {
      return SWIG_NewPointerObj(new Slice(*a), $descriptor(casadi::Slice *), SWIG_POINTER_OWN);
    }

  } // namespace casadi
}

%fragment("casadi_map", "header", fragment="casadi_aux") {
  namespace casadi {
    template<typename M> bool to_ptr(GUESTOBJECT *p, std::map<std::string, M>** m) {
#ifdef SWIGWASMJS
      if (is_null(p)) return true;
      {
        emscripten::val v = emscripten::val::take_ownership(
            reinterpret_cast< ::emscripten::EM_VAL >(p));
        std::string ty = v.typeOf().as<std::string>();
        bool is_obj = (ty == "object") && !v.isArray() && !v.hasOwnProperty("_ptr");
        if (!is_obj) { v.release_ownership(); return false; }
        emscripten::val keys = emscripten::val::global("Object").call<emscripten::val>("keys", v);
        unsigned n = keys["length"].as<unsigned>();
        bool ok = true;
        for (unsigned i = 0; i < n && ok; ++i) {
          emscripten::val key = keys[i];
          emscripten::val val = v[key];
          ::emscripten::EM_VAL vh = val.as_handle();
          if (m) {
            std::string ks = key.as<std::string>();
            M tmp; M* tmp_p = &tmp;
            if (!casadi::to_ptr(reinterpret_cast<GUESTOBJECT*>(vh), &tmp_p)) {
              ok = false;
            } else {
              (**m)[ks] = *tmp_p;
            }
          } else {
            if (!casadi::to_ptr(reinterpret_cast<GUESTOBJECT*>(vh), (M**)0)) ok = false;
          }
        }
        v.release_ownership();
        return ok;
      }
#endif
#ifdef SWIGPYTHON
      if (PyDict_Check(p)) {
        PyObject *key, *value;
        Py_ssize_t pos = 0;
        while (PyDict_Next(p, &pos, &key, &value)) {
          if (!(PyString_Check(key) || PyUnicode_Check(key))) return false;
          if (m) {
            M *v=&(**m)[python_string_to_std_string(key)], *v2=v;
            if (!casadi::to_ptr(value, &v)) return false;
            if (v!=v2) *v2=*v; // if only pointer changed
          } else {
            if (!casadi::to_ptr(value, static_cast<M**>(0))) return false;
          }
        }
        return true;
      }
#elif defined(SWIGMATLAB)
      if (mxIsStruct(p) && mxGetM(p)==1 && mxGetN(p)==1) {
	casadi_int len = mxGetNumberOfFields(p);
	for (casadi_int k=0; k<len; ++k) {
	  mxArray *value = mxGetFieldByNumber(p, 0, k);
          if (m) {
	    M *v=&(**m)[std::string(mxGetFieldNameByNumber(p, k))], *v2=v;
            if (!casadi::to_ptr(value, &v)) return false;
            if (v!=v2) *v2=*v; // if only pointer changed
	  } else {
            if (!casadi::to_ptr(value, static_cast<M**>(0))) return false;
	  }
	}
        return true;
      }
#endif
      return false;
    }

    template<typename M> GUESTOBJECT* from_ptr(const std::map<std::string, M> *a) {
#ifdef SWIGPYTHON
      PyObject *p = PyDict_New();
      for (typename std::map<std::string, M>::const_iterator it=a->begin(); it!=a->end(); ++it) {
        PyObject * e = from_ptr(&it->second);
        if (!e) {
          Py_DECREF(p);
          return 0;
        }
        PyDict_SetItemString(p, it->first.c_str(), e);
        Py_DECREF(e);
      }
      return p;
#elif defined(SWIGMATLAB)
      // Get vectors of the field names and mxArrays
      std::vector<const char*> fieldnames;
      std::vector<mxArray*> fields;
      for (typename std::map<std::string, M>::const_iterator it=a->begin(); it!=a->end(); ++it) {
	fieldnames.push_back(it->first.c_str());
	mxArray* f = from_ptr(&it->second);
	if (!f) {
	  // Deallocate elements created up to now
	  for (casadi_int k=0; k<fields.size(); ++k) mxDestroyArray(fields[k]);
	  return 0;
	}
	fields.push_back(f);
      }

      // Create return object
      mxArray *p = mxCreateStructMatrix(1, 1, fields.size(),
					fieldnames.empty() ? 0 : &fieldnames[0]);
      for (casadi_int k=0; k<fields.size(); ++k) mxSetFieldByNumber(p, 0, k, fields[k]);
      return p;
#elif defined(SWIGWASMJS)
      // Pack as EM_VAL of a JS object {key: value, ...}.  For casadi
      // value types (DM/MX/SX/...), from_ptr<M> calls
      // SWIG_NewPointerObj with the type descriptor, which routes
      // through SWIG_WASMJS_NewPointerObj's `type->clientdata` lookup
      // (set by SWIG_WASMJS_init_cast_chains in wasm_js.cxx) to
      // M.__wrap_<JsName>(ptr) -- yielding a fully-wrapped proxy.
      // So the assembled dict has typed-proxy values, no JS-side
      // post-wrap dance required.  Mirrors matlab.cxx where
      // SWIG_Matlab_NewPointerObj produces a full <Class> object
      // via mexCallMATLAB("<Class>_construct", ptr).
      emscripten::val obj = emscripten::val::object();
      for (typename std::map<std::string, M>::const_iterator it=a->begin(); it!=a->end(); ++it) {
        GUESTOBJECT* e = from_ptr(&it->second);
        if (!e) return 0;
        emscripten::val v = emscripten::val::take_ownership(
            reinterpret_cast< ::emscripten::EM_VAL >(e));
        obj.set(it->first, v);
      }
      return reinterpret_cast<GUESTOBJECT*>(obj.release_ownership());
#else
      return 0;
#endif
    }
  } // namespace casadi
}

%fragment("casadi_pair", "header", fragment="casadi_aux") {
  namespace casadi {
#ifdef SWIGMATLAB
    bool to_ptr(GUESTOBJECT *p, std::pair<casadi_int, casadi_int>** m) {
      // (casadi_int,casadi_int) mapped to 2-by-1 double matrix
      if (mxIsDouble(p) && mxGetNumberOfDimensions(p)==2 && !mxIsSparse(p)
          && mxGetM(p)==1 && mxGetN(p)==2) {
        double* data = static_cast<double*>(mxGetData(p));
        casadi_int first = static_cast<casadi_int>(data[0]);
        casadi_int second = static_cast<casadi_int>(data[1]);
        if (data[0]==first && data[1]==second) {
          if (m) **m = std::make_pair(first, second);
          return true;
        } else {
          return false;
        }
      }

      // No match
      return false;
    }
#endif // SWIGMATLAB

    template<typename M1, typename M2> bool to_ptr(GUESTOBJECT *p, std::pair<M1, M2>** m) {
#ifdef SWIGPYTHON
      if (PyTuple_Check(p) && PyTuple_Size(p)==2) {
        PyObject *p_first = PyTuple_GetItem(p, 0);
        PyObject *p_second = PyTuple_GetItem(p, 1);
	return to_val(p_first, m ? &(**m).first : 0)
	  && to_val(p_second, m ? &(**m).second : 0);
      }
#elif defined(SWIGMATLAB)
      // Other overloads mapped to 2-by-1 cell array
      if (mxGetClassID(p)==mxCELL_CLASS && mxGetM(p)==1 && mxGetN(p)==2) {
        mxArray *p_first = mxGetCell(p, 0);
        mxArray *p_second = mxGetCell(p, 1);
        return to_val(p_first, m ? &(**m).first : 0)
          && to_val(p_second, m ? &(**m).second : 0);
      }
#endif
      // No match
      return false;
    }

#ifdef SWIGMATLAB
    GUESTOBJECT* from_ptr(const std::pair<casadi_int, casadi_int>* a) {
      // (casadi_int,casadi_int) mapped to 2-by-1 double matrix
      mxArray* ret = mxCreateDoubleMatrix(1, 2, mxREAL);
      double* data = static_cast<double*>(mxGetData(ret));
      data[0] = a->first;
      data[1] = a->second;
      return ret;
    }
#endif // SWIGMATLAB

    template<typename M1, typename M2> GUESTOBJECT* from_ptr(const std::pair<M1, M2>* a) {
#ifdef SWIGPYTHON
      PyObject* ret = PyTuple_New(2);
      PyTuple_SetItem(ret, 0, from_ref(a->first));
      PyTuple_SetItem(ret, 1, from_ref(a->second));
      return ret;
#elif defined(SWIGMATLAB)
      // Other overloads mapped to 2-by-1 cell array
      mxArray* ret = mxCreateCellMatrix(1, 2);
      mxSetCell(ret, 0, from_ref(a->first));
      mxSetCell(ret, 1, from_ref(a->second));
      return ret;
#else
      return 0;
#endif // SWIGPYTHON
    }
  } // namespace casadi
 }

%fragment("casadi_sx", "header", fragment="casadi_aux") {
  namespace casadi {
    bool to_ptr(GUESTOBJECT *p, SX** m) {
      // Treat Null
      if (is_null(p)) return false;

      // SX already?
      if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                    $descriptor(casadi::Matrix<casadi::SXElem>*), 0))) {
        return true;
      }

#ifdef SWIGPYTHON
      // Object has __SX__ method: honour it BEFORE the scalar/array/vector
      // heuristics below.  A custom iterable-numeric object (e.g. the
      // experimental NumpyArray) that exposes a 1-D `shape` + `__iter__` +
      // `__float__` would otherwise be dragged into the vector<double>
      // fallback, which calls float() on each (possibly symbolic) element.
      if (PyErr_Occurred()) PyErr_Clear(); // Clear pending exception before type check
      if (PyObject_HasAttrString(p,"__SX__")) {
        PyObject *cr = PyObject_CallMethod(p, (char*) "__SX__", 0);
        if (!cr) return false;
        if (cr == Py_None) { Py_DECREF(cr); }
        else {
          // to_val COPIES into *m; to_ptr would only point *m INTO cr,
          // which the Py_DECREF then frees (use-after-free).
          casadi_int flag = to_val(cr, m ? *m : 0);
          Py_DECREF(cr);
          if (flag) return true;
        }
      }
#endif // SWIGPYTHON

      // Try first converting to a temporary DM
      {
        DM tmp;
        if(to_val(p, m? &tmp: 0)) {
          if (m) **m = tmp;
          return true;
        }
      }

#ifdef SWIGPYTHON
      // Clear any error raised by a failed DM conversion (e.g. a __DM__
      // that threw) so it doesn't leak past this no-match path.
      if (PyErr_Occurred()) PyErr_Clear();
      // Numpy arrays will be cast to dense SX
      if (SX_from_array(p, m)) return true;
#endif // SWIGPYTHON

      // No match
      return false;
    }

    GUESTOBJECT* from_ptr(const SX *a) {
      return SWIG_NewPointerObj(new SX(*a), $descriptor(casadi::Matrix<casadi::SXElem> *), SWIG_POINTER_OWN);
    }
  } // namespace casadi
 }

%fragment("casadi_sxelem", "header", fragment="casadi_aux") {
  namespace casadi {
    bool to_ptr(GUESTOBJECT *p, SXElem** m) {
      // Treat Null
      if (is_null(p)) return false;

      // Try first converting to a temporary SX
      {
        SX tmp, *mt=&tmp;
        if(casadi::to_ptr(p, m ? &mt : 0)) {
          if (m && !mt->is_scalar()) return false;
          if (m) **m = mt->scalar();
          return true;
        }
      }

      // No match
      return false;
    }

    GUESTOBJECT* from_ptr(const SXElem *a) {
      return from_ref(SX(*a));
    }
  } // namespace casadi
 }

%fragment("casadi_mx", "header", fragment="casadi_decl") {
  namespace casadi {
    bool to_ptr(GUESTOBJECT *p, MX** m) {
      // Treat Null
      if (is_null(p)) return false;

      // MX already?
      if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                    $descriptor(casadi::MX*), 0))) {
        return true;
      }

#ifdef SWIGPYTHON
      // Object has __MX__ method: honour it BEFORE the DM/scalar/vector
      // heuristics (see the __SX__ note in to_ptr(SX**)).
      if (PyErr_Occurred()) PyErr_Clear(); // Clear pending exception before type check
      if (PyObject_HasAttrString(p,"__MX__")) {
        PyObject *cr = PyObject_CallMethod(p, (char*) "__MX__", 0);
        if (!cr) return false;
        if (cr == Py_None) { Py_DECREF(cr); }
        else {
          // to_val COPIES into *m (see the __SX__ note in to_ptr(SX**)).
          casadi_int flag = to_val(cr, m ? *m : 0);
          Py_DECREF(cr);
          if (flag) return true;
        }
      }
#endif // SWIGPYTHON

      // Try first converting to a temporary DM
      {
        DM tmp;
        if(to_val(p, m ? &tmp : 0)) {
          if (m) **m = tmp;
          return true;
        }
      }

#ifdef SWIGPYTHON
      // Clear any error raised by a failed DM conversion (e.g. a __DM__
      // that threw) so it doesn't leak past this no-match path.
      if (PyErr_Occurred()) PyErr_Clear();
#endif // SWIGPYTHON

      // No match
      return false;
    }

    GUESTOBJECT* from_ptr(const MX *a) {
      return SWIG_NewPointerObj(new MX(*a), $descriptor(casadi::MX*), SWIG_POINTER_OWN);
    }
  } // namespace casadi
 }

%fragment("casadi_dmatrix", "header", fragment="casadi_aux") {
  namespace casadi {
#ifdef SWIGPYTHON
    /** Check PyObjects by class name */
    bool PyObjectHasClassName(PyObject* p, const char * name) {
      PyObject * classo = PyObject_GetAttrString( p, "__class__");
      PyObject * classname = PyObject_GetAttrString( classo, "__name__");

      bool ret = python_string_to_std_string(classname) == name;
      Py_DECREF(classo);Py_DECREF(classname);
      return ret;
    }
#endif // SWIGPYTHON

    bool to_ptr(GUESTOBJECT *p, DM** m) {
      // Treat Null
      if (is_null(p)) return false;

      // DM already?
      if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                    $descriptor(casadi::Matrix<double>*), 0))) {
        return true;
      }

      // Object is a sparsity pattern
      {
        Sparsity *m2;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2),
                                      $descriptor(casadi::Sparsity*), 0))) {
          if (m) **m=DM::ones(*m2);
          return true;
        }
      }

      // Double scalar
      {
        double tmp;
        if (to_val(p, m? &tmp: 0)) {
          if (m) **m=tmp;
          return true;
        }
      }

#ifdef SWIGPYTHON
      // Object has __DM__ method.  A returned None means "I am explicitly
      // not a DM" (e.g. a symbolic NumpyArray) -- stop here rather than
      // falling into the float()-per-element vector fallback below.
      if (PyErr_Occurred()) PyErr_Clear(); // Clear pending exception before type check
      if (PyObject_HasAttrString(p,"__DM__")) {
        char name[] = "__DM__";
        PyObject *cr = PyObject_CallMethod(p, name, 0);
        if (!cr) return false;
        if (cr == Py_None) { Py_DECREF(cr); return false; }
        casadi_int result = to_val(cr, m ? *m : 0);
        Py_DECREF(cr);
        return result;
      }

      if (DM_from_array(p, m)) return true;

      if (DM_from_csc(p,m)) return true;

      {
        std::vector <double> t;
        casadi_int res = to_val(p, &t);
        if (t.size()>0) {
          if (m) **m = casadi::Matrix<double>(t);
        } else {
          if (m) **m = casadi::Matrix<double>(0,0);
        }
        return res;
      }
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
      // MATLAB double matrix (sparse or dense)
      if (mxIsDouble(p) && mxGetNumberOfDimensions(p)==2) {
        if (m) {
          **m = casadi::DM(get_sparsity(p));
          double* data = static_cast<double*>(mxGetData(p));
          casadi_copy(data, (*m)->nnz(), (*m)->ptr());
        }
        return true;
      }
#endif // SWIGMATLAB
#ifdef SWIGWASMJS
      // Phase 3.4: JS Array of numbers -> column-vector DM.
      // 2D nested array [[a,b], [c,d]] -> 2x2 DM (row-major in JS).
      // Borrow the val (take + release) to preserve refcount.
      {
        emscripten::val v = emscripten::val::take_ownership(
            reinterpret_cast< ::emscripten::EM_VAL >(p));
        bool is_arr = v.isArray();
        if (is_arr) {
          unsigned n = v["length"].as<unsigned>();
          // Sniff: nested-array if first element is itself an array.
          bool nested = false;
          if (n > 0) {
            emscripten::val first = v[0];
            nested = first.isArray();
          }
          if (nested) {
            // 2D: each top-level element is a row.  Build a flat
            // column-major data buffer (casadi DM stores cols major)
            // then construct a dense matrix.
            unsigned ncols = 0;
            bool ok = true;
            for (unsigned i = 0; i < n && ok; ++i) {
              emscripten::val row = v[i];
              if (!row.isArray()) { ok = false; break; }
              unsigned w = row["length"].as<unsigned>();
              if (i == 0) ncols = w;
              else if (w != ncols) { ok = false; break; }
            }
            if (ok) {
              std::vector<double> data((size_t)n * (size_t)ncols);
              for (unsigned i = 0; i < n && ok; ++i) {
                emscripten::val row = v[i];
                for (unsigned j = 0; j < ncols && ok; ++j) {
                  emscripten::val el = row[j];
                  std::string ty = el.typeOf().as<std::string>();
                  if (ty == "number" || ty == "bigint") {
                    /* Column-major: data[j*n + i]. */
                    data[(size_t)j * (size_t)n + (size_t)i] = el.as<double>();
                  } else ok = false;
                }
              }
              v.release_ownership();
              if (ok) {
                if (m) **m = casadi::Matrix<double>::reshape(
                    casadi::Matrix<double>(data), n, ncols);
                return true;
              }
              return false;
            }
            v.release_ownership();
            return false;
          }
          // 1D: column vector.
          std::vector<double> data(n);
          bool ok = true;
          for (unsigned i = 0; i < n && ok; ++i) {
            emscripten::val el = v[i];
            std::string ty = el.typeOf().as<std::string>();
            if (ty == "number" || ty == "bigint") data[i] = el.as<double>();
            else ok = false;
          }
          v.release_ownership();
          if (ok) {
            if (m) **m = casadi::Matrix<double>(data);
            return true;
          }
          return false;
        }
        v.release_ownership();
      }
#endif // SWIGWASMJS

      // No match
      return false;
    }

    GUESTOBJECT* from_ptr(const DM *a) {
      return SWIG_NewPointerObj(new DM(*a), $descriptor(casadi::Matrix<double>*), SWIG_POINTER_OWN);
    }
  } // namespace casadi
}

%fragment("casadi_sparsity", "header", fragment="casadi_aux") {
  namespace casadi {
    bool to_ptr(GUESTOBJECT *p, Sparsity** m) {
      // Treat Null
      if (is_null(p)) return false;

      // Sparsity already?
      if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                    $descriptor(casadi::Sparsity*), 0))) {
        return true;
      }

      // No match
      return false;
    }

    GUESTOBJECT* from_ptr(const Sparsity *a) {
      return SWIG_NewPointerObj(new Sparsity(*a), $descriptor(casadi::Sparsity*), SWIG_POINTER_OWN);
    }
  } // namespace casadi
}

%fragment("casadi_imatrix", "header", fragment="casadi_aux", fragment=SWIG_AsVal_frag(int)) {
  namespace casadi {
    bool to_ptr(GUESTOBJECT *p, IM** m) {
      // Treat Null
      if (is_null(p)) return false;

      // Object is a sparsity pattern
      {
        Sparsity *m2;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2),
                                      $descriptor(casadi::Sparsity*), 0))) {
          if (m) **m=IM::ones(*m2);
          return true;
        }
      }

      // First convert to integer
      {
        casadi_int tmp;
        if (to_val(p, m? &tmp: 0)) {
          if (m) **m=tmp;
          return true;
        }
      }

#ifdef SWIGPYTHON
      // Numpy arrays will be cast to dense Matrix<casadi_int>
      if (IM_from_array(p, m)) return true;

      if (PyErr_Occurred()) PyErr_Clear(); // Clear pending exception before type check
      if (PyObject_HasAttrString(p,"__IM__")) {
        PyObject *cr = PyObject_CallMethod(p, (char*) "__IM__", 0);
        if (!cr) return false;
        casadi_int result = to_val(cr, m ? *m : 0);
        Py_DECREF(cr);
        return result;
      }

      {
        std::vector <casadi_int> t;
        if (to_val(p, &t)) {
          if (m) **m = casadi::Matrix<casadi_int>(t);
          return true;
        }
      }
#endif // SWIGPYTHON

#ifdef SWIGMATLAB
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
            **m = casadi::IM(get_sparsity(p));
            for (size_t i=0; i<sz; ++i) {
              (**m)->at(i) = casadi_int(data[i]);
            }
          }
          return true;
        }
      }
#endif // SWIGMATLAB

      // Convert from DM
      {
        DM tmp;
        if (to_val(p, m? &tmp: 0)) {
          // Check integrality
          for (double d : tmp.nonzeros()) {
            if (d!=casadi_int(d)) return false;
          }
          // Convert
          if (m) {
            **m = casadi::Matrix<double>(tmp);
          }
          return true;
        }
      }

      // No match
      return false;
    }
    GUESTOBJECT* from_ptr(const IM *a) {
      DM tmp(*a);
      return from_ref(tmp);
    }

  } // namespace casadi
 }

// Can be overloaded by specifying before importing casadi.i
%fragment("casadi_extra", "header") {}

// Collect all fragments
%fragment("casadi_all", "header", fragment="casadi_aux,casadi_extra,casadi_bool,casadi_int,casadi_double,casadi_vector,casadi_vectorvector,casadi_function,casadi_generictype,casadi_string,casadi_slice,casadi_map,casadi_pair,casadi_sx,casadi_sxelem,casadi_mx,casadi_dmatrix,casadi_sparsity,casadi_imatrix") { }

#endif // SWIGXML

 // Input typemaps.
 //   xStubIn  becomes the pystub_in=  PEP-484 annotation (for `-stubs`
 //           emission by python.cxx).
 //   xTsStub  becomes the tsstub_in=  TypeScript annotation (for `-stubs`
 //           emission by wasm_js.cxx).  Pasted verbatim into the .d.ts;
 //           no translation layer.  Use TS syntax (e.g. "bigint",
 //           "string", "MX[]", "Record<string, GenericType>").
%define %casadi_input_typemaps(xName, xStubIn, xTsStub, xPrec, xType...)
#ifdef SWIGWASMJS
 // wasm_js: all wrapped types cross the wasm boundary as EM_VAL -- an
 // i32 handle into Embind's value table, the wasm-js analog of mxArray*
 // and PyObject*.  to_ptr / SWIG_ConvertPtr unwrap the JS-side carrier
 // (a `{_ptr: <int>}` proxy) via emscripten::val on the wasm side.
%typemap(ctype) xType, const xType&, xType& "EM_VAL"
#endif
 // Pass input by value, check if matches
%typemap(typecheck, noblock=1, precedence=xPrec, fragment="casadi_all") xType {
  $1 = casadi::to_ptr($input, static_cast< xType **>(0));
 }

 // Directorout typemap; as input by value.
 //
 // Works for all targets: wasm-js, matlab, python.  `casadi::to_val`
 // routes through `to_ptr` which has per-target SWIGWASMJS / SWIGPYTHON /
 // SWIGMATLAB branches.  `%dirout_fail` is defined for wasm-js in
 // `Lib/wasm_js/wasm_js.swg` and for matlab/python via
 // `<typemaps/swigtypemaps.swg>`.
 //
 // Per-target specialisations that need a richer body (e.g. wasm-js
 // wrapping class returns via M.__wrap_<C>) override this via
 // last-registered-wins; see `%casadi_director_class` below.
%typemap(directorout, noblock=1, fragment="casadi_all") xType {
    if (!casadi::to_val($input, &$result)) {
      %dirout_fail(SWIG_TypeError,"$type");
    }
 }

 // Pass input by value, convert argument
%typemap(in, doc=xName, pystub_in=xStubIn, tsstub_in=xTsStub, noblock=1, fragment="casadi_all") xType {
  if (!casadi::to_val($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input $argnum to type '" xName "'.");
 }

 // Pass input by value, cleanup
%typemap(freearg, noblock=1) xType {}

 // Pass input by reference, check if matches
%typemap(typecheck, noblock=1, precedence=xPrec, fragment="casadi_all") const xType& {
  $1 = casadi::to_ptr($input, static_cast< xType **>(0));
 }

 // Pass input by reference, convert argument
%typemap(in, doc=xName, pystub_in=xStubIn, tsstub_in=xTsStub, noblock=1, fragment="casadi_all") const xType & (xType m) {
  $1 = &m;
  if (!casadi::to_ptr($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input $argnum to type '" xName "'.");
 }

 // Pass input by reference, cleanup
%typemap(freearg, noblock=1) const xType & {}

#ifdef SWIGWASMJS
 // Non-const reference: use the standard SWIG convention.  SWIG_ConvertPtr
 // (Lib/wasm_js/wasm_jsrun.swg) wraps the EM_VAL via emscripten::val,
 // reads the `_ptr` property and writes the C++ pointer through &$1.
 // Mirrors MATLAB's default SWIGTYPE & typemap.
%typemap(in, doc=xName, pystub_in=xStubIn, tsstub_in=xTsStub, noblock=1) xType & {
  if (!SWIG_IsOK(SWIG_ConvertPtr($input, (void**)&$1, $1_descriptor, 0))) {
    SWIG_exception_fail(SWIG_TypeError, "Failed to convert input $argnum to type '" xName "'.");
  }
 }
%typemap(freearg, noblock=1) xType & {}
#endif

%enddef

 // Output typemaps.
 //   xStubOut becomes pystub_out= on return typemaps and &OUTPUT/&INOUT.
 //   xTsStub  becomes tsstub_out= -- single arg used for both input
 //           (e.g. INOUT pass-through) and output positions.  wasm_js
 //           reads tsstub_out for return types directly.
%define %casadi_output_typemaps(xName, xStubIn, xStubOut, xTsStub, xType...)

#ifdef SWIGWASMJS
 // wasm_js: returns also flow as EM_VAL.  from_ref / from_ptr build the
 // JS-side carrier via SWIG_NewPointerObj; js_marshal_return on the JS
 // side unwraps it and constructs the proxy class.
%typemap(ctype) xType, const xType "EM_VAL"
%typemap(ctype) const xType& "EM_VAL"
#endif

 // Return-by-value
%typemap(out, doc=xName, pystub_out=xStubOut, tsstub_out=xTsStub, noblock=1, fragment="casadi_all") xType, const xType {
  if(!($result = casadi::from_ref($1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type '" xName "'.");
}

// Return a const-ref behaves like return-by-value
%typemap(out, doc=xName, pystub_out=xStubOut, tsstub_out=xTsStub, noblock=1, fragment="casadi_all") const xType& {
  if(!($result = casadi::from_ptr($1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type '" xName "'.");
}

// Inputs marked OUTPUT are also returned by the function, ...
%typemap(argout, noblock=1,fragment="casadi_all") xType &OUTPUT {
  %append_output(casadi::from_ptr($1));
 }

// ... and the corresponding inputs are ignored.
%typemap(in, doc=xName, pystub_out=xStubOut, tsstub_out=xTsStub, noblock=1, numinputs=0) xType &OUTPUT (xType m) {
 $1 = &m;
}

 // Directorin typemap; as output.
 //
 // Works for all targets: wasm-js, matlab, python.  Routes through
 // `casadi::from_ref` / `casadi::from_ptr` which have per-target
 // SWIGWASMJS / SWIGPYTHON / SWIGMATLAB implementations.  Targets where
 // the resulting wire form needs domain-specific wrapping (e.g. wasm-js
 // class types want a real JS proxy via `M.__wrap_<C>`, not the raw
 // `{_ptr}` carrier `from_ref` builds) override this body via
 // last-registered-wins.  See `%casadi_director_class` /
 // `%casadi_director_vec` below.
%typemap(directorin, noblock=1, fragment="casadi_all") xType, const xType {
    if(!($input = casadi::from_ref($1))) %dirout_fail(SWIG_TypeError,"For director inputs, failed to convert input to " xName ".");
 }

%typemap(directorin, noblock=1, fragment="casadi_all") const xType& {
    if(!($input = casadi::from_ptr(&$1))) %dirout_fail(SWIG_TypeError,"For director inputs, failed to convert input to " xName ".");
 }

 // Enable dynamic dispatch
%typemap(typecheck, noblock=1, fragment="casadi_all") xType &OUTPUT {
  $1 = casadi::to_ptr($input, static_cast< xType **>(0));
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
  %append_output(casadi::from_ptr($1));
 }

// ... but kept as inputs
%typemap(in, doc=xName, pystub_in=xStubIn, tsstub_in=xTsStub, noblock=1, fragment="casadi_all") xType &INOUT (xType m) {
  $1 = &m;
  if (!casadi::to_ptr($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input to type '" xName "'.");
 }

 // ... also for dynamic dispatch
%typemap(typecheck, noblock=1, fragment="casadi_all") xType& INOUT {
  $1 = casadi::to_ptr($input, static_cast< xType **>(0));
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

 // All-in-one template instantiation + typemaps.
 // xTsStub is a TypeScript type literal pasted verbatim into the .d.ts;
 // single arg used for both input and output positions.
%define %casadi_template(xName, xStubIn, xStubOut, xTsStub, xPrec, xType...)
%template() xType;
%casadi_input_typemaps(xName, xStubIn, xTsStub, xPrec, xType)
%casadi_output_typemaps(xName, xStubIn, xStubOut, xTsStub, %arg(xType))
%enddef

// Variant for std::vector<E> -- adds a NAMED template instantiation
// (xJSName) plus the JS-side array<->vector autoconversion typemaps
// (%wasm_vec) when SWIGWASMJS is defined.  Python / MATLAB get the
// usual %casadi_template behaviour (anonymous typecheck registration);
// wasm-js additionally needs the named class because JS callers
// occasionally construct the proxy directly and because the array
// marshaling typemaps reference it by name.
//
// Use at every `%casadi_template(..., std::vector<E>)` site instead of
// declaring a parallel `%template(<Foo>Vector) std::vector<E>` block
// elsewhere.
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
// ----------------------------------------------------------------------------
// Director (C++ <-> JS) typemaps for casadi value types.
//
// The wasm-js target produces wrapped JS proxy instances (e.g. real `DM`
// objects with `.nonzeros()`, `.size1()`) rather than raw `{_ptr: N}`
// carriers, so that a user's `eval(args)` override receives a useful API.
// Achieved by calling `M.__wrap_<JsName>(ptr)` from C++ via embind --
// each JS proxy class registers a wrap helper in f_js_classes:
//   M.__wrap_DM = (p) => new DM(__PRIVATE_CTOR, p);
//
// %casadi_director_class(xJsName, xType...)
//   directorin / directorout for the value type and (in-direction)
//   const reference to it.
//
// %casadi_director_vec(xElemJsName, xElemType...)
//   directorin / directorout for std::vector<xElemType>.  Builds a JS
//   array of wrapped element proxies; reverse drains it back into a C++
//   vector via SWIG_WASMJS_ConvertPtr.
//
// %casadi_director_vec_str
//   Specialisation for std::vector<std::string> -- a JS array of
//   primitive strings, no per-element wrap helper needed.
//
// All bodies follow the wire contract documented in wasm_js.swg:
//   directorin  : $input <- EM_VAL produced by .release_ownership()
//   directorout : $result <- C++ value extracted from $input (EM_VAL)
// ----------------------------------------------------------------------------

%define %casadi_director_class(xJsName, xType...)
// SWIG_NewPointerObj reads $descriptor's clientdata (set in
// init_cast_chains from swig_js_classes[]) and routes to
// M.__wrap_<JsName>(ptr) automatically -- no need for an explicit
// module_property call here.  xJsName is unused but kept for
// callsite-arity compatibility.
%typemap(directorin, noblock=1) xType, const xType {
  {
    xType *__p = new xType($1);
    $input = SWIG_NewPointerObj(__p, $descriptor(xType *), 0);
  }
}
%typemap(directorin, noblock=1) const xType& {
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

%define %casadi_director_vec(xElemJsName, xElemType...)
%typemap(directorin, noblock=1) std::vector< xElemType >,
                                 const std::vector< xElemType >&,
                                       std::vector< xElemType >& {
  // SWIG_NewPointerObj routes through clientdata to the JS wrap
  // automatically (see %casadi_director_class above).  xElemJsName
  // is unused here, kept for callsite-arity compatibility.
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

// std::vector<std::string>: each element is a JS primitive string, not
// a proxy class -- so no __wrap_ helper, plain val::as<std::string>().
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

// ----------------------------------------------------------------------------
// jsout for std::map<std::string, X>: drain the wasm-returned EM_VAL handle
// into a native JS object via `M.__swig_release_handle($call)`.  Without
// this, the regular return path (e.g. `Function.stats()`) leaks the raw
// EM_VAL handle (an i32) to the user.  Pairs with the wasm-js branch of
// `from_ptr<std::map<std::string, M>>` (casadi.i above) which builds the
// JS object on the C++ side.
//
// For dict types whose VALUE type is a casadi value class (DM/MX/SX/...),
// from_ptr<map> recurses into from_ptr<value_type>, which for class types
// produces a `{_ptr: <int>}` carrier object.  So dicts of casadi-class
// values render JS-side as `{key: {_ptr: N}, ...}` -- the user can post-
// wrap with `new <C>(__PRIVATE_CTOR, val._ptr)` if they need the proxy
// methods.  Same trade-off as the (intentional) wasm-js return convention
// for `std::vector<casadi::T>` -> XVector wrapping.
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
%define PREC_PAIR_SLICE_SLICE 93 %enddef
%define PREC_SLICE 94 %enddef
%define PREC_PAIR_IVector_IVector 96 %enddef
%define PREC_IM 97 %enddef
%define PREC_DVector 99 %enddef
%define PREC_DM 100 %enddef
%define PREC_DMVector 101 %enddef
%define PREC_DMVectorVector 101 %enddef
%define PREC_SX 103 %enddef
%define PREC_SXVector 103 %enddef
%define PREC_SXVectorVector 103 %enddef
%define PREC_MX 104 %enddef
%define PREC_MXVector 105 %enddef
%define PREC_MXVectorVector 106 %enddef
%define PREC_CREATOR 150 %enddef
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

#define L_INT "int"
#define L_BOOL "bool"
#define LPAIR(A,B) "(" A "," B ")"

#if defined(SWIGMATLAB) || defined(SWIGOCTAVE)
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
// wasm_js has direct ctype/in/out typemaps for std::string in
// Lib/wasm_js/wasm_js.swg (const char* across the wasm boundary,
// malloc+memcpy for the return path).  The casadi_typemaps version
// would route through casadi::from_ref(string) returning GUESTOBJECT*,
// which doesn't match our typed ABI.
%casadi_typemaps(L_STR, "str", "str", "string", PREC_STRING, std::string)
#endif
%casadi_template_vec(LL L_STR LR, "Sequence[str]", "list[str]", "string[]", PREC_VECTOR, StringVector, std::string)
%casadi_template(LL LL L_STR LR LR, "Sequence[Sequence[str]]", "list[list[str]]", "string[][]", PREC_VECTOR, std::vector<std::vector<std::string> >)

/* -------- PEP-484 stub aliases and preamble.  Active with -stubs. --
 * Each %stub_alias_in(NAME, TYPES) emits `_NAME = TYPES` into the
 * preamble and registers NAME -> _NAME for pystub token substitution.
 * Input aliases widen (to_ptr<T> accepts more than the C++ signature);
 * output aliases stay narrow unless the wrapper unwraps at from_ptr.
 * Do not cross-chain matrix aliases -- narrowest-first overloads want
 * DM subset-of SX and DM subset-of MX but SX and MX must stay independent.
 * ----------------------------------------------------------------- */

#ifdef SWIG_STUBS_ENABLED

%insert("stubs_preamble") %{
import numpy as np
from numpy.typing import NDArray
from collections.abc import Callable
from typing import Protocol, runtime_checkable
_T = TypeVar("_T", "DM", "SX", "MX")

@runtime_checkable
class _SupportsDM(Protocol):
    """Any object exposing a ``__DM__()`` method -- the conversion hook
    used by ``to_ptr<DM>()`` for user-defined numeric types."""
    def __DM__(self) -> "DM": ...

@runtime_checkable
class _SupportsSX(Protocol):
    def __SX__(self) -> "SX": ...

@runtime_checkable
class _SupportsMX(Protocol):
    def __MX__(self) -> "MX": ...

%}

/* issue #2959: casadi.ArrayInterface (a %pythoncode class, invisible to
 * SWIG's -stubs).  Hand-written here via %stubcode (= %insert("stubs")),
 * the body section SWIG also scans into __all__ -- unlike stubs_preamble,
 * so `from casadi import *` re-exports it.  The three backing-typed
 * subclasses each carry exactly one conversion hook, so an SX-backed array
 * satisfies _SupportsSX only (accepted where _SX is) while a numeric
 * DM-backed array satisfies _SupportsDM (hence _DM, and via the unions _SX
 * and _MX too) -- matching casadi's runtime to_ptr rules. */
#ifdef SWIG_STUBS_ENABLED
%stubcode %{class ArrayInterface(Generic[_T]):
    """numpy-semantics array view over a casadi DM/SX/MX.  Constructing one
    dispatches to the backing-typed subclass ArrayInterfaceDM/SX/MX."""
    __array_priority__: float
    @overload
    def __new__(cls, value: "SX", ndim: int = ...) -> "ArrayInterfaceSX": ...
    @overload
    def __new__(cls, value: "MX", ndim: int = ...) -> "ArrayInterfaceMX": ...
    @overload
    def __new__(cls, value: object = ..., ndim: int = ...) -> "ArrayInterfaceDM": ...
    @property
    def shape(self) -> tuple[int, ...]: ...
    @property
    def ndim(self) -> int: ...
    @property
    def size(self) -> int: ...
    def to_casadi(self) -> _T: ...
    @property
    def T(self) -> "ArrayInterface[_T]": ...
    def reshape(self, *shape: int) -> "ArrayInterface[_T]": ...
    def transpose(self, *axes: int) -> "ArrayInterface[_T]": ...
    def squeeze(self, axis: "int | tuple[int, ...] | None" = ...) -> "ArrayInterface[_T]": ...
    def flatten(self) -> "ArrayInterface[_T]": ...
    def ravel(self) -> "ArrayInterface[_T]": ...
    @property
    def flat(self) -> Any: ...
    def to_DM(self) -> "DM": ...
    def sum(self, axis: "int | tuple[int, ...] | None" = ...) -> "ArrayInterface[_T]": ...
    def mean(self, axis: "int | tuple[int, ...] | None" = ...) -> "ArrayInterface[_T]": ...
    def __getitem__(self, idx: Any) -> "ArrayInterface[_T]": ...
    def __setitem__(self, idx: Any, value: Any) -> None: ...
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator["ArrayInterface[_T]"]: ...
    def __array__(self, dtype: Any = ...) -> NDArray[Any]: ...
    def __float__(self) -> float: ...
    def __int__(self) -> int: ...
    def __add__(self, o: Any) -> "ArrayInterface[_T]": ...
    def __radd__(self, o: Any) -> "ArrayInterface[_T]": ...
    def __sub__(self, o: Any) -> "ArrayInterface[_T]": ...
    def __rsub__(self, o: Any) -> "ArrayInterface[_T]": ...
    def __mul__(self, o: Any) -> "ArrayInterface[_T]": ...
    def __rmul__(self, o: Any) -> "ArrayInterface[_T]": ...
    def __truediv__(self, o: Any) -> "ArrayInterface[_T]": ...
    def __pow__(self, o: Any) -> "ArrayInterface[_T]": ...
    def __matmul__(self, o: Any) -> "ArrayInterface[_T]": ...
    def __neg__(self) -> "ArrayInterface[_T]": ...
    @classmethod
    def DM(cls, value: Any, shape: Any = ...) -> "ArrayInterfaceDM": ...
    @classmethod
    def SX(cls, value: Any, shape: Any = ...) -> "ArrayInterfaceSX": ...
    @classmethod
    def MX(cls, value: Any, shape: Any = ...) -> "ArrayInterfaceMX": ...

class ArrayInterfaceDM(ArrayInterface["DM"]):
    def __DM__(self) -> "DM": ...

class ArrayInterfaceSX(ArrayInterface["SX"]):
    def __SX__(self) -> "SX": ...

class ArrayInterfaceMX(ArrayInterface["MX"]):
    def __MX__(self) -> "MX": ...

%}
#endif

/* Matrix input aliases.  SX/MX both extend DM; they don't extend each
 * other -- to_ptr<SX>(MX) and to_ptr<MX>(SX) fail at runtime.
 *
 * NDArray dtype is intentionally `Any`: numpy stubs evolve rapidly
 * (dtype changes between 1.x / 2.x / 2.3+ using `_AnyShape` vs
 * `tuple[int,...]`; literals like `np.linspace(...)` resolve to
 * different concrete dtypes across versions).  Since casadi's C++
 * runtime accepts any numeric ndarray, false-precision dtype narrowing
 * just caused overload-resolution failures that differed between
 * Python/numpy-version combinations.  `NDArray[Any]` is
 * version-agnostic. */
%stub_alias_in(DM, bool | int | float | DM | Sparsity | _SupportsDM | Sequence[bool | int | float] | Sequence[Sequence[bool | int | float]] | NDArray[Any])
%stub_alias_in(SX, SX | _SupportsSX | _DM)
%stub_alias_in(MX, MX | _SupportsMX | _DM)
%stub_alias_in(IM, IM | int | Sequence[int] | Sequence[Sequence[int]] | NDArray[Any])
%stub_alias_in(SXElem, SXElem | float)

/* Bracketed doc tokens -- not valid SWIG identifiers. */
%stub_alias_in_key([float], Sequence[bool | int | float] | NDArray[Any])
%stub_alias_in_key([int],   Sequence[bool | int] | NDArray[Any])
%stub_alias_in_key([bool],  Sequence[bool | int] | NDArray[Any])

%stub_alias_in(Slice, Slice | int | slice)

/* Input: precise recursive union -- matches what to_ptr<GenericType>()
 * accepts, guides users constructing options dicts like
 * `{"ipopt": {"tol": 1e-6, ...}, "print_time": True}`.
 *
 * Output: `Any`, deliberately.  Reasons:
 * (1) On return paths (`stats()`, `info()`, etc.) users compare against
 *     plain scalars (`stats["nsteps"] < 2500`).  A recursive union
 *     would disallow operators because not all members support them.
 * (2) pyright's cross-module resolution of our recursive alias is
 *     typeshed/target-dependent: on Py3.11 it silently falls back to
 *     `Unknown` (hiding the real errors), on Py3.8 it resolves and
 *     flags them.  `Any` is stable across all targets.
 * (3) Users inspect stats output freely; precise typing buys nothing. */
%stub_alias_in(GenericType, bool | int | float | str | "Function" | Sequence["_GenericType"] | %arg(Mapping[str, "_GenericType"]))
%stub_alias_out_key(GenericType, Any)

/* __getitem__/__setitem__ axis index.  MX additionally accepts MX
 * indices, added via a narrower %extend overload on MX itself. */
%stub_alias_in(MIndex, int | slice | Sequence[bool | int] | NDArray[Any] | Sparsity | DM)

/* CasadiMatrix operator overloads.  Emitted on GenericExpressionCommon,
 * which DM/SX/MX all inherit -- per-class overloads narrow the `other`
 * argument to what each class actually accepts at runtime.
 *
 *   DM.__add__(other: _DM) -> DM
 *   SX.__add__(other: _SX) -> SX        (_SX = SX | _SupportsSX | _DM)
 *   MX.__add__(other: _MX) -> MX        (_MX = MX | _SupportsMX | _DM)
 *
 * Crucially, MX + SX is NOT allowed: SX isn't in _MX, and MX isn't in
 * _SX.  Both `__add__` on MX and `__radd__` on SX fail to match,
 * mirroring the runtime behaviour. */
%define %stub_CasadiMatrix_binop(NAME)
%stub_overload_method_selftyped(NAME, DM, DM, other: _DM)
%stub_overload_method_selftyped(NAME, SX, SX, other: _SX)
%stub_overload_method_selftyped(NAME, MX, MX, other: _MX)
%enddef

%define %stub_CasadiMatrix_unop(NAME)
%stub_method0(NAME, Self)
%enddef

%define %stub_CasadiMatrix_cmp(NAME)
%stub_overload_method_selftyped(NAME, DM, DM, other: _DM)
%stub_overload_method_selftyped(NAME, SX, SX, other: _SX)
%stub_overload_method_selftyped(NAME, MX, MX, other: _MX)
%enddef

#else

%define %stub_CasadiMatrix_binop(NAME) %enddef
%define %stub_CasadiMatrix_unop(NAME)  %enddef
%define %stub_CasadiMatrix_cmp(NAME)   %enddef

#endif

%casadi_typemaps("Sparsity", "Sparsity", "Sparsity", "Sparsity", PREC_SPARSITY, casadi::Sparsity)
%casadi_template_vec(LL "Sparsity" LR, "Sequence[Sparsity]", "list[Sparsity]", "Sparsity[]", PREC_SPARSITY, SparsityVector, casadi::Sparsity)
%casadi_template(LL LL "Sparsity"  LR  LR, "Sequence[Sequence[Sparsity]]", "list[list[Sparsity]]", "Sparsity[][]", PREC_SPARSITY, std::vector<std::vector< casadi::Sparsity> >)
%casadi_template(LDICT("Sparsity"), "Mapping[str, Sparsity]", "dict[str, Sparsity]", "Record<string, Sparsity>", PREC_SPARSITY, std::map<std::string, casadi::Sparsity >)
%casadi_template(LDICT(LL "Sparsity" LR), "Mapping[str, Sequence[Sparsity]]", "dict[str, list[Sparsity]]", "Record<string, Sparsity[]>", PREC_SPARSITY, std::map<std::string, std::vector<casadi::Sparsity > >)
%casadi_template(LPAIR(LDICT("Sparsity"),"[" L_STR "]"), "tuple[Mapping[str, Sparsity], Sequence[str]]", "tuple[dict[str, Sparsity], list[str]]", "[Record<string, Sparsity>, string[]]", PREC_SPARSITY, std::pair<std::map<std::string, casadi::Sparsity >, std::vector<std::string> >)
#ifndef SWIGWASMJS
// wasm_js: bool / casadi_int handled by direct ctype/in/out typemaps
// in Lib/wasm_js/wasm_js.swg; vector instantiations go through
// %template(...) + %wasm_vec(...) for JS-side Array marshaling.
%casadi_typemaps(L_BOOL, "bool", "bool", "boolean", SWIG_TYPECHECK_BOOL, bool)
%casadi_template("[" L_BOOL "]", "Sequence[bool | int] | NDArray[Any]", "list[bool]", "boolean[]", SWIG_TYPECHECK_BOOL, std::vector<bool>)
%casadi_template("[[" L_BOOL "]]", "Sequence[Sequence[bool | int]] | Sequence[NDArray[Any]]", "list[list[bool]]", "boolean[][]", SWIG_TYPECHECK_BOOL, std::vector<std::vector<bool> >)
%casadi_typemaps( L_INT , "int", "int", "bigint", SWIG_TYPECHECK_INTEGER, casadi_int)
#endif

#ifdef MATLABSTYLE
#define LABEL "[int,int]"
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
%casadi_typemaps("SXElem", "SXElem", "SXElem", "SXElem", PREC_SX, casadi::SXElem)
%casadi_template(LL "SXElem" LR, "Sequence[SXElem]", "list[SXElem]", "SXElem[]", PREC_SXVector, std::vector<casadi::SXElem>)
%casadi_typemaps("SX", "SX", "SX", "SX", PREC_SX, casadi::Matrix<casadi::SXElem>)
%casadi_template_vec(LL "SX" LR, "Sequence[SX]", "list[SX]", "SX[]", PREC_SXVector, SXVector, casadi::Matrix<casadi::SXElem>)
%casadi_template(LL LL "SX" LR LR, "Sequence[Sequence[SX]]", "list[list[SX]]", "SX[][]", PREC_SXVectorVector, std::vector<std::vector< casadi::Matrix<casadi::SXElem> > >)
%casadi_template(LDICT("SX"), "Mapping[str, SX]", "dict[str, SX]", "Record<string, SX>", PREC_SX, std::map<std::string, casadi::Matrix<casadi::SXElem> >)
%casadi_typemaps("MX", "MX", "MX", "MX", PREC_MX, casadi::MX)
%casadi_template_vec(LL "MX" LR, "Sequence[MX]", "list[MX]", "MX[]", PREC_MXVector, MXVector, casadi::MX)
%casadi_template(LL LL "MX" LR LR, "Sequence[Sequence[MX]]", "list[list[MX]]", "MX[][]", PREC_MXVectorVector, std::vector<std::vector<casadi::MX> >)
%casadi_template(LDICT("MX"), "Mapping[str, MX]", "dict[str, MX]", "Record<string, MX>", PREC_MX, std::map<std::string, casadi::MX>)
%casadi_template(LPAIR("MX","MX"), "tuple[MX, MX]", "tuple[MX, MX]", "[MX, MX]", PREC_MXVector, std::pair<casadi::MX, casadi::MX>)
%casadi_typemaps("DM", "DM", "DM", "DM", PREC_DM, casadi::Matrix<double>)
%casadi_template_vec(LL "DM" LR, "Sequence[DM]", "list[DM]", "DM[]", PREC_DMVector, DMVector, casadi::Matrix<double>)
%casadi_template(LL LL "DM" LR LR, "Sequence[Sequence[DM]]", "list[list[DM]]", "DM[][]", PREC_DMVectorVector, std::vector<std::vector< casadi::Matrix<double> > >)
%casadi_template(LDICT("DM"), "Mapping[str, DM]", "dict[str, DM]", "Record<string, DM>", PREC_DM, std::map<std::string, casadi::Matrix<double> >)
%casadi_typemaps("IM", "IM", "IM", "IM", PREC_IM, casadi::Matrix<casadi_int>)
// Without CASADI_INT_TYPE, you get SwigValueWrapper
// With it, docstrings are screwed
%casadi_typemaps("GenericType", "GenericType", "GenericType", "GenericType", PREC_GENERICTYPE, casadi::GenericType)
%casadi_template(LL "GenericType" LR, "Sequence[GenericType]", "list[GenericType]", "GenericType[]", PREC_GENERICTYPE, std::vector<casadi::GenericType>)
%casadi_template(LL LL "GenericType" LR LR, "Sequence[Sequence[GenericType]]", "list[list[GenericType]]", "GenericType[][]", PREC_GENERICTYPE, std::vector<std::vector<casadi::GenericType> >)
%casadi_typemaps("Slice", "Slice", "Slice", "Slice", PREC_SLICE, casadi::Slice)
%casadi_typemaps("Function", "Function", "Function", "Function", PREC_FUNCTION, casadi::Function)
%casadi_template_vec(LL "Function" LR, "Sequence[Function]", "list[Function]", "Function[]", PREC_FUNCTION, FunctionVector, casadi::Function)
%casadi_template(LPAIR("Function","Function"), "tuple[Function, Function]", "tuple[Function, Function]", "[Function, Function]", PREC_FUNCTION, std::pair<casadi::Function, casadi::Function>)
%casadi_template(L_DICT, "Mapping[str, GenericType]", "dict[str, GenericType]", "Record<string, GenericType>", PREC_DICT, std::map<std::string, casadi::GenericType>)
%casadi_template(LDICT(LL L_STR LR), "Mapping[str, Sequence[str]]", "dict[str, list[str]]", "Record<string, string[]>", PREC_DICT, std::map<std::string, std::vector<std::string> >)

#ifdef SWIGWASMJS
// ----------------------------------------------------------------------------
// Director (C++ <-> JS) typemap registrations for the wasm-js target.
//
// These OVERRIDE the unconditional %typemap(directorin)/(directorout)
// emitted by %casadi_input_typemaps / %casadi_output_typemaps above (which
// use casadi::from_ref / casadi::to_val and produce raw `{_ptr}` carriers
// on the JS side).  For wasm-js we want the JS subclass to see real
// proxy class instances (with `.size1()`, `.nonzeros()`, etc.) so we
// emit per-class typemaps that wrap via M.__wrap_<JsName>.
//
// Ordering matters: must appear AFTER the %casadi_typemaps invocations
// so that the wasm-js typemap is the LAST-registered for each pattern
// and wins typemap lookup.
// ----------------------------------------------------------------------------

%casadi_director_class("Sparsity",    casadi::Sparsity)
%casadi_director_class("DM",          casadi::Matrix<double>)
%casadi_director_class("MX",          casadi::MX)
%casadi_director_class("SX",          casadi::Matrix<casadi::SXElem>)
%casadi_director_class("IM",          casadi::Matrix<casadi_int>)
%casadi_director_class("Function",    casadi::Function)
%casadi_director_class("Slice",       casadi::Slice)
%casadi_director_class("SXElem",      casadi::SXElem)
// GenericType has no `M.__wrap_GenericType` proxy on the JS side
// (callers see plain JS values via `casadi::to_ptr`-style coercion),
// so it's intentionally NOT registered as a director class here.

%casadi_director_vec("Sparsity",      casadi::Sparsity)
%casadi_director_vec("DM",            casadi::Matrix<double>)
%casadi_director_vec("MX",            casadi::MX)
%casadi_director_vec("SX",            casadi::Matrix<casadi::SXElem>)
%casadi_director_vec("Function",      casadi::Function)

%casadi_director_vec_str

// jsout for every dict-shaped return type the casadi API surfaces.
// Drains the wasm-returned EM_VAL handle into a native JS object.
// Affects Function.stats() (returns Dict), all opts:Dict accessors,
// the result.x_dict etc. fields of solver returns, and so on.
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

// Named std::vector instantiations + array<->vector autoconversion
// typemaps are now folded into the relevant %casadi_template_vec
// calls above (one site per type instead of two parallel blocks).

#undef L_INT
#undef L_BOOL
#undef LPAIR
#undef L_DOUBLE
#undef L_DICT
#undef LL
#undef LR
#undef L_STR
#undef MATLABSTYLE

/* TypeScript .d.ts preamble: aliases for casadi-side typedefs that
   appear in signatures by name (typedef alias kept by SWIG, not
   resolved to the underlying type before stub emission).  Without
   these the generated .d.ts has dangling references like `Dict`,
   `MXDict`, etc.  wasm_js.cxx is intentionally casadi-agnostic --
   these mappings live here, in casadi.i, as TS-equivalent typedefs. */
#ifdef SWIGWASMJS
%insert("stubs") %{
// --- casadi-specific TypeScript type aliases ---
export type Dict         = Record<string, GenericType>;
export type DMDict       = Record<string, DM>;
export type MXDict       = Record<string, MX>;
export type SXDict       = Record<string, SX>;
export type SparsityDict = Record<string, Sparsity>;
export type DMVector     = DM[];
export type MXVector     = MX[];
export type SXVector     = SX[];
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

%}

// ---------------------------------------------------------------------------
// JS-side ergonomics: post-class-definition overrides.  The "js" slot is
// emitted by wasm_js.cxx after all class definitions and before the module
// `return { ... }`, so the proxy classes (DM, DMVector, Function, ...) are
// in lexical scope here.
// ---------------------------------------------------------------------------
%insert("js") %{
  /* Function.call vs callable: two distinct conventions, mirroring
     Python's casadi.
       f.call([a, b])          -- ALWAYS returns the full list of outputs
                                  (`[y0, y1, ...]`).  Programmatically
                                  stable.  Element type (DM/MX/SX) is
                                  detected from the input list and the
                                  matching overload is invoked, so the
                                  return is `DM[]`, `MX[]`, or `SX[]`
                                  respectively (mirrors Python's
                                  type-overloaded f.call).
       f.call({name: value})   -- dict input -> dict output keyed by
                                  the function's output names; value
                                  type detected from the input values
                                  (DMDict/MXDict/SXDict overloads).
       f(a, b)                 -- syntax sugar (variadic).  Return shape
                                  depends on n_out:
                                     n_out == 0  -> null
                                     n_out == 1  -> bare output
                                     n_out  > 1  -> list
                                  Matches Python's `f(a, b)`.
       f({i0: x})              -- dict-input sugar; returns dict.

     Implementation note: the SWIG-emit-time Function.prototype.call
     dispatcher had the "first-wins-per-arity" trap -- only the DM
     overload was generated for arity 3, so `f.call([MX])` errored
     with `Failed to convert input 1 to type 'DM'`.  We REPLACE the
     dispatcher here with one that probes the input type and routes
     to one of the 6 wasm exports
       _swig_casadi_Function_call_0  vector<DM>   -> DMVector
       _swig_casadi_Function_call_1  vector<SX>   -> SXVector
       _swig_casadi_Function_call_2  vector<MX>   -> MXVector
       _swig_casadi_Function_call_3  DMDict       -> DMDict
       _swig_casadi_Function_call_4  SXDict       -> SXDict
       _swig_casadi_Function_call_5  MXDict       -> MXDict
     all of which the SWIG-emit already generates on the C++ side. */
  Function.prototype.call = function(arg, always_inline, never_inline) {
    if (always_inline === undefined) always_inline = false;
    if (never_inline  === undefined) never_inline  = false;
    const ail = always_inline ? 1 : 0;
    const nil = never_inline  ? 1 : 0;
    // Detect positional (list-shaped) vs dict input.
    const is_list = Array.isArray(arg) || (arg instanceof DMVector)
                 || (arg instanceof MXVector) || (arg instanceof SXVector);
    const is_dict = !is_list && arg && typeof arg === 'object'
                 && !(arg instanceof DM) && !(arg instanceof MX)
                 && !(arg instanceof SX) && !(arg instanceof Sparsity);
    if (is_list) {
      // Detect element type to pick the overload.
      // Default: DM (matches the existing wasm-js convention for raw
      // arrays of numbers via casadi's array-to-DM auto-coercion).
      let mode = 'DM';
      if (arg instanceof MXVector) mode = 'MX';
      else if (arg instanceof SXVector) mode = 'SX';
      else if (Array.isArray(arg)) {
        for (const e of arg) {
          if (e instanceof MX) { mode = 'MX'; break; }
          if (e instanceof SX) { mode = 'SX'; break; }
        }
      }
      let vec, own = false;
      if (Array.isArray(arg)) {
        const VecCls = mode === 'MX' ? MXVector : (mode === 'SX' ? SXVector : DMVector);
        vec = new VecCls();
        own = true;
        for (const x of arg) vec.push_back(x);
      } else {
        vec = arg;
      }
      const fn = mode === 'DM' ? M._swig_casadi_Function_call_0 :
                 mode === 'SX' ? M._swig_casadi_Function_call_1 :
                                M._swig_casadi_Function_call_2;
      const OutVec = mode === 'DM' ? DMVector : (mode === 'SX' ? SXVector : MXVector);
      // Vector returns: from_ptr<std::vector<M>> in casadi.i passes
      // type=0 to SWIG_NewPointerObj, so __from_handle yields a bare
      // {_ptr} carrier rather than a fully-wrapped XVector.  Reach
      // into ._ptr to recover the raw pointer; both bare carriers and
      // typed proxies expose it, so this is safe either way.
      const result = __vec_to_arr(new OutVec(__PRIVATE_CTOR,
        __from_handle(__chk(fn(__unwrap(this), __unwrap(vec), ail, nil)))._ptr));
      if (own && vec.delete) vec.delete();
      return result;
    }
    if (is_dict) {
      // Detect value type from the dict's actual values.
      let mode = 'DM';
      for (const k of Object.keys(arg)) {
        const v = arg[k];
        if (v instanceof MX) { mode = 'MX'; break; }
        if (v instanceof SX) { mode = 'SX'; break; }
      }
      // Route via the DICT wasm exports (_call_3 = DMDict, _call_4 =
      // SXDict, _call_5 = MXDict).  Those honor casadi's "missing key
      // = use default" convention -- nlpsol's lbx/ubx default to
      // -inf/+inf when the key is absent, etc.  Going via the LIST
      // path would force every slot to receive an empty-DM, which the
      // solver interprets as a length-0 user-provided bound (defaults
      // x to 0).  Dict values come back as fully-wrapped proxies via
      // SWIG_WASMJS_NewPointerObj (which consults swig_type_info::
      // clientdata and calls M.__wrap_<JsName> on the C++ side), so
      // no JS-side rewrap is needed.
      const fn  = mode === 'DM' ? M._swig_casadi_Function_call_3
               : mode === 'SX' ? M._swig_casadi_Function_call_4
               :                  M._swig_casadi_Function_call_5;
      return M.__swig_release_handle(__chk(fn(__unwrap(this),
        __unwrap(arg), ail, nil)));
    }
    throw new Error('Function.call: arg must be a list (DM/MX/SX vector or JS array) or a dict object');
  };

  /* x.T as a getter property -- matches Python's `x.T`.  The original
     auto-generated method stays available as x.T_method() should
     anyone hold a stale binding; the getter is the canonical form. */
  for (const Cls of [DM, MX, SX, Sparsity]) {
    const __orig_T = Cls.prototype.T;
    if (!__orig_T) continue;
    Object.defineProperty(Cls.prototype, "T", {
      get() { return __orig_T.call(this); },
      configurable: true,
    });
  }

  /* Variadic + list aliases for the concat family, mirroring Python:
       horzcat(a, b, c, d)   variadic           -- 1xN row (stack)
       hcat([a, b, c, d])    list form          -- 1xN row (stack)
       horzcat([1, 2, 3, 4]) 1-arg list of nums -- 4x1 column (identity)
       (analogous for vertcat/vcat, diagcat/dcat).
     Python's `def horzcat(*args): return _horzcat(args)` always passes
     the variadic rest-tuple as the single C++ "vector input" arg.  We
     mirror that exactly: always pass the rest-array.  When the rest
     has 1 element AND that element is a JS array of primitives, the
     C++ to_ptr<DM> nested-array handler builds a column DM, and the
     C++ horzcat is an identity on the single column.  When the rest
     has 1 element that is a list of matrices (e.g. `vertcat([x, y])`),
     it FAILS in Python too -- callers should use `vertcat(x, y)`
     variadic or `vcat([x, y])` list-form instead. */
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

  /* Factory + callable Function instances.  Two ergonomic wins together:
       1. CONSTRUCTION:
            const f = M.Function('name', [x], [y]);   // factory (preferred)
            const f = new M.Function('name', [x], [y]); // also works
          Established convention in test/javascript/*.js is the factory
          form (no `new`).  Internally `Function(...)` is just a regular
          function that delegates to `new __OrigFunction(...)`; since the
          body returns an object (the callable), `new Function(...)` and
          `Function(...)` produce the same callable per JS semantics.
       2. INVOCATION:
            f([x, y])         <- direct call (sugar for f.call([x, y]))
            f({i0: x})        <- dict-style call (sugar for f.call({...}))
            f.call([x, y])    <- explicit
       Mirrors Python's `f([x, y])` and matches the natural JS calling
       convention.

     Implementation:  the constructor-returns-a-callable pattern.  Build
     a real JS function whose `[[Prototype]]` is our `Function.prototype`,
     so all methods (call, n_in, delete, name_in, expand, ...) still work
     via prototype lookup AND the function-invocation slot delegates to
     `.call`.  See the inner-binding-immutability workaround below for
     why every prototype method needs to be re-wrapped at registration.

     Lifecycle: the wasm pointer was originally registered with
     `__fr_casadi_Function` against the bare inst (held internally by the
     class constructor).  We transfer registration to the callable so its
     destruction (rather than the orphan inst's) triggers
     `_swig_casadi_Function_delete`.  Setting `inst._ptr = 0` prevents
     inst.delete() from running on the now-stale handle.

     Why not extend Function: `class Wrapped extends Function { ... }`
     would require `super(<body>)` with a string body that can't capture
     enclosing scope -- forcing global state plumbing.  The setPrototypeOf
     approach is cleaner and equally well-supported.

     SWIG-emitted internal call-sites that build Function instances via
     `new Function(__PRIVATE_CTOR, ptr)` (e.g. Function.find_function,
     f.reverse()) bypass the OUTER `Function` binding because ES6 class
     method bodies bind to the INNER (immutable) class binding -- so we
     ALSO patch every method on Function.prototype to post-wrap any
     bare Function instance in its return value.  See the
     `__wrapFunctionReturning` block below. */
  const __OrigFunction = Function;
  // Function methods that collide with own properties of every JS function
  // instance (`name`, `length`, `caller`, `arguments`, `prototype`).  Of
  // the 103 methods on casadi::Function only `name()` actually collides.
  // The collision matters because the JS engine sets `callable.name` as an
  // own property at function-creation time -- the prototype-chain lookup
  // that would otherwise resolve to our SWIG method is shadowed.
  // Re-define the own property to point at the SWIG method so `f.name()`
  // still returns the function's casadi-name.  (`f.name` thus returns the
  // method function rather than the JS-engine's debug name -- a minor
  // trade-off; consoles will display "[Function: name]" instead of "".)
  const __builtinFnOwn = new Set(Object.getOwnPropertyNames(function(){}));
  const __FunctionMethodCollisions = Object.getOwnPropertyNames(
      __OrigFunction.prototype).filter(n =>
        n !== 'constructor' && __builtinFnOwn.has(n));

  // Turn a freshly-constructed plain Function instance into a callable.
  // Transfers ownership of the wasm ptr's FinalizationRegistry registration
  // from `inst` to `callable` and zeros `inst._ptr` so the orphan inst's
  // GC doesn't double-free.
  //
  // Calling convention (mirrors Python's casadi):
  //   f(a, b)        variadic positional, n_out-dependent return shape
  //   f({i0: a})     dict input; returns dict (always)
  //
  // The args passed to `f(...)` ARE the positional inputs, one per
  // C++ parm.  We do NOT flatten a single-array arg into multiple
  // positionals -- `f([1, 2])` for an n_in=1 function passes the
  // ENTIRE `[1, 2]` array as the first input (where casadi's
  // `to_ptr<DM>` auto-coerces it to a column DM, matching Python's
  // `f([1, 2])` semantics for a vector-valued input).  Calling
  // `f([a, b])` against an n_in=2 function errors loudly -- use
  // `f(a, b)` for variadic-positional.
  //
  // For positional input the return shape is:
  //   n_out == 0     null      (no outputs; useful for side-effect-only fns)
  //   n_out == 1     bare      (single output unwrapped)
  //   n_out  > 1     list      (full output array)
  //
  // Programmatic users wanting a stable return shape should use
  // `f.call([...])` -- the patched .call above always returns a list.
  function __makeCallableFunction(inst) {
    const callable = function (...invokeArgs) {
      // dict-input branch: route through .call which returns a dict.
      // The dict-shape check must EXCLUDE casadi matrix types (DM, MX,
      // SX, Sparsity, IM) which are objects too -- otherwise
      // `f(M.DM(5))` would be misinterpreted as a {...}-dict input.
      if (invokeArgs.length === 1) {
        const a = invokeArgs[0];
        if (a && typeof a === 'object' && !Array.isArray(a)
            && !(a instanceof DMVector)
            && !(a instanceof MXVector)
            && !(a instanceof SXVector)
            && !(a instanceof DM)
            && !(a instanceof MX)
            && !(a instanceof SX)
            && !(a instanceof Sparsity)
            && (typeof IM === 'undefined' || !(a instanceof IM))) {
          return __OrigFunction.prototype.call.call(callable, a);
        }
      }
      // positional: invokeArgs IS the positional list (one entry per parm).
      const result = __OrigFunction.prototype.call.call(callable, invokeArgs);
      // After the un-wrap removal from .call, result is always a list.
      const n_out = Number(callable.n_out());
      if (n_out === 0) return null;
      if (n_out === 1) return result[0];
      return result;
    };
    Object.setPrototypeOf(callable, __OrigFunction.prototype);
    callable._ptr = inst._ptr;
    __fr_casadi_Function.unregister(inst);
    __fr_casadi_Function.register(callable, inst._ptr, callable);
    inst._ptr = 0;
    for (const __k of __FunctionMethodCollisions) {
      Object.defineProperty(callable, __k, {
        value: __OrigFunction.prototype[__k],
        writable: false,
        configurable: true,
      });
    }
    return callable;
  }

  function __WrappedFunction(...args) {
    return __makeCallableFunction(new __OrigFunction(...args));
  }
  __WrappedFunction.prototype = __OrigFunction.prototype;
  __WrappedFunction.prototype.constructor = __WrappedFunction;
  // Copy static class members (e.g. Function.expand, Function.if_else,
  // Function.deserialize, ...) onto the wrapper so M.Function.expand(...)
  // keeps working.
  for (const __k of Object.getOwnPropertyNames(__OrigFunction)) {
    if (__k === 'length' || __k === 'name' || __k === 'prototype') continue;
    const __d = Object.getOwnPropertyDescriptor(__OrigFunction, __k);
    if (__d) Object.defineProperty(__WrappedFunction, __k, __d);
  }
  Function = __WrappedFunction;
  __m.Function = __WrappedFunction;

  // Class-method bodies refer to the INNER class binding (immutable per
  // ES6 semantics) -- so `new Function(__PRIVATE_CTOR, ptr)` inside e.g.
  // `Function.prototype.expand` always invokes the ORIGINAL constructor,
  // never `__WrappedFunction`.  As a result `f.expand()` would return a
  // plain non-callable Function.  Workaround: monkey-patch every method
  // on Function.prototype so its return value is post-wrapped if it's a
  // bare Function instance.  Also patch static methods (`Function.expand`,
  // `Function.deserialize`, ...) for the same reason.
  const __wrapFunctionReturning = (orig) => function (...args) {
    const r = orig.apply(this, args);
    if (r && typeof r === 'object' && r instanceof __OrigFunction
            && typeof r !== 'function' && r._ptr !== 0) {
      return __makeCallableFunction(r);
    }
    return r;
  };
  for (const __k of Object.getOwnPropertyNames(__OrigFunction.prototype)) {
    if (__k === 'constructor') continue;
    const __d = Object.getOwnPropertyDescriptor(__OrigFunction.prototype, __k);
    if (!__d || typeof __d.value !== 'function') continue;
    Object.defineProperty(__OrigFunction.prototype, __k,
      { ...__d, value: __wrapFunctionReturning(__d.value) });
  }
  for (const __k of Object.getOwnPropertyNames(__WrappedFunction)) {
    if (__k === 'length' || __k === 'name' || __k === 'prototype') continue;
    const __d = Object.getOwnPropertyDescriptor(__WrappedFunction, __k);
    if (!__d || typeof __d.value !== 'function') continue;
    Object.defineProperty(__WrappedFunction, __k,
      { ...__d, value: __wrapFunctionReturning(__d.value) });
  }
%}
#endif

// Matlab is index-1 based
#ifdef SWIGMATLAB
%typemap(in, doc="index", pystub_in="int", noblock=1) casadi_index {
  if (!casadi::to_val($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input $argnum to type ' index '.");
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

/* `inf` and `pi` exist at runtime (casadi/core const doubles) but are
 * deliberately NOT declared in the stub.  numpy declares both as
 * `Final[float]`; a parallel casadi declaration triggers "declared as
 * Final and cannot be reassigned" on any user file that does both
 * `from casadi import *` and `from numpy import *`.  Users who need
 * these names in a typed context should import them from numpy. */

/* Stubs for the numpy-flavoured aliases above.  Each dispatches to
 * the corresponding casadi.<atan/asin/...> overload; overloads mirror
 * the underlying auto-generated ones so the return type narrows. */
%define %stub_unary_math(NAME)
%stub_overload_func(NAME, DM, x: _DM)
%stub_overload_func(NAME, SX, x: _SX)
%stub_overload_func(NAME, MX, x: _MX)
%enddef
%stub_unary_math(arcsin)
%stub_unary_math(arccos)
%stub_unary_math(arctan)
%stub_unary_math(arctanh)
%stub_unary_math(arcsinh)
%stub_unary_math(arccosh)
%stub_overload_func(arctan2, DM, y: _DM, x: _DM)
%stub_overload_func(arctan2, SX, y: _SX, x: _SX)
%stub_overload_func(arctan2, MX, y: _MX, x: _MX)
#endif // SWIGPYTHON

// Strip leading casadi_ unless followed by ML/int
%rename("%(regex:/casadi_(?!ML|int\\b)(.*)/\\1/)s") "";
%rename(casadi_int) "casadi_int";

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

/* Stubs for the %pythoncode *cat wrappers above.  Overloads narrow
 * the return to the dominant input type.  Order matters: pyright
 * picks the first matching overload, so narrowest (Sparsity) first,
 * then DM, then SX, then MX -- otherwise `vertcat(1, 1)` (all inputs
 * in _DM) resolves to the MX overload because _MX is a superset
 * of _DM. */
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
%feature("python:maybecall") casadi_min;
%feature("python:maybecall") casadi_max;
%feature("python:maybecall") casadi_and;
%feature("python:maybecall") casadi_or;
%feature("python:maybecall") casadi_mod;
%feature("python:maybecall") casadi_copysign;
%feature("python:maybecall") casadi_constpow;
#endif // SWIGPYTHON

#ifdef SWIGMATLAB
%rename(uminus) operator-;
%rename(uplus) operator+;
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

/* NZproxy has no C++ counterpart; it's defined via %pythoncode.
 * Parameterised on the owning matrix type so `x.nz[i]` narrows to
 * the same matrix type as `x`.  MX.nz additionally accepts MX-typed
 * masks -- the _selftyped overload matches only when _T=MX, keeping
 * DM.nz / SX.nz narrow. */
#ifdef SWIG_STUBS_ENABLED
%stubcode %{class NZproxy(Generic[_T]):
%}
%stub_method(__init__,    None, matrix: _T)
%stub_overload_method_selftyped(__getitem__, MX, "NZproxy[MX]", s: _MIndex | MX)
%stub_overload_method(__getitem__, _T,  s: _MIndex)
%stub_overload_method_selftyped(__setitem__, None, "NZproxy[MX]", s: _MIndex | MX, val: bool | int | float | MX | Sequence[bool | int | float])
%stub_overload_method(__setitem__, None, s: _MIndex, val: bool | int | float | _T | Sequence[bool | int | float])
%stub_method0(__len__,  int)
%stub_method0(__iter__, Iterator[_T])
#endif

%define %matrix_helpers(Type)
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
      raise Exception("""CasADi matrices are not iterable by design.
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
/* __getitem__ / __setitem__: _MIndex is the common (DM-compatible)
 * index set shared by DM/SX/MX.  MX additionally accepts MX as index
 * (see the MX %extend block below), registered there via a widening
 * overload so Self stays narrow. */
%stub_overload_method(__getitem__, Self, s: _MIndex | tuple[_MIndex, _MIndex])
%stub_overload_method(__setitem__, None, s: _MIndex | tuple[_MIndex, _MIndex], val: bool | int | float | DM | SX | MX | Sequence[bool | int | float] | Sequence[Sequence[bool | int | float]] | NDArray[Any])
%stub_method0(__iter__, %arg(Iterator[Self]))
%enddef

%define %python_array_wrappers(arraypriority)
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
/* __array__ makes DM / SX / MX acceptable where numpy / scipy expect
 * an ArrayLike.  Runtime path: `.full()` for DM (dense numeric),
 * scalar-object boxing for symbolic.  Keeps scipy.linalg.solve(A,b)
 * etc. typing without forcing users to wrap calls in np.array().  */
%stub_method(__array__, %arg(NDArray[np.float64]), *args: Any, **kwargs: Any)
%enddef
#endif // SWIGPYTHON

#if defined(SWIGXML) || defined(SWIGWASMJS)
// Both metadata-only (XML) and pointer-ABI (wasm_js) targets skip the
// per-language ergonomics layer; the raw class methods are enough.
%define %matrix_helpers(Type)
%enddef
#endif

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

%define %matrix_helpers(Type)
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
    const Type paren(char rr, char cc) const {
      Type m;
      $self->get(m, true, casadi::char2Slice(rr), casadi::char2Slice(cc));
      return m;
    }
    const Type paren(char rr, const Matrix<casadi_int>& cc) const {
      Type m;
      $self->get(m, true, casadi::char2Slice(rr), cc);
      return m;
    }
    const Type paren(const Matrix<casadi_int>& rr, char cc) const {
      Type m;
      $self->get(m, true, rr, casadi::char2Slice(cc));
      return m;
    }
    const Type paren(const Matrix<casadi_int>& rr, const Matrix<casadi_int>& cc) const {
      Type m;
      $self->get(m, true, rr, cc);
      return m;
    }

    // Set a submatrix (index-1)
    void paren_asgn(const Type& m, char rr) {
      casadi_assert_dev(rr==':');
      $self->set(m, false, casadi::IM(casadi::range($self->numel())));
    }
    void paren_asgn(const Type& m, const Matrix<casadi_int>& rr) { $self->set(m, true, rr);}
    void paren_asgn(const Type& m, const Sparsity& sp) { $self->set(m, true, sp);}
    void paren_asgn(const Type& m, char rr, char cc) { $self->set(m, true, casadi::char2Slice(rr), casadi::char2Slice(cc));}
    void paren_asgn(const Type& m, char rr, const Matrix<casadi_int>& cc) { $self->set(m, true, casadi::char2Slice(rr), cc);}
    void paren_asgn(const Type& m, const Matrix<casadi_int>& rr, char cc) { $self->set(m, true, rr, casadi::char2Slice(cc));}
    void paren_asgn(const Type& m, const Matrix<casadi_int>& rr, const Matrix<casadi_int>& cc) { $self->set(m, true, rr, cc);}

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

%enddef
#endif

#if defined(SWIGMATLAB) || defined(SWIGOCTAVE)
// Mark the empty-struct mixin bases as abstract for the MATLAB/Octave
// generator. Effect: they inherit from handle instead of SwigRef, emit no
// delete or swig_this, and the concrete leaf appends "& SwigRef" so the
// diamond collapses to a single SwigRef path. Breaks the diamond that
// trips Octave 10+ (Savannah bug #50011 / #66930); no effect on MATLAB
// correctness but removes the fragile per-class delete chain.
// NOTE: these %feature directives must appear before the corresponding
// %include that brings each class into the SWIG parse.
%feature("abstract") casadi::PrintableCommon;
%feature("abstract") casadi::GenericExpressionCommon;
%feature("abstract") casadi::GenericMatrixCommon;
%feature("abstract") casadi::SparsityInterfaceCommon;
%feature("abstract") casadi::MatrixCommon;
// GenSharedObject is a %template(GenSharedObject) instantiation of
// GenericShared<...> (see the %template directive further down). The
// feature has to be applied to the template base class, not the
// instantiated short name, so it takes effect when the template expands.
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
  /* Runtime-accessible attributes exposed via C++ getters but
   * without an explicit typemap doc= annotation.  */
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
}

} // namespace casadi
#endif // SWIGPYTHON

/* There is no reason to expose the Slice class to e.g. Python or MATLAB. Only if an interfaced language
   lacks a slice type, the type should be exposed here */
// #if !(defined(SWIGPYTHON) || defined(SWIGMATLAB))
%include <casadi/core/slice.hpp>
 //#endif


%include <casadi/core/generic_matrix.hpp>

%template(GenDM)        casadi::GenericMatrix<casadi::Matrix<double> >;
%template(GenSX)             casadi::GenericMatrix<casadi::Matrix<casadi::SXElem> >;
%template(GenMX)             casadi::GenericMatrix<casadi::MX>;

%include <casadi/core/generic_expression.hpp>

// Flags to allow differentiating the wrapping by type
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

 DECL M casadi_horzcat(const std::vector< M > &v) {
  return horzcat(v);
 }
 DECL M casadi_vertcat(const std::vector< M > &v) {
 return vertcat(v);
 }
 DECL std::vector< M >
 casadi_horzsplit(const M& v, const std::vector<casadi_int>& offset) {
 return horzsplit(v, offset);
 }
 DECL std::vector< M > casadi_horzsplit(const M& v, casadi_int incr=1) {
 return horzsplit(v, incr);
 }
 DECL std::vector< M > casadi_horzsplit_n(const M& v, casadi_int n) {
 return horzsplit_n(v, n);
 }
 DECL std::vector< M >
 casadi_vertsplit(const M& v, const std::vector<casadi_int>& offset) {
 return vertsplit(v, offset);
 }
 DECL std::vector<casadi_int >
 casadi_offset(const std::vector< M > &v, bool vert=true) {
 return offset(v, vert);
 }
 DECL std::vector< M >
 casadi_vertsplit(const M& v, casadi_int incr=1) {
 return vertsplit(v, incr);
 }
 DECL std::vector< M >
 casadi_vertsplit_n(const M& v, casadi_int n) {
 return vertsplit_n(v, n);
 }
 DECL M casadi_blockcat(const M& A, const M& B, const M& C, const M& D) {
 return vertcat(horzcat(A, B), horzcat(C, D));
 }
 DECL std::vector< std::vector< M > >
 casadi_blocksplit(const M& x, const std::vector<casadi_int>& vert_offset,
 const std::vector<casadi_int>& horz_offset) {
 return blocksplit(x, vert_offset, horz_offset);
 }
 DECL std::vector< std::vector< M > >
 casadi_blocksplit(const M& x, casadi_int vert_incr=1, casadi_int horz_incr=1) {
 return blocksplit(x, vert_incr, horz_incr);
 }
 DECL M casadi_diagcat(const std::vector< M > &A) {
 return diagcat(A);
 }
 DECL std::vector< M >
 casadi_diagsplit(const M& x, const std::vector<casadi_int>& output_offset1,
 const std::vector<casadi_int>& output_offset2) {
 return diagsplit(x, output_offset1, output_offset2);
 }
 DECL std::vector< M >
 casadi_diagsplit(const M& x, const std::vector<casadi_int>& output_offset) {
 return diagsplit(x, output_offset);
 }
 DECL std::vector< M > casadi_diagsplit(const M& x, casadi_int incr=1) {
 return diagsplit(x, incr);
 }
 DECL std::vector< M >
 casadi_diagsplit(const M& x, casadi_int incr1, casadi_int incr2) {
 return diagsplit(x, incr1, incr2);
 }
 DECL M casadi_veccat(const std::vector< M >& x) {
 return veccat(x);
 }
 DECL M casadi_mtimes(const M& x, const M& y,
                      const std::string& blas = "reference") {
 return mtimes(x, y, blas);
 }
 DECL M casadi_mtimes(const std::vector< M > &args) {
 return mtimes(args);
 }
 DECL M casadi_mtimes(const std::vector< M > &args, const std::string& blas) {
 return mtimes(args, blas);
 }
 DECL M casadi_mac(const M& X, const M& Y, const M& Z,
                   const std::string& blas = "reference") {
 return mac(X, Y, Z, blas);
 }
 DECL M casadi_transpose(const M& X) {
 return X.T();
 }
 DECL M casadi_vec(const M& a) {
 return vec(a);
 }
 DECL M casadi_reshape(const M& a, casadi_int nrow, casadi_int ncol) {
 return reshape(a, nrow, ncol);
 }
 DECL M casadi_reshape(const M& a, std::pair<casadi_int, casadi_int> rc) {
 return reshape(a, rc.first, rc.second);
 }
 DECL M casadi_reshape(const M& a, const Sparsity& sp) {
 return reshape(a, sp);
 }
 DECL M casadi_sparsity_cast(const M& a, const Sparsity& sp) {
 return sparsity_cast(a, sp);
 }
 DECL casadi_int casadi_sprank(const M& A) {
 return sprank(A);
 }
 DECL casadi_int casadi_norm_0_mul(const M& x, const M& y) {
 return norm_0_mul(x, y);
 }
 DECL M casadi_triu(const M& a, bool includeDiagonal=true) {
 return triu(a, includeDiagonal);
 }
 DECL M casadi_tril(const M& a, bool includeDiagonal=true) {
 return tril(a, includeDiagonal);
 }
 DECL M casadi_kron(const M& a, const M& b) {
 return kron(a, b);
 }
 DECL M casadi_repmat(const M& A, casadi_int n, casadi_int m=1) {
 return repmat(A, n, m);
 }
 DECL M casadi_repmat(const M& A, const std::pair<casadi_int, casadi_int>& rc) {
 return repmat(A, rc.first, rc.second);
 }
 DECL M casadi_sum2(const M& x) {
 return sum2(x);
 }
 DECL M casadi_sum1(const M& x) {
 return sum1(x);
 }
#endif
%enddef

%define SPARSITY_INTERFACE_ALL(DECL, FLAG)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_SPARSITY), Sparsity)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_MX), MX)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_DMATRIX), Matrix<double>)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_SX), Matrix<SXElem>)
%enddef

#ifdef SWIGMATLAB
  %define SPARSITY_INTERFACE_FUN(DECL, FLAG, M)
    SPARSITY_INTERFACE_FUN_BASE(DECL, FLAG, M)
    #if FLAG & IS_MEMBER
      DECL casadi_int casadi_length(const M &v) {
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
    #endif
  %enddef
#else
  %define SPARSITY_INTERFACE_FUN(DECL, FLAG, M)
    SPARSITY_INTERFACE_FUN_BASE(DECL, FLAG, M)
    #if FLAG & IS_MEMBER
      DECL M casadi_sum(const M& x, casadi_int dim) {
        if (dim==0) return sum1(x);
        if (dim==1) return sum2(x);
        casadi_error(
          "Expected sum(A,0), sum(A,1), sum(A,\"all\") got " + casadi::str(dim) + " instead.");
      }
      DECL M casadi_sum(const M& x) {
        return sum(x);
      }
    #endif
  %enddef
#endif

%define GENERIC_MATRIX_FUN(DECL, FLAG, M)
#if FLAG & IS_MEMBER
DECL M casadi_mpower(const M& x, const M& n) {
  return mpower(x, n);
}

DECL M casadi_mrdivide(const M& x, const M& y) {
  return mrdivide(x, y);
}

DECL M casadi_mldivide(const M& x, const M& y) {
  return mldivide(x, y);
}

DECL std::vector< M > casadi_symvar(const M& x) {
  return symvar(x);
}

DECL M casadi_transform(const M& x, const casadi::Dict& opts = casadi::Dict()) {
  return transform(x, opts);
}

DECL M casadi_transform(const M& x, const std::vector< std::vector< casadi::GenericType > >& passes, const casadi::Dict& opts = casadi::Dict()) {
  return transform(x, passes, opts);
}

DECL std::vector< M > casadi_transform(const std::vector< M >& x, const casadi::Dict& opts = casadi::Dict()) {
  return transform(x, opts);
}

DECL std::vector< M > casadi_transform(const std::vector< M >& x, const std::vector< std::vector< casadi::GenericType > >& passes, const casadi::Dict& opts = casadi::Dict()) {
  return transform(x, passes, opts);
}

DECL M casadi_bilin(const M& A, const M& x, const M& y) {
  return bilin(A, x, y);
}

DECL M casadi_bilin(const M& A, const M& x) {
  return bilin(A, x);
}

DECL M casadi_rank1(const M& A, const M& alpha, const M& x, const M& y) {
  return rank1(A, alpha, x, y);
}

DECL M casadi_sumsqr(const M& X) {
  return sumsqr(X);
}

DECL M casadi_linspace(const M& a, const M& b, casadi_int nsteps) {
  return linspace(a, b, nsteps);
}

DECL M casadi_logsumexp(const M& a) {
  return logsumexp(a);
}

DECL M casadi_logsumexp(const M& a, const M& margin) {
  return logsumexp(a, margin);
}

DECL M casadi_interp1d(const std::vector<double>& x, const M&v,
        const std::vector<double>& xq, const std::string& mode="linear", bool equidistant=false) {
  return interp1d(x, v, xq, mode, equidistant);
}

DECL M casadi_soc(const M& x, const M& y) {
  return soc(x, y);
}

DECL M casadi_cross(const M& a, const M& b, casadi_int dim = -1) {
  return cross(a, b, dim);
}

DECL M casadi_skew(const M& a) {
  return skew(a);
}

DECL M casadi_inv_skew(const M& a) {
  return inv_skew(a);
}

DECL M casadi_det(const M& A) {
  return det(A);
}

DECL M casadi_det(const M& A, const std::string& lsolver,
                      const casadi::Dict& opts = casadi::Dict()) {
  return det(A, lsolver, opts);
}

DECL M casadi_inv_minor(const M& A) {
  return inv_minor(A);
}

DECL M casadi_inv(const M& A) {
  return inv(A);
}

DECL M casadi_inv(const M& A, const std::string& lsolver,
                      const casadi::Dict& opts = casadi::Dict()) {
  return inv(A, lsolver, opts);
}

DECL M casadi_trace(const M& a) {
  return trace(a);
}

DECL M casadi_tril2symm(const M& a) {
  return tril2symm(a);
}

DECL M casadi_triu2symm(const M& a) {
  return triu2symm(a);
}

DECL M casadi_norm_fro(const M& x) {
  return norm_fro(x);
}

DECL M casadi_norm_2(const M& x) {
  return norm_2(x);
}

DECL M casadi_norm_1(const M& x) {
  return norm_1(x);
}

DECL M casadi_norm_inf(const M& x) {
  return norm_inf(x);
}

DECL M casadi_dot(const M& x, const M& y) {
  return dot(x, y);
}

DECL M casadi_nullspace(const M& A) {
  return nullspace(A);
}

DECL M casadi_polyval(const M& p, const M& x) {
  return polyval(p, x);
}

DECL M casadi_diag(const M& A) {
  return diag(A);
}

DECL M casadi_unite(const M& A, const M& B) {
  return unite(A, B);
}

DECL M casadi_densify(const M& x) {
  return densify(x);
}

DECL M casadi_project(const M& A, const Sparsity& sp, bool intersect=false) {
  return project(A, sp, intersect);
}

DECL M casadi_if_else(const M& cond, const M& if_true,
                    const M& if_false, bool short_circuit=false) {
  return if_else(cond, if_true, if_false, short_circuit);
}

DECL M casadi_conditional(const M& ind, const std::vector< M > &x,
                        const M& x_default, bool short_circuit=false) {
  return conditional(ind, x, x_default, short_circuit);
}

DECL bool casadi_depends_on(const M& f, const M& arg) {
  return depends_on(f, arg);
}

DECL bool casadi_contains(const std::vector<M>& v, const M& n) {
  return contains(v, n);
}

DECL bool casadi_contains_all(const std::vector<M>& v, const std::vector<M>& n) {
  return contains_all(v, n);
}

DECL bool casadi_contains_any(const std::vector<M>& v, const std::vector<M>& n) {
  return contains_any(v, n);
}

DECL M casadi_solve(const M& A, const M& b) {
  return solve(A, b);
}

DECL M casadi_solve(const M& A, const M& b,
                       const std::string& lsolver,
                       const casadi::Dict& opts = casadi::Dict()) {
  return solve(A, b, lsolver, opts);
}

DECL M casadi_pinv(const M& A) {
  return pinv(A);
}

DECL M casadi_pinv(const M& A, const std::string& lsolver,
                      const casadi::Dict& opts = casadi::Dict()) {
  return pinv(A, lsolver, opts);
}

DECL M casadi_expm_const(const M& A, const M& t) {
  return expm_const(A, t);
}

DECL M casadi_expm(const M& A) {
  return expm(A);
}

DECL M casadi_jacobian(const M &ex, const M &arg, const Dict& opts=Dict()) {
  return jacobian(ex, arg, opts);
}

DECL M casadi_jtimes(const M& ex, const M& arg, const M& v, bool tr=false, const Dict& opts=Dict()) {
  return jtimes(ex, arg, v, tr);
}

DECL M casadi_linearize(const M& f, const M& x, const M& x0) {
  return linearize(f, x, x0);
}

DECL std::vector<bool> casadi_which_depends(const M& expr, const M& var,
                                            casadi_int order=1, bool tr=false) {
  return which_depends(expr, var, order, tr);
}

DECL Sparsity casadi_jacobian_sparsity(const M& f, const M& x) {
  return jacobian_sparsity(f, x);
}

DECL bool casadi_is_linear(const M& expr, const M& var) {
  return is_linear(expr, var);
}

DECL bool casadi_is_quadratic(const M& expr, const M& var) {
  return is_quadratic(expr, var);
}

DECL M casadi_gradient(const M &ex, const M &arg, const Dict& opts=Dict()) {
  return gradient(ex, arg, opts);
}

DECL M casadi_tangent(const M &ex, const M &arg, const Dict& opts=Dict()) {
  return tangent(ex, arg, opts);
}

DECL M casadi_hessian(const M& ex, const M& arg, M& OUTPUT1, const casadi::Dict& opts = casadi::Dict()) {
  return hessian(ex, arg, OUTPUT1, opts);
}

DECL void casadi_quadratic_coeff(const M& ex, const M& arg, M& OUTPUT1, M& OUTPUT2, M& OUTPUT3, bool check=true) {
  quadratic_coeff(ex, arg, OUTPUT1, OUTPUT2, OUTPUT3, check);
}

DECL void casadi_linear_coeff(const M& ex, const M& arg, M& OUTPUT1, M& OUTPUT2, bool check=true) {
  linear_coeff(ex, arg, OUTPUT1, OUTPUT2, check);
}

DECL casadi_int casadi_n_nodes(const M& A) {
  return n_nodes(A);
}

DECL std::string casadi_print_operator(const M& xb,
                                                  const std::vector<std::string>& args) {
  return print_operator(xb, args);
}
DECL M casadi_repsum(const M& A, casadi_int n, casadi_int m=1) {
  return repsum(A, n, m);
}
DECL M casadi_diff(const M& A, casadi_int n=1, casadi_index axis=-1) {
  return diff(A, n, axis);
}
DECL M casadi_cumsum(const M& A, casadi_index axis=-1) {
  return cumsum(A, axis);
}
DECL M casadi_einstein(const M& A, const M& B, const M& C,
  const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b, const std::vector<casadi_int>& dim_c,
  const std::vector<casadi_int>& a, const std::vector<casadi_int>& b, const std::vector<casadi_int>& c) {
  return einstein(A, B, C, dim_a, dim_b, dim_c, a, b, c);
}
DECL M casadi_einstein(const M& A, const M& B,
  const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b, const std::vector<casadi_int>& dim_c,
  const std::vector<casadi_int>& a, const std::vector<casadi_int>& b, const std::vector<casadi_int>& c) {
  return einstein(A, B, dim_a, dim_b, dim_c, a, b, c);
}
DECL M casadi_mmin(const M& x) { return mmin(x); }
DECL M casadi_mmax(const M& x) { return mmax(x); }
DECL casadi::DM casadi_evalf(const M& x) {
  return evalf(x);
}
DECL void casadi_separate_linear(const M &expr,
      const M &sym_lin, const M &sym_const,
      M& OUTPUT1, M& OUTPUT2, M& OUTPUT3) {
  separate_linear(expr, sym_lin, sym_const, OUTPUT1, OUTPUT2, OUTPUT3);
}
DECL void casadi_separate_linear(const M &expr,
  const std::vector<M> &sym_lin, const std::vector<M> &sym_const,
  M& OUTPUT1, M& OUTPUT2, M& OUTPUT3) {
separate_linear(expr, sym_lin, sym_const, OUTPUT1, OUTPUT2, OUTPUT3);
}
#endif // FLAG & IS_MEMBER

#if FLAG & IS_GLOBAL
DECL std::vector<M> casadi_cse(const std::vector<M>& e) {
  return cse(e);
}
DECL M casadi_cse(const M& e) {
  return cse(e);
}

DECL void casadi_extract_parametric(const M &expr, const M& par,
        M& OUTPUT1, std::vector<M>& OUTPUT2, std::vector<M>& OUTPUT3, const Dict& opts=Dict()) {
  extract_parametric(expr, par, OUTPUT1, OUTPUT2, OUTPUT3, opts);
}
DECL void casadi_extract_parametric(const M &expr, const std::vector<M>& par,
        M& OUTPUT1, std::vector<M>& OUTPUT2, std::vector<M>& OUTPUT3, const Dict& opts=Dict()) {
  extract_parametric(expr, par, OUTPUT1, OUTPUT2, OUTPUT3, opts);
}
DECL void casadi_extract_parametric(const std::vector<M> &expr, const M& par,
        std::vector<M>& OUTPUT1, std::vector<M>& OUTPUT2, std::vector<M>& OUTPUT3, const Dict& opts=Dict()) {
  extract_parametric(expr, par, OUTPUT1, OUTPUT2, OUTPUT3, opts);
}
DECL void casadi_extract_parametric(const std::vector<M> &expr, const std::vector<M>& par,
        std::vector<M>& OUTPUT1, std::vector<M>& OUTPUT2, std::vector<M>& OUTPUT3, const Dict& opts=Dict()) {
  extract_parametric(expr, par, OUTPUT1, OUTPUT2, OUTPUT3, opts);
}

DECL std::vector<std::vector< M > >
casadi_forward(const std::vector< M > &ex, const std::vector< M > &arg,
               const std::vector<std::vector< M > > &v,
               const Dict& opts = Dict()) {
  return forward(ex, arg, v, opts);
}

DECL std::vector<std::vector< M > >
casadi_reverse(const std::vector< M > &ex, const std::vector< M > &arg,
               const std::vector<std::vector< M > > &v,
               const Dict& opts = Dict()) {
  return reverse(ex, arg, v, opts);
}

DECL M casadi_substitute(const M& ex, const M& v, const M& vdef) {
  return substitute(ex, v, vdef);
}

DECL std::vector< M > casadi_substitute(const std::vector< M >& ex,
                                         const std::vector< M >& v,
                                         const std::vector< M >& vdef) {
  return substitute(ex, v, vdef);
}

DECL void casadi_substitute_inplace(const std::vector< M >& v,
                                      std::vector< M >& INOUT1,
                                      std::vector< M >& INOUT2,
                                      bool reverse=false) {
  return substitute_inplace(v, INOUT1, INOUT2, reverse);
}

DECL void casadi_extract(const std::vector< M >& ex,
    std::vector< M >& OUTPUT1,
    std::vector< M >& OUTPUT2,
    std::vector< M >& OUTPUT3,
    const Dict& opts = Dict()) {
  OUTPUT1 = ex;
  extract(OUTPUT1, OUTPUT2, OUTPUT3, opts);
}

DECL void casadi_shared(const std::vector< M >& ex,
                               std::vector< M >& OUTPUT1,
                               std::vector< M >& OUTPUT2,
                               std::vector< M >& OUTPUT3,
                               const std::string& v_prefix="v_",
                               const std::string& v_suffix="") {
  OUTPUT1 = ex;
  shared(OUTPUT1, OUTPUT2, OUTPUT3, v_prefix, v_suffix);
}

DECL M casadi_blockcat(const std::vector< std::vector< M > > &v) {
 return blockcat(v);
}
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
DECL M casadi_all(const M& x) {
  return all(x);
}

DECL M casadi_any(const M& x) {
  return any(x);
}

DECL M casadi_adj(const M& A) {
  return adj(A);
}

DECL M casadi_minor(const M& x, casadi_int i, casadi_int j) {
  return minor(x, i, j);
}

DECL M casadi_cofactor(const M& x, casadi_int i, casadi_int j) {
  return cofactor(x, i, j);
}

DECL void casadi_qr(const M& A, M& OUTPUT1, M& OUTPUT2) {
  return qr(A, OUTPUT1, OUTPUT2);
}

DECL void casadi_qr_sparse(const M& A, M& OUTPUT1, M& OUTPUT2, M& OUTPUT3,
          std::vector<casadi_int>& OUTPUT4, std::vector<casadi_int>& OUTPUT5, bool amd=true) {
  return qr_sparse(A, OUTPUT1, OUTPUT2, OUTPUT3, OUTPUT4, OUTPUT5, amd);
}

DECL M casadi_qr_solve(const M& b, const M& v, const M& r, const M& beta,
                       const std::vector<casadi_int>& prinv,
                       const std::vector<casadi_int>& pc, bool tr=false) {
  return qr_solve(b, v, r, beta, prinv, pc, tr);
}

DECL void casadi_ldl(const M& A, M& OUTPUT1, M& OUTPUT2, std::vector<casadi_int>& OUTPUT3, bool amd=true) {
  return ldl(A, OUTPUT1, OUTPUT2, OUTPUT3, amd);
}

DECL M casadi_ldl_solve(const M& b, const M& D, const M& LT, const std::vector<casadi_int>& p) {
  return ldl_solve(b, D, LT, p);
}

DECL M casadi_chol(const M& A) {
  return chol(A);
}

DECL M casadi_norm_inf_mul(const M& x, const M& y) {
  return norm_inf_mul(x, y);
}

DECL M casadi_sparsify(const M& A, double tol=0) {
  return sparsify(A, tol);
}

DECL void casadi_expand(const M& ex, M& OUTPUT1, M& OUTPUT2) {
  expand(ex, OUTPUT1, OUTPUT2);
}

DECL M casadi_pw_const(const M &t, const M& tval, const M& val) {
  return pw_const(t, tval, val);
}

DECL M casadi_pw_lin(const M& t, const M& tval, const M& val) {
  return pw_lin(t, tval, val);
}

DECL M casadi_heaviside(const M& x) {
  return heaviside(x);
}

DECL M casadi_rectangle(const M& x) {
  return rectangle(x);
}

DECL M casadi_triangle(const M& x) {
  return triangle(x);
}

DECL M casadi_ramp(const M& x) {
  return ramp(x);
}

DECL M casadi_gauss_quadrature(const M& f, const M& x,
                               const M& a, const M& b,
                               casadi_int order=5) {
  return gauss_quadrature(f, x, a, b, order);
}

DECL M casadi_gauss_quadrature(const M& f, const M& x,
                               const M& a, const M& b,
                               casadi_int order, const M& w) {
  return gauss_quadrature(f, x, a, b, order, w);
}

DECL M casadi_taylor(const M& ex, const M& x, const M& a=0, casadi_int order=1) {
  return taylor(ex, x, a, order);
}

DECL M casadi_mtaylor(const M& ex, const M& x, const M& a, casadi_int order=1) {
  return mtaylor(ex, x, a, order);
}

DECL M casadi_mtaylor(const M& ex, const M& x, const M& a, casadi_int order,
                      const std::vector<casadi_int>& order_contributions) {
  return mtaylor(ex, x, a, order, order_contributions);
}

DECL M casadi_poly_coeff(const M& ex,
                         const M&x) {
  return poly_coeff(ex, x);
}

DECL M casadi_poly_roots(const M& p) {
  return poly_roots(p);
}

DECL M casadi_eig_symbolic(const M& m) {
  return eig_symbolic(m);
}

#endif
%enddef

%define MATRIX_ALL(DECL, FLAG)
MATRIX_FUN(DECL, (FLAG | IS_DMATRIX), Matrix<double>)
MATRIX_FUN(DECL, (FLAG | IS_SX), Matrix<SXElem>)
%enddef

%define MX_FUN(DECL, FLAG, M)
#if FLAG & IS_MEMBER
DECL M casadi_find(const M& x) {
  return find(x);
}
DECL M casadi_low(const M& v, const M& p, const Dict& options = Dict()) {
  return low(v, p, options);
}
DECL M casadi_inv_node(const M& x) {
  return inv_node(x);
}
#endif // FLAG & IS_MEMBER

#if FLAG & IS_GLOBAL
DECL std::vector< M >
casadi_matrix_expand(const std::vector< M >& e,
                     const std::vector< M > &boundary = std::vector< M >(),
                     const Dict& options = Dict()) {
  return matrix_expand(e, boundary, options);
}

DECL M casadi_matrix_expand(const M& e,
                            const std::vector< M > &boundary = std::vector< M >(),
                            const Dict& options = Dict()) {
  return matrix_expand(e, boundary, options);
}

DECL M casadi_graph_substitute(const M& ex, const std::vector< M >& v,
                         const std::vector< M > &vdef) {
  return graph_substitute(ex, v, vdef);
}

DECL std::vector< M >
casadi_graph_substitute(const std::vector< M > &ex,
                 const std::vector< M > &v,
                 const std::vector< M > &vdef) {
  return graph_substitute(ex, v, vdef);
}
DECL M casadi_bspline(const M& x,
        const DM& coeffs,
        const std::vector< std::vector<double> >& knots,
        const std::vector<casadi_int>& degree,
        casadi_int m,
        const Dict& opts = Dict()) {
  return bspline(x, coeffs, knots, degree, m, opts);
}
DECL M casadi_bspline(const M& x,
        const M& coeffs,
        const std::vector< std::vector<double> >& knots,
        const std::vector<casadi_int>& degree,
        casadi_int m,
        const Dict& opts = Dict()) {
  return bspline(x, coeffs, knots, degree, m, opts);
}
DECL MX casadi_bspline(const MX& x,
        const MX& coeffs,
        const std::vector<MX>& knots,
        const std::vector<casadi_int>& degree,
        casadi_int m,
        const Dict& opts = Dict()) {
  return bspline(x, coeffs, knots, degree, m, opts);
}
DECL M casadi_convexify(const M& H,
        const Dict& opts = Dict()) {
  return convexify(H, opts);
}
DECL M casadi_stop_diff(const M& expr, casadi_int order) {
  return stop_diff(expr, order);
}
DECL M casadi_stop_diff(const M& expr, const M& var, casadi_int order) {
  return stop_diff(expr, var, order);
}
DECL std::vector< M > casadi_difference(const std::vector< M >& a, const std::vector< M >& b) {
  return difference(a, b);
}
DECL M casadi_no_hess(const M& expr) {
  return no_hess(expr);
}
DECL M casadi_no_grad(const M& expr) {
  return no_grad(expr);
}
DECL M casadi_kron_contract(const M& m, const M& x, bool inner) {
  return kron_contract(m, x, inner);
}

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

%template(DM) casadi::Matrix<double>;
%extend casadi::Matrix<double> {
   %template(DM) Matrix<SXElem>;
};


namespace casadi{
  %extend Matrix<double> {
#ifndef SWIGWASMJS
    void assign(const casadi::Matrix<double>&rhs) { (*$self)=rhs; }
#endif
    %matrix_helpers(casadi::Matrix<double>)

  }

}

#ifdef SWIGPYTHON
  %feature("nothread") casadi::Matrix<double>::full;
  %feature("nothread") casadi::Matrix<double>::sparse;
#endif

// Extend DM with SWIG unique features
#ifndef SWIGWASMJS
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

%python_array_wrappers(999.0)

%pythoncode %{
  def tocsc(self):
    import numpy as np
    import warnings
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      from scipy.sparse import csc_matrix
    return csc_matrix( (self.nonzeros(),self.row(),self.colind()), shape = self.shape, dtype=np.double )
  def toarray(self,simplify=False):
    import numpy as np
    if simplify:
      if self.is_scalar():
        return float(self)
      elif self.is_vector():
        return np.array(self.T.elements())
    return np.array(self.T.elements()).reshape(self.shape)
  def to_numpy(self):
    # Pandas-style densification hook.  matplotlib's
    # cbook._unpack_to_numpy looks for `to_numpy()` before falling back
    # to iteration, so defining it here makes plt.plot(DM) work without
    # forcing densification through NEP-18 shape ops like np.atleast_1d.
    return self.full()
%}
/* `tocsc` returns a scipy.sparse.csc_matrix; not type-importing scipy.sparse
 * here keeps the stub independent of an optional dependency.  `Any` is the
 * least-bad option until we also bridge scipy.sparse types. */
%stub_method0(tocsc,    Any)
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

}


%extend Function {

  %pythoncode %{
    def __setstate__(self, state):
        self.__init__(Function.deserialize(state["serialization"]))

    def __getstate__(self):
        return {"serialization": self.serialize()}
  %}

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

#endif // SWIGPYTHON

#ifdef SWIGPYTHON
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
    %matrix_helpers(casadi::Matrix<casadi::SXElem>)

  #ifdef SWIGPYTHON
  %python_array_wrappers(1001.0)
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
  %matrix_helpers(casadi::MX)
  #ifdef SWIGPYTHON
  %python_array_wrappers(1002.0)
  /* MX is unique among CasadiMatrix types in accepting MX-valued
   * indices (runtime: MX_get / MX_set).  Additional overloads on top
   * of the DM-compatible ones from %matrix_helpers. */
  %stub_overload_method(__getitem__, MX, s: MX | tuple[_MIndex | MX, _MIndex | MX])
  %stub_overload_method(__setitem__, None, s: MX | tuple[_MIndex | MX, _MIndex | MX], val: bool | int | float | DM | SX | MX | Sequence[bool | int | float] | Sequence[Sequence[bool | int | float]] | NDArray[Any])
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
  /* buffer() returns (FunctionBuffer, callable) -- the callable is a
   * functools.partial around the raw eval function, typed loosely. */
  %stub_method0(buffer, %arg(tuple[FunctionBuffer, Any]))
  /* Function.__call__ -- narrowing overloads per homogeneous input
   * matrix type so pyright resolves f(mx) -> MX.  Runtime semantics
   * (see the %pythoncode above):
   *   f(pos_args)   -> single matrix if 1 output,  tuple[matrix,...] if >1
   *   f(**kwargs)   -> dict[str, matrix]
   *   f()           -> dict (no-arg path routes to `self.call({})`)
   * Output arity is not statically known.  Returning the single-matrix
   * case covers the common pattern (one-output Function); multi-output
   * sites use `f.call([...])` which returns the tuple explicitly.
   * The positional overloads require at least one positional arg
   * (arg0 pos-only) so f() / f(**kw) routes to the dict overloads. */
  %stub_overload_method(__call__, DM, arg0: _DM, /, *args: _DM)
  %stub_overload_method(__call__, SX, arg0: _SX, /, *args: _SX)
  %stub_overload_method(__call__, MX, arg0: _MX, /, *args: _MX)
  %stub_overload_method(__call__, %arg(dict[str, DM]), **kwargs: _DM)
  %stub_overload_method(__call__, %arg(dict[str, SX]), **kwargs: _SX)
  %stub_overload_method(__call__, %arg(dict[str, MX]), **kwargs: _MX)
  %stub_overload_method(__call__, %arg(dict[str, DM | SX | MX]), **kwargs: bool | int | float | DM | SX | MX)
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
      if numel(s)==1 && strcmp(s.type,'()')
        [varargout{1}] = paren(self, s.subs{:});
      elseif numel(s)==1 && strcmp(s.type,'{}')
        [varargout{1}] = brace(self, s.subs{:});
      else
        [varargout{1:nargout}] = builtin('subsref',self,s);
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
%include <casadi/core/rootfinder.hpp>
%include <casadi/core/linsol.hpp>
%include <casadi/core/dple.hpp>
%include <casadi/core/expm.hpp>
%include <casadi/core/interpolant.hpp>
%include <casadi/core/blazing_spline.hpp>

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
    const MX paren(char rr, const casadi::MX& cc) const {
      casadi::MX m;
      $self->get(m, true, casadi::char2Slice(rr), cc);
      return m;
    }
    const MX paren(const casadi::MX& rr, char cc) const {
      casadi::MX m;
      $self->get(m, true, rr, casadi::char2Slice(cc));
      return m;
    }
    const MX paren(const casadi::MX& rr, const casadi::MX& cc) const {
      casadi::MX m;
      $self->get(m, true, rr, cc);
      return m;
    }
    /*
    Not yet implemeted in core
    set(const MX& m, bool ind1, const MX&, const MX&); does not seem to exist
    void paren_asgn(const MX& m, char rr, const casadi::MX& cc) {
      $self->set(m, true, casadi::char2Slice(rr), cc);
    }
    void paren_asgn(const MX& m, const casadi::MX& rr, char cc) {
      $self->set(m, true, rr, casadi::char2Slice(cc));
    }
    void paren_asgn(const MX& m, const casadi::MX& rr, const casadi::MX& cc) {
      $self->set(m, true, rr, cc);
    }*/
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
    %stub_method0(__hash__, int)
    %stub_CasadiMatrix_binop(__matmul__)
    %stub_CasadiMatrix_binop(__rmatmul__)
    %stub_CasadiMatrix_binop(__imatmul__)
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
      def fmin(x, y): return _casadi.fmin(x, y)
      def fmax(x, y): return _casadi.fmax(x, y)
      def __fmin__(x, y): return _casadi.fmin(x, y)
      def __rfmin__(y, x): return _casadi.fmin(x, y)
      def __fmax__(x, y): return _casadi.fmax(x, y)
      def __rfmax__(y, x): return _casadi.fmax(x, y)
      def logic_and(x, y): return _casadi.logic_and(x, y)
      def logic_or(x, y): return _casadi.logic_or(x, y)
      def fabs(x): return _casadi.fabs(x)
      def sqrt(x): return _casadi.sqrt(x)
      def sin(x): return _casadi.sin(x)
      def cos(x): return _casadi.cos(x)
      def tan(x): return _casadi.tan(x)
      def arcsin(x): return _casadi.asin(x)
      def arccos(x): return _casadi.acos(x)
      def arctan(x): return _casadi.atan(x)
      def sinh(x): return _casadi.sinh(x)
      def cosh(x): return _casadi.cosh(x)
      def tanh(x): return _casadi.tanh(x)
      def arcsinh(x): return _casadi.asinh(x)
      def arccosh(x): return _casadi.acosh(x)
      def arctanh(x): return _casadi.atanh(x)
      def exp(x): return _casadi.exp(x)
      def log(x): return _casadi.log(x)
      def log10(x): return _casadi.log10(x)
      def log1p(x): return _casadi.log1p(x)
      def expm1(x): return _casadi.expm1(x)
      def floor(x): return _casadi.floor(x)
      def ceil(x): return _casadi.ceil(x)
      def erf(x): return _casadi.erf(x)
      def sign(x): return _casadi.sign(x)
      def fmod(x, y): return _casadi.mod(x, y)
      def hypot(x, y): return _casadi.hypot(x, y)
      def remainder(x, y): return _casadi.remainder(x, y)
      def __copysign__(x, y): return _casadi.copysign(x, y)
      def __rcopysign__(y, x): return _casadi.copysign(x, y)
      def copysign(x, y): return _casadi.copysign(x, y)
      def rcopysign(y, x): return _casadi.copysign(x, y)
      def __constpow__(x, y): return _casadi.constpow(x, y)
      def __rconstpow__(y, x): return _casadi.constpow(x, y)
      def constpow(x, y): return _casadi.constpow(x, y)
      def rconstpow(y, x): return _casadi.constpow(x, y)
    %}
    /* PEP-484 stubs for the %pythoncode operators above. */
    %stub_CasadiMatrix_binop(__add__)
    %stub_CasadiMatrix_binop(__radd__)
    %stub_CasadiMatrix_binop(__sub__)
    %stub_CasadiMatrix_binop(__rsub__)
    %stub_CasadiMatrix_binop(__mul__)
    %stub_CasadiMatrix_binop(__rmul__)
    %stub_CasadiMatrix_binop(__truediv__)
    %stub_CasadiMatrix_binop(__rtruediv__)
    %stub_CasadiMatrix_binop(__floordiv__)
    %stub_CasadiMatrix_binop(__rfloordiv__)
    %stub_CasadiMatrix_binop(__mod__)
    %stub_CasadiMatrix_binop(__rmod__)
    %stub_CasadiMatrix_binop(__pow__)
    %stub_CasadiMatrix_binop(__rpow__)
    %stub_CasadiMatrix_binop(__iadd__)
    %stub_CasadiMatrix_binop(__isub__)
    %stub_CasadiMatrix_binop(__imul__)
    %stub_CasadiMatrix_binop(__itruediv__)
    %stub_CasadiMatrix_binop(__ifloordiv__)
    %stub_CasadiMatrix_binop(__imod__)
    %stub_CasadiMatrix_binop(__ipow__)
    %stub_CasadiMatrix_unop(__neg__)
    %stub_CasadiMatrix_unop(__pos__)
    %stub_CasadiMatrix_unop(__abs__)
    %stub_CasadiMatrix_unop(__trunc__)
    %stub_CasadiMatrix_unop(__floor__)
    %stub_CasadiMatrix_unop(__ceil__)
    %stub_CasadiMatrix_cmp(__lt__)
    %stub_CasadiMatrix_cmp(__le__)
    %stub_CasadiMatrix_cmp(__gt__)
    %stub_CasadiMatrix_cmp(__ge__)
    %stub_CasadiMatrix_cmp(__eq__)
    %stub_CasadiMatrix_cmp(__ne__)
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
    %stub_CasadiMatrix_binop(__mldivide__)
    %stub_CasadiMatrix_binop(__rmldivide__)
    %stub_CasadiMatrix_binop(__mrdivide__)
    %stub_CasadiMatrix_binop(__rmrdivide__)
    %stub_CasadiMatrix_binop(__mpower__)
    %stub_CasadiMatrix_binop(__rmpower__)
  }

} // namespace casadi
#endif // SWIGPYTHON

%feature("director") casadi::Callback;

/* The Python-visible Callback interface differs from the C++ signatures
 * because SWIG's director magic reshapes eval_buffer into a 2-tuple
 * form and hands the user-defined eval() a list of DM.  Suppress the
 * auto-generated .pyi entries for these methods and emit hand-crafted
 * ones below so pyright guides users toward correct overrides. */
#ifdef SWIGPYTHON
%feature("python:stub_skip") casadi::Callback::eval;
%feature("python:stub_skip") casadi::Callback::eval_buffer;
%feature("python:stub_skip") casadi::Callback::get_sparsity_in;
%feature("python:stub_skip") casadi::Callback::get_sparsity_out;
#endif

%include <casadi/core/importer.hpp>
%include <casadi/core/callback.hpp>

#ifdef SWIGPYTHON
%extend casadi::Callback {
  /* eval: Python-visible signature is `eval(self, arg, /)` where arg
   * is an iterable of DM.  Positional-only (arg0 can have any name
   * in subclasses).  Return type Sequence[_DM] so subclasses may
   * return anything that to_ptr<DM> promotes -- int / float / DM /
   * ndarray / nested sequences -- each gets converted at the C++
   * boundary. */
  %stub_method(eval, %arg(Sequence[_DM]), %arg(arg: Sequence[DM], /))
  /* eval_buffer: director reshapes the 4-arg C++ interface (arg,
   * sizes_arg, res, sizes_res) into a 2-tuple Python interface: two
   * tuples of memoryview objects exposing the structural nonzeros. */
  %stub_method(eval_buffer, int, %arg(arg: tuple[memoryview, ...], /, res: tuple[memoryview, ...]))
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

/* Pure-pythoncode classes -- SWIG doesn't see them, so the .pyi has to
 * be hand-rolled.  Both are referenced from error messages
 * (`Use ca.global_pickle_context(): ...`) and from test/python/serialize.py. */
%stub_class_begin(global_pickle_context)
%stub_method0(__enter__, StringSerializer)
%stub_method(__exit__, None, *args: Any)

%stub_class_begin(global_unpickle_context)
%stub_method0(__enter__, StringDeserializer)
%stub_method(__exit__, None, *args: Any)


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
/* JS-side analog of the Python / MATLAB %extend blocks above.  Same
   idea exactly: pop the type tag, name it via
   SerializerBase::type_to_string, dispatch to blind_unpack_<name> by
   dynamic property lookup.  Patching DeserializerBase.prototype covers
   both StringDeserializer and FileDeserializer (subclasses share the
   prototype chain).  Today %jscode (= %insert("js")) still lands in
   the global js slot, not the class body -- wasm_js.cxx does not yet
   wire %insert("js") inside %extend back into the per-class output --
   so the prototype assignment is explicit.  Once that wiring exists,
   this can drop the `DeserializerBase.prototype.` prefix and read
   `unpack = function() { ... }`. */
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

// Return-by-value
#ifdef SWIGWASMJS
/* wasm-js: full() / sparse() helpers aren't implemented for the wasm-js
   target (they return 0).  Treat native_DM as a regular DM-by-value
   so OptiSol::value(MX) etc. work: from_ref produces a {_ptr} EM_VAL
   carrier, the JS-side wrapping in Opti/OptiSol's `value` method
   then wraps it as `new DM(__PRIVATE_CTOR, ptr)`. */
%typemap(out, doc="DM", tsstub_out="DM", noblock=1, fragment="casadi_all") casadi::native_DM {
  if(!($result = casadi::from_ref($1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type 'DM'.");
}
#else
%typemap(out, doc="double", pystub_out="float", noblock=1, fragment="casadi_all") casadi::native_DM {
  if(!($result = full_or_sparse($1, true))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type 'double'.");
}
#endif


%apply casadi_int &OUTPUT { Opti::ConstraintType &OUTPUT };

%typemap(argout, noblock=1,fragment="casadi_all") casadi::Opti::ConstraintType &OUTPUT {
  %append_output(casadi::from_ptr((casadi_int *) $1));
}

%typemap(in, doc="Opti.ConstraintType", pystub_in="Any", noblock=1, numinputs=0) casadi::Opti::ConstraintType &OUTPUT (casadi::Opti::ConstraintType m) {
 $1 = &m;
}


#ifdef SWIGPYTHON

%define make_property(class, name, type)
  %rename(_ ## name) class ## :: ## name;
  %extend class {
    %pythoncode %{
      @property
      def name(self):
        return self._ ## name()
    %}
    %stub_attr(name, type)
  }
%enddef


make_property(casadi::Opti, debug, OptiAdvanced);
make_property(casadi::Opti, advanced, OptiAdvanced);
make_property(casadi::OptiSol, opti, Opti);

%define make_property_opti(name, type)
  make_property(casadi::Opti, name, type);
%enddef

make_property(casadi::OptiSol, debug, OptiAdvanced);
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

make_property(casadi::Opti, casadi_solver, Function);
%define opti_metadata_modifiers(class)
  %rename(_variable) class ## :: variable;
  %rename(_parameter) class ## :: parameter;
  %rename(_subject_to) class ## :: subject_to;
  %extend class {
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
    %stub_overload_method(variable, MX, n: int = ..., m: int = ..., attribute: str = ...)
    %stub_overload_method(variable, MX, sp: Sparsity | DM, attribute: str = ...)
    %stub_overload_method(variable, MX, symbol: MX | SX | DM, attribute: str = ...)
    %stub_overload_method(parameter, MX, n: int = ..., m: int = ..., attribute: str = ...)
    %stub_overload_method(parameter, MX, sp: Sparsity | DM, attribute: str = ...)
    %stub_overload_method(parameter, MX, symbol: MX | SX | DM, attribute: str = ...)
    /* subject_to(g[, linear_scale][, options]) -- the %pythoncode
     * above inspects args[1] to decide whether it's an options dict
     * or a linear_scale; the 3-arg form always has options last. */
    %stub_overload_method0(subject_to, None)
    %stub_overload_method(subject_to, None, g: MX | SX | DM | bool | int | float)
    %stub_overload_method(subject_to, None, g: Sequence[MX | SX | DM | bool | int | float])
    %stub_overload_method(subject_to, None, g: MX | SX | DM | bool | int | float, options: dict)
    %stub_overload_method(subject_to, None, g: Sequence[MX | SX | DM | bool | int | float], options: dict)
    %stub_overload_method(subject_to, None, g: MX | SX | DM | bool | int | float, linear_scale: _DM)
    %stub_overload_method(subject_to, None, g: Sequence[MX | SX | DM | bool | int | float], linear_scale: _DM)
    %stub_overload_method(subject_to, None, g: MX | SX | DM | bool | int | float, linear_scale: _DM, options: dict)
    %stub_overload_method(subject_to, None, g: Sequence[MX | SX | DM | bool | int | float], linear_scale: _DM, options: dict)
  }
%enddef

opti_metadata_modifiers(casadi::Opti);

#endif


#ifdef SWIGMATLAB
%define opti_metadata_modifiers(class)
  %rename(internal_variable) class ## ::variable;
  %rename(internal_parameter) class ## ::parameter;
  %rename(internal_subject_to) class ## ::subject_to;
  %extend class {
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
%enddef

opti_metadata_modifiers(casadi::Opti)
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
  /* Opti.callback registers a per-iteration Python callable; the %pythoncode
   * above defines it, so the SWIG C++ parser never sees it and pyright loses
   * visibility.  `fh` is invoked with the iteration index (int). */
  %stub_method(callback, None, %arg(fh: "Callable[[int], Any] | None" = ...))

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
#endif

#ifdef SWIGMATLAB
%{
#ifdef HAVE_OCTAVE
  // Mandatory as of Octave 10
  // Null effect for prior versions
  extern "C" const int __octave_mex_soversion__ = 1;
#endif
%}

#endif

%include <casadi/core/resource.hpp>

// Cleanup for dependent modules
%exception {
  $action
}
