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

  /// Data structure in the target language holding data
#ifdef SWIGPYTHON
#define GUESTOBJECT PyObject
#elif defined(SWIGMATLAB)
#define GUESTOBJECT mxArray
#else
#define GUESTOBJECT void
#endif

// Define printing routine

#ifdef SWIGPYTHON
%{
  namespace casadi {
    // Redirect printout
    static void pythonlogger(const char* s, std::streamsize num, bool error) {
      if (error) {
        PySys_WriteStderr("%.*s", static_cast<int>(num), s);
      } else {
        PySys_WriteStdout("%.*s", static_cast<int>(num), s);
      }
    }

    static bool pythoncheckinterrupted() {
      if (!casadi::InterruptHandler::is_main_thread()) return false;
      return PyErr_CheckSignals();
    }

    void handle_director_exception() {
	    std::string msg = "Exception in SWIG director ";
      SWIG_PYTHON_THREAD_BEGIN_BLOCK;
      if (PyErr_ExceptionMatches(PyExc_KeyboardInterrupt)) {
        PyErr_Clear();
        SWIG_PYTHON_THREAD_END_BLOCK;
        throw casadi::KeyboardInterruptException();
      }
      PyObject *ptype, *pvalue, *ptraceback;
      PyErr_Fetch(&ptype, &pvalue, &ptraceback);
      PyObject* msg_py = PyObject_Str(pvalue);
      char *msg_char = SWIG_Python_str_AsChar(msg_py);
      msg = msg_char;
      SWIG_Python_str_DelForPy3(msg_char);
      Py_DECREF(msg_py);
      PyErr_Restore(ptype, pvalue, ptraceback);
      PyErr_Print();
      SWIG_PYTHON_THREAD_END_BLOCK;
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

class _copyableObject(object):
  def __copy__(self):
    return self.__class__(self)

  def __deepcopy__(self,dummy=None):
    return self.__class__(self)

_object = object = _copyableObject

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
    if m.dtype!=np.object: return None
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
  %feature("customdoc:main", "  $brief\n\n$overview\n$main");
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


// Note: Only from 3.0.0 onwards,
// DirectorException inherits from std::exception
#if SWIG_VERSION >= 0x030000
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

#ifdef WITH_PYTHON3
// See https://github.com/casadi/casadi/issues/701
// Recent numpys will only catch TypeError or ValueError in printing logic
%exception __bool__ {
 try {
    $action
  } catch (const std::exception& e) {
   SWIG_exception(SWIG_TypeError, e.what());
  }
}
#else
%exception __nonzero__ {
 try {
    $action
  } catch (const std::exception& e) {
   SWIG_exception(SWIG_TypeError, e.what());
  }
}
#endif
#else
// Exceptions handling
%include "exception.i"
%exception {
  try {
    $action
   } catch(const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
   } catch (const Swig::DirectorException& e) {
    SWIG_exception(SWIG_TypeError, e.getMessage());
   }
}

// Python sometimes takes an approach to not check, but just try.
// It expects a python error to be thrown.
%exception __int__ {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  } catch (const Swig::DirectorException& e) {
    SWIG_exception(SWIG_TypeError, e.getMessage());
  }
}

#ifdef WITH_PYTHON3
// See https://github.com/casadi/casadi/issues/701
// Recent numpys will only catch TypeError or ValueError in printing logic
%exception __bool__ {
 try {
    $action
  } catch (const std::exception& e) {
   SWIG_exception(SWIG_TypeError, e.what());
  } catch (const Swig::DirectorException& e) {
    SWIG_exception(SWIG_TypeError, e.getMessage());
  }
}
#else
%exception __nonzero__ {
 try {
    $action
  } catch (const std::exception& e) {
   SWIG_exception(SWIG_TypeError, e.what());
  }
  catch (const Swig::DirectorException& e) {
    SWIG_exception(SWIG_TypeError, e.getMessage());
  }
}
#endif
#endif

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
    bool to_ptr(GUESTOBJECT *p, const std::vector<bool> **m);
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
      if (PyBool_Check(cr)) {
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
#ifdef WITH_PYTHON3
      return PyLong_FromLongLong(*a);
#else
      // For python on Windows
      if (*a > PyInt_GetMax() || *a < -(PyInt_GetMax()-1)) return PyLong_FromLongLong(*a);
      return PyInt_FromLong(*a);
#endif

#elif defined(SWIGMATLAB)
      return mxCreateDoubleScalar(static_cast<double>(*a));
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
#ifdef SWIGPYTHON

      // Some built-in types are iterable
      if (PyDict_Check(p) || PyString_Check(p) || PySet_Check(p) || PyUnicode_Check(p)) return false;

      // Make sure shape is 1D, if defined.
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

    template<typename M> GUESTOBJECT* from_ptr(const std::vector<M> *a) {
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
      if (to_generic<casadi_int>(p, m)
          || to_generic<double>(p, m)
          || to_generic<std::string>(p, m)
          || to_generic<std::vector<casadi_int> >(p, m)
          || to_generic<std::vector<double> >(p, m)
          || to_generic<std::vector<bool> >(p, m)
          || to_generic<std::vector<std::string> >(p, m)
          || to_generic<std::vector<std::vector<casadi_int> > >(p, m)
          || to_generic<std::vector<std::vector<double> > >(p, m)
          || to_generic<casadi::Function>(p, m)
          || to_generic<std::vector<casadi::Function> >(p, m)
          || to_generic<casadi::GenericType::Dict>(p, m)) {
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
      case OT_DICT: return from_tmp(a->as_dict());
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

#ifdef SWIGPYTHON
      if (PyString_Check(p) || PyUnicode_Check(p)) {
        if (m) (*m)->clear();
        char* my_char = SWIG_Python_str_AsChar(p);
        if (m) (*m)->append(my_char);
        SWIG_Python_str_DelForPy3(my_char);
        return true;
      }
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
      if (mxIsChar(p) && mxGetM(p)<=1) {
        if (m) {
          if (mxGetM(p)==0) return true;
          size_t len=mxGetN(p);
          std::vector<char> s(len+1);
          if (mxGetString(p, &s[0], (len+1)*sizeof(char))) return false;
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
      // Python slice
      if (PySlice_Check(p)) {
        PySliceObject *r = (PySliceObject*)(p);
        if (m) {
          (**m).start = (r->start == Py_None || PyNumber_AsSsize_t(r->start, NULL) <= std::numeric_limits<int>::min())
            ? std::numeric_limits<casadi_int>::min() : PyInt_AsLong(r->start);
          (**m).stop  = (r->stop ==Py_None || PyNumber_AsSsize_t(r->stop, NULL)>= std::numeric_limits<int>::max())
            ? std::numeric_limits<casadi_int>::max() : PyInt_AsLong(r->stop);
          if(r->step !=Py_None) (**m).step  = PyInt_AsLong(r->step);
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
#ifdef SWIGPYTHON
      if (PyDict_Check(p)) {
        PyObject *key, *value;
        Py_ssize_t pos = 0;
        while (PyDict_Next(p, &pos, &key, &value)) {
          if (!(PyString_Check(key) || PyUnicode_Check(key))) return false;
          if (m) {
            char* c_key = SWIG_Python_str_AsChar(key);
            M *v=&(**m)[std::string(c_key)], *v2=v;
            SWIG_Python_str_DelForPy3(c_key);
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

      // Try first converting to a temporary DM
      {
        DM tmp;
        if(to_val(p, m? &tmp: 0)) {
          if (m) **m = tmp;
          return true;
        }
      }

#ifdef SWIGPYTHON
      // Numpy arrays will be cast to dense SX
      if (SX_from_array(p, m)) return true;
      // Object has __SX__ method
      if (PyObject_HasAttrString(p,"__SX__")) {
        PyObject *cr = PyObject_CallMethod(p, (char*) "__SX__", 0);
        if (!cr) return false;
        casadi_int flag = to_ptr(cr, m);
        Py_DECREF(cr);
        return flag;
      }
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

      // Try first converting to a temporary DM
      {
        DM tmp;
        if(to_val(p, m ? &tmp : 0)) {
          if (m) **m = tmp;
          return true;
        }
      }

#ifdef SWIGPYTHON
      if (PyObject_HasAttrString(p,"__MX__")) {
        PyObject *cr = PyObject_CallMethod(p, (char*) "__MX__", 0);
        if (!cr) return false;
        casadi_int flag = to_ptr(cr, m);
        Py_DECREF(cr);
        return flag;
      }
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

      char* c_classname = SWIG_Python_str_AsChar(classname);
      bool ret = strcmp(c_classname, name)==0;

      Py_DECREF(classo);Py_DECREF(classname);
      SWIG_Python_str_DelForPy3(c_classname);
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
      // Object has __DM__ method
      if (PyObject_HasAttrString(p,"__DM__")) {
        char name[] = "__DM__";
        PyObject *cr = PyObject_CallMethod(p, name, 0);
        if (!cr) return false;
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

 // Define all input typemaps
%define %casadi_input_typemaps(xName, xPrec, xType...)
 // Pass input by value, check if matches
%typemap(typecheck, noblock=1, precedence=xPrec, fragment="casadi_all") xType {
  $1 = casadi::to_ptr($input, static_cast< xType **>(0));
 }

 // Directorout typemap; as input by value
%typemap(directorout, noblock=1, fragment="casadi_all") xType {
    if (!casadi::to_val($input, &$result)) {
      %dirout_fail(SWIG_TypeError,"$type");
    }
 }

 // Pass input by value, convert argument
%typemap(in, doc=xName, noblock=1, fragment="casadi_all") xType {
  if (!casadi::to_val($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input $argnum to type '" xName "'.");
 }

 // Pass input by value, cleanup
%typemap(freearg, noblock=1) xType {}

 // Pass input by reference, check if matches
%typemap(typecheck, noblock=1, precedence=xPrec, fragment="casadi_all") const xType& {
  $1 = casadi::to_ptr($input, static_cast< xType **>(0));
 }

 // Pass input by reference, convert argument
%typemap(in, doc=xName, noblock=1, fragment="casadi_all") const xType & (xType m) {
  $1 = &m;
  if (!casadi::to_ptr($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input $argnum to type '" xName "'.");
 }

 // Pass input by reference, cleanup
%typemap(freearg, noblock=1) const xType & {}

%enddef

 // Define all output typemaps
%define %casadi_output_typemaps(xName, xType...)

 // Return-by-value
%typemap(out, doc=xName, noblock=1, fragment="casadi_all") xType, const xType {
  if(!($result = casadi::from_ref($1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type '" xName "'.");
}

// Return a const-ref behaves like return-by-value
%typemap(out, doc=xName, noblock=1, fragment="casadi_all") const xType& {
  if(!($result = casadi::from_ptr($1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type '" xName "'.");
}

// Inputs marked OUTPUT are also returned by the function, ...
%typemap(argout, noblock=1,fragment="casadi_all") xType &OUTPUT {
  %append_output(casadi::from_ptr($1));
 }

// ... and the corresponding inputs are ignored
%typemap(in, doc=xName, noblock=1, numinputs=0) xType &OUTPUT (xType m) {
 $1 = &m;
}

 // Directorin typemap; as output
%typemap(directorin, noblock=1, fragment="casadi_all") xType, const xType {
    if(!($input = casadi::from_ref($1))) %dirout_fail(SWIG_TypeError,"For director inputs, failed to convert input to " xName ".");
 }

 // Directorin typemap; as output
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
%typemap(in, doc=xName, noblock=1, fragment="casadi_all") xType &INOUT (xType m) {
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

 // Define all typemaps for a template instantiation without proxy classes
%define %casadi_template(xName, xPrec, xType...)
%template() xType;
%casadi_input_typemaps(xName, xPrec, xType)
%casadi_output_typemaps(xName, %arg(xType))
%enddef

 // Define all input and ouput typemaps
%define %casadi_typemaps(xName, xPrec, xType...)
%casadi_input_typemaps(xName, xPrec, xType)
%casadi_output_typemaps(xName, xType)
%enddef

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
%typemap(in, doc="memoryview(ro)", noblock=1, fragment="casadi_all") (const double * a, casadi_int size) (Py_buffer* buffer) {
  if (!PyMemoryView_Check($input)) SWIG_exception_fail(SWIG_TypeError, "Must supply a MemoryView.");
  buffer = PyMemoryView_GET_BUFFER($input);
  $1 = static_cast<double*>(buffer->buf); // const double cast comes later
  $2 = buffer->len;
 }

%typemap(in, doc="memoryview(rw)", noblock=1, fragment="casadi_all") (double * a, casadi_int size)  (Py_buffer* buffer) {
  if (!PyMemoryView_Check($input)) SWIG_exception_fail(SWIG_TypeError, "Must supply a writable MemoryView.");
  buffer = PyMemoryView_GET_BUFFER($input);
  if (buffer->readonly) SWIG_exception_fail(SWIG_TypeError, "Must supply a writable MemoryView.");
  $1 = static_cast<double*>(buffer->buf);
  $2 = buffer->len;
 }

// Directorin typemap; as output
%typemap(directorin, noblock=1, fragment="casadi_all") (const double** arg, const std::vector<casadi_int>& sizes_arg) (PyObject* my_tuple) {
  PyObject * arg_tuple = PyTuple_New($2.size());
  for (casadi_int i=0;i<$2.size();++i) {
    
#ifdef WITH_PYTHON3
    PyObject* buf = $1[i] ? PyMemoryView_FromMemory(reinterpret_cast<char*>(const_cast<double*>($1[i])), $2[i]*sizeof(double), PyBUF_READ) : SWIG_Py_Void();
#else
    PyObject* buf = $1[i] ? PyBuffer_FromMemory(const_cast<double*>($1[i]), $2[i]*sizeof(double)) : SWIG_Py_Void();
#endif
    PyTuple_SET_ITEM(arg_tuple, i, buf);
  }
  $input = arg_tuple;
}

%typemap(directorin, noblock=1, fragment="casadi_all") (double** res, const std::vector<casadi_int>& sizes_res) {
  PyObject* res_tuple = PyTuple_New($2.size());
  for (casadi_int i=0;i<$2.size();++i) {
#ifdef WITH_PYTHON3
    PyObject* buf = $1[i] ? PyMemoryView_FromMemory(reinterpret_cast<char*>(const_cast<double*>($1[i])), $2[i]*sizeof(double), PyBUF_WRITE) : SWIG_Py_Void();
#else
    PyObject* buf = $1[i] ? PyBuffer_FromReadWriteMemory($1[i], $2[i]*sizeof(double)) : SWIG_Py_Void();
#endif
    PyTuple_SET_ITEM(res_tuple, i, buf);
  }
  $input = res_tuple;
}

%typemap(in, doc="void*", noblock=1, fragment="casadi_all") void* raw {
  $1 = PyCapsule_GetPointer($input, NULL);
}

%typemap(out, doc="void*", noblock=1, fragment="casadi_all") void* {
  $result = PyCapsule_New($1, NULL,NULL);
}
#endif

%casadi_typemaps(L_STR, PREC_STRING, std::string)
%casadi_template(LL L_STR LR, PREC_VECTOR, std::vector<std::string>)
%casadi_typemaps("Sparsity", PREC_SPARSITY, casadi::Sparsity)
%casadi_template(LL "Sparsity" LR, PREC_SPARSITY, std::vector< casadi::Sparsity>)
%casadi_template(LL LL "Sparsity"  LR  LR, PREC_SPARSITY, std::vector<std::vector< casadi::Sparsity> >)
%casadi_template(LDICT("Sparsity"), PREC_SPARSITY, std::map<std::string, casadi::Sparsity >)
%casadi_template(LDICT(LL "Sparsity" LR), PREC_SPARSITY, std::map<std::string, std::vector<casadi::Sparsity > >)
%casadi_template(LPAIR(LDICT("Sparsity"),"[" L_STR "]"), PREC_SPARSITY, std::pair<std::map<std::string, casadi::Sparsity >, std::vector<std::string> >)
%casadi_typemaps(L_BOOL, SWIG_TYPECHECK_BOOL, bool)
%casadi_template("[" L_BOOL "]", SWIG_TYPECHECK_BOOL, std::vector<bool>)
%casadi_template("[[" L_BOOL "]]", SWIG_TYPECHECK_BOOL, std::vector<std::vector<bool> >)
%casadi_typemaps( L_INT , SWIG_TYPECHECK_INTEGER, casadi_int)

#ifdef MATLABSTYLE
#define LABEL "[int,int]"
#else
#define LABEL LPAIR("int","int")
#endif
%casadi_template(LABEL, SWIG_TYPECHECK_INTEGER, std::pair<casadi_int,casadi_int>)
#undef LABEL
%casadi_template("[" L_INT "]", PREC_IVector, std::vector<casadi_int>)
%casadi_template(LL "[" L_INT "]" LR, PREC_IVectorVector, std::vector<std::vector<casadi_int> >)
%casadi_typemaps(L_DOUBLE, SWIG_TYPECHECK_DOUBLE, double)
%casadi_template("[" L_DOUBLE "]", SWIG_TYPECHECK_DOUBLE, std::vector<double>)
%casadi_template(LL "[" L_DOUBLE "]" LR, SWIG_TYPECHECK_DOUBLE, std::vector<std::vector<double> >)
%casadi_typemaps("SXElem", PREC_SX, casadi::SXElem)
%casadi_template(LL "SXElem" LR, PREC_SXVector, std::vector<casadi::SXElem>)
%casadi_typemaps("SX", PREC_SX, casadi::Matrix<casadi::SXElem>)
%casadi_template(LL "SX" LR, PREC_SXVector, std::vector< casadi::Matrix<casadi::SXElem> >)
%casadi_template(LL LL "SX" LR LR, PREC_SXVectorVector, std::vector<std::vector< casadi::Matrix<casadi::SXElem> > >)
%casadi_template(LDICT("SX"), PREC_SX, std::map<std::string, casadi::Matrix<casadi::SXElem> >)
%casadi_typemaps("MX", PREC_MX, casadi::MX)
%casadi_template(LL "MX" LR, PREC_MXVector, std::vector<casadi::MX>)
%casadi_template(LL LL "MX" LR LR, PREC_MXVectorVector, std::vector<std::vector<casadi::MX> >)
%casadi_template(LDICT("MX"), PREC_MX, std::map<std::string, casadi::MX>)
%casadi_template(LPAIR("MX","MX"), PREC_MXVector, std::pair<casadi::MX, casadi::MX>)
%casadi_typemaps("DM", PREC_DM, casadi::Matrix<double>)
%casadi_template(LL "DM" LR, PREC_DMVector, std::vector< casadi::Matrix<double> >)
%casadi_template(LL LL "DM" LR LR, PREC_DMVectorVector, std::vector<std::vector< casadi::Matrix<double> > >)
%casadi_template(LDICT("DM"), PREC_DM, std::map<std::string, casadi::Matrix<double> >)
%casadi_typemaps("IM", PREC_IM, casadi::Matrix<casadi_int>)
// Without CASADI_INT_TYPE, you get SwigValueWrapper
// With it, docstrings are screwed
%casadi_typemaps("GenericType", PREC_GENERICTYPE, casadi::GenericType)
%casadi_template(LL "GenericType" LR, PREC_GENERICTYPE, std::vector<casadi::GenericType>)
%casadi_typemaps("Slice", PREC_SLICE, casadi::Slice)
%casadi_typemaps("Function", PREC_FUNCTION, casadi::Function)
%casadi_template(LL "Function" LR, PREC_FUNCTION, std::vector<casadi::Function>)
%casadi_template(LPAIR("Function","Function"), PREC_FUNCTION, std::pair<casadi::Function, casadi::Function>)
%casadi_template(L_DICT, PREC_DICT, std::map<std::string, casadi::GenericType>)
%casadi_template(LDICT(LL L_STR LR), PREC_DICT, std::map<std::string, std::vector<std::string> >)

#undef L_INT
#undef L_BOOL
#undef LPAIR
#undef L_DOUBLE
#undef L_DICT
#undef LL
#undef LR
#undef L_STR
#undef MATLABSTYLE

// Matlab is index-1 based
#ifdef SWIGMATLAB
%typemap(in, doc="index", noblock=1) casadi_index {
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

%pythoncode%{
try:
  from numpy import pi, inf
except:
  pass

arcsin = lambda x: _casadi.asin(x)
arccos = lambda x: _casadi.acos(x)
arctan = lambda x: _casadi.atan(x)
arctan2 = lambda x,y: _casadi.atan2(x, y)
arctanh = lambda x: _casadi.atanh(x)
arcsinh = lambda x: _casadi.asinh(x)
arccosh = lambda x: _casadi.acosh(x)
%}
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
def veccat(*args): return _veccat(args)
def vertcat(*args): return _vertcat(args)
def horzcat(*args): return _horzcat(args)
def diagcat(*args): return _diagcat(args)
def vvcat(args): return _veccat(args)
def vcat(args): return _vertcat(args)
def hcat(args): return _horzcat(args)
def dcat(args): return _diagcat(args)
%}

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
#endif // SWIGMATLAB

#ifdef WITH_PYTHON3
%rename(__bool__) __nonzero__;
#endif

#ifdef SWIGPYTHON

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
                      """)

    def __setitem__(self,s,val):
          if isinstance(s,tuple) and len(s)==2:
            return self.set(val, False, s[0], s[1])
          return self.set(val, False, s)

    @property
    def nz(self):
      return NZproxy(self)

%}
%enddef

%define %python_array_wrappers(arraypriority)
%pythoncode %{

  __array_priority__ = arraypriority

  def __array_wrap__(self,out_arr,context=None):
    if context is None:
      return out_arr
    name = context[0].__name__
    args = list(context[1])

    if len(context[1])==3:
      raise Exception("Error with %s. Looks like you are using an assignment operator, such as 'a+=b' where 'a' is a numpy type. This is not supported, and cannot be supported without changing numpy." % name)

    if "vectorized" in name:
        name = name[:-len(" (vectorized)")]

    conversion = {"multiply": "mul", "divide": "div", "true_divide": "div", "subtract":"sub","power":"pow","greater_equal":"ge","less_equal": "le", "less": "lt", "greater": "gt"}
    if name in conversion:
      name = conversion[name]
    if len(context[1])==2 and context[1][1] is self and not(context[1][0] is self):
      name = 'r' + name
      args.reverse()
    if not(hasattr(self,name)) or ('mul' in name):
      name = '__' + name + '__'
    fun=getattr(self, name)
    return fun(*args[1:])

  def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
    conversion = {"multiply": "mul", "divide": "div", "true_divide": "div", "subtract":"sub","power":"pow","greater_equal":"ge","less_equal": "le", "less": "lt", "greater": "gt"}
    name = ufunc.__name__
    inputs = list(inputs)
    if len(inputs)==3:
      import warnings
      warnings.warn("Error with %s. Looks like you are using an assignment operator, such as 'a+=b' where 'a' is a numpy type. This is not supported, and cannot be supported without changing numpy." % name, RuntimeWarning)
      return NotImplemented
    if "vectorized" in name:
        name = name[:-len(" (vectorized)")]
    if name in conversion:
      name = conversion[name]
    if len(inputs)==2 and inputs[1] is self and not(inputs[0] is self):
      name = 'r' + name
      inputs.reverse()
    if not(hasattr(self,name)) or ('mul' in name):
      name = '__' + name + '__'
    try:
      assert method=="__call__"
      fun=getattr(self, name)
      return fun(*inputs[1:])
    except:
      # Fall back to numpy conversion
      new_inputs = list(inputs)
      try:
        new_inputs[0] = new_inputs[0].full()
      except:
        import warnings
        warnings.warn("Implicit conversion of symbolic CasADi type to numeric matrix not supported.\n"
                               + "This may occur when you pass a CasADi object to a numpy function.\n"
                               + "Use an equivalent CasADi function instead of that numpy function.", RuntimeWarning)
        return NotImplemented
      return new_inputs[0].__array_ufunc__(ufunc, method, *new_inputs, **kwargs)


  def __array__(self,*args,**kwargs):
    import numpy as n
    if len(args) > 1 and isinstance(args[1],tuple) and isinstance(args[1][0],n.ufunc) and isinstance(args[1][0],n.ufunc) and len(args[1])>1 and args[1][0].nin==len(args[1][1]):
      if len(args[1][1])==3:
        raise Exception("Error with %s. Looks like you are using an assignment operator, such as 'a+=b'. This is not supported when 'a' is a numpy type, and cannot be supported without changing numpy itself. Either upgrade a to a CasADi type first, or use 'a = a + b'. " % args[1][0].__name__)
      return n.array([n.nan])
    else:
      if hasattr(self,'__array_custom__'):
        return self.__array_custom__(*args,**kwargs)
      else:
        try:
          return self.full()
        except:
          raise Exception("Implicit conversion of symbolic CasADi type to numeric matrix not supported.\n"
                     + "This may occur when you pass a CasADi object to a numpy function.\n"
                     + "Use an equivalent CasADi function instead of that numpy function.")

%}
%enddef
#endif // SWIGPYTHON

#ifdef SWIGXML
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

%include <casadi/core/printable.hpp>

namespace casadi{
%extend PrintableCommon {
#ifdef SWIGPYTHON
  %pythoncode %{
    def __str__(self): return self.str()
    def repr(self): return self.type_name() + '(' + self.str() + ')'
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
 DECL M casadi_mtimes(const M& x, const M& y) {
 return mtimes(x, y);
 }
 DECL M casadi_mtimes(const std::vector< M > &args) {
 return mtimes(args);
 }
 DECL M casadi_mac(const M& X, const M& Y, const M& Z) {
 return mac(X, Y, Z);
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
    #endif
  %enddef
#else
  %define SPARSITY_INTERFACE_FUN(DECL, FLAG, M)
    SPARSITY_INTERFACE_FUN_BASE(DECL, FLAG, M)
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

DECL M casadi_bilin(const M& A, const M& x, const M& y) {
  return bilin(A, x, y);
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

DECL M casadi_jtimes(const M& ex, const M& arg, const M& v, bool tr=false) {
  return jtimes(ex, arg, v, tr);
}

DECL M casadi_linearize(const M& f, const M& x, const M& x0) {
  return linearize(f, x, x0);
}

DECL std::vector<bool> casadi_which_depends(const M& expr, const M& var,
                                            casadi_int order=1, bool tr=false) {
  return which_depends(expr, var, order, tr);
}

DECL bool casadi_is_linear(const M& expr, const M& var) {
  return is_linear(expr, var);
}

DECL bool casadi_is_quadratic(const M& expr, const M& var) {
  return is_quadratic(expr, var);
}

DECL M casadi_gradient(const M &ex, const M &arg) {
  return gradient(ex, arg);
}

DECL M casadi_tangent(const M &ex, const M &arg) {
  return tangent(ex, arg);
}

DECL M casadi_hessian(const M& ex, const M& arg, M& OUTPUT1) {
  return hessian(ex, arg, OUTPUT1);
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
#endif // FLAG & IS_MEMBER

#if FLAG & IS_GLOBAL
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

DECL void casadi_shared(const std::vector< M >& ex,
                               std::vector< M >& OUTPUT1,
                               std::vector< M >& OUTPUT2,
                               std::vector< M >& OUTPUT3,
                               const std::string& v_prefix="v_",
                               const std::string& v_suffix="") {
  shared(ex, OUTPUT1, OUTPUT2, OUTPUT3, v_prefix, v_suffix);
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
DECL M casadi_floor(const M& x) { return floor(x); }
DECL M casadi_ceil(const M& x) { return ceil(x); }
DECL M casadi_erf(const M& x) { return erf(x); }
DECL M casadi_erfinv(const M& x) { using casadi::erfinv; return erfinv(x); }
DECL M casadi_sign(const M& x) { using casadi::sign; return sign(x); }
DECL M casadi_power(const M& x, const M& n) { return pow(x, n); }
DECL M casadi_mod(const M& x, const M& y) { return fmod(x, y); }
DECL M casadi_fmod(const M& x, const M& y) { return fmod(x, y); }
DECL M casadi_atan2(const M& x, const M& y) { return atan2(x, y); }
DECL M casadi_fmin(const M& x, const M& y) { return fmin(x, y); }
DECL M casadi_fmax(const M& x, const M& y) { return fmax(x, y); }
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
DECL M casadi_convexify(const M& H,
        const Dict& opts = Dict()) {
  return convexify(H, opts);
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
    void assign(const casadi::Matrix<double>&rhs) { (*$self)=rhs; }
    %matrix_helpers(casadi::Matrix<double>)

  }

}

// Extend DM with SWIG unique features
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


#ifdef SWIGPYTHON
namespace casadi{
%extend Matrix<double> {

%python_array_wrappers(999.0)

// The following code has some trickery to fool numpy ufunc.
// Normally, because of the presence of __array__, an ufunctor like nump.sqrt
// will unleash its activity on the output of __array__
// However, we wish DM to remain a DM
// So when we receive a call from a functor, we return a dummy empty array
// and return the real result during the postprocessing (__array_wrap__) of the functor.
%pythoncode %{
  def __array_custom__(self,*args,**kwargs):
    if "dtype" in kwargs and not(isinstance(kwargs["dtype"],n.double)):
      return n.array(self.full(),dtype=kwargs["dtype"])
    else:
      return self.full()
%}

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
%}


#ifdef WITH_PYTHON3
%pythoncode %{
  def __bool__(self):
    if self.numel()!=1:
      raise Exception("Only a scalar can be cast to a float")
    if self.nnz()==0:
      return False
    return float(self)!=0
%}
#else
%pythoncode %{
  def __nonzero__(self):
    if self.numel()!=1:
      raise Exception("Only a scalar can be cast to a float")
    if self.nnz()==0:
      return False
    return float(self)!=0
%}
#endif

%pythoncode %{
  def __abs__(self):
    return abs(float(self))
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
    function out = sum(self,varargin)
      narginchk(1,3);
      if nargin==1
        if is_vector(self)
          if is_column(self)
            out = sum1(self);
          else
            out = sum2(self);
          end
        else
          out = sum1(self);
        end
      else
        i = varargin{1};
        if i==1
          out = sum1(self);
        elseif i==2
          out = sum2(self);
        else
          error('sum argument (if present) must be 1 or 2');
        end
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
#ifdef WITH_PYTHON3
namespace casadi {
  %extend GenericExpressionCommon {
    %pythoncode %{
      def __hash__(self):
        try:
          return self.element_hash()
        except:
          return SharedObject.__hash__(self)
      def __matmul__(x, y): return _casadi.mtimes(x, y)
      def __rmatmul__(x, y): return _casadi.mtimes(y, x)
    %}
  }
}
#endif
namespace casadi {
  %extend GenericExpressionCommon {
    %pythoncode %{
      def __add__(x, y): return _casadi.plus(x, y)
      def __radd__(x, y): return _casadi.plus(y, x)
      def __sub__(x, y): return _casadi.minus(x, y)
      def __rsub__(x, y): return _casadi.minus(y, x)
      def __mul__(x, y): return _casadi.times(x, y)
      def __rmul__(x, y): return _casadi.times(y, x)
      def __div__(x, y): return _casadi.rdivide(x, y)
      def __rdiv__(x, y): return _casadi.rdivide(y, x)
      def __truediv__(x, y): return _casadi.rdivide(x, y)
      def __rtruediv__(x, y): return _casadi.rdivide(y, x)
      def __lt__(x, y): return _casadi.lt(x, y)
      def __rlt__(x, y): return _casadi.lt(y, x)
      def __le__(x, y): return _casadi.le(x, y)
      def __rle__(x, y): return _casadi.le(y, x)
      def __gt__(x, y): return _casadi.lt(y, x)
      def __rgt__(x, y): return _casadi.lt(x, y)
      def __ge__(x, y): return _casadi.le(y, x)
      def __rge__(x, y): return _casadi.le(x, y)
      def __eq__(x, y): return _casadi.eq(x, y)
      def __req__(x, y): return _casadi.eq(y, x)
      def __ne__(x, y): return _casadi.ne(x, y)
      def __rne__(x, y): return _casadi.ne(y, x)
      def __pow__(x, n): return _casadi.power(x, n)
      def __rpow__(n, x): return _casadi.power(x, n)
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
      def floor(x): return _casadi.floor(x)
      def ceil(x): return _casadi.ceil(x)
      def erf(x): return _casadi.erf(x)
      def sign(x): return _casadi.sign(x)
      def fmod(x, y): return _casadi.mod(x, y)
      def __copysign__(x, y): return _casadi.copysign(x, y)
      def __rcopysign__(y, x): return _casadi.copysign(x, y)
      def copysign(x, y): return _casadi.copysign(x, y)
      def rcopysign(y, x): return _casadi.copysign(x, y)
      def __constpow__(x, y): return _casadi.constpow(x, y)
      def __rconstpow__(y, x): return _casadi.constpow(x, y)
      def constpow(x, y): return _casadi.constpow(x, y)
      def rconstpow(y, x): return _casadi.constpow(x, y)
    %}
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
  }

} // namespace casadi
#endif // SWIGPYTHON

%feature("director") casadi::Callback;

%include <casadi/core/importer.hpp>
%include <casadi/core/callback.hpp>
%include <casadi/core/global_options.hpp>
%include <casadi/core/casadi_meta.hpp>
%include <casadi/core/integration_tools.hpp>
%include <casadi/core/nlp_tools.hpp>
%include <casadi/core/nlp_builder.hpp>
%include <casadi/core/variable.hpp>
%include <casadi/core/dae_builder.hpp>
%include <casadi/core/xml_file.hpp>

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
  %pythoncode %{
    def unpack(self):
      type = SerializerBase.type_to_string(self._pop_type())
      f = getattr(self, "blind_unpack_"+type)
      return f()
  %}
}
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

%feature("director") casadi::OptiCallback;

// Return-by-value
%typemap(out, doc="double", noblock=1, fragment="casadi_all") casadi::native_DM {
  if(!($result = full_or_sparse($1, true))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type 'double'.");
}


%apply casadi_int &OUTPUT { Opti::ConstraintType &OUTPUT };

%typemap(argout, noblock=1,fragment="casadi_all") casadi::Opti::ConstraintType &OUTPUT {
  %append_output(casadi::from_ptr((casadi_int *) $1));
}

%typemap(in, doc="Opti.ConstraintType", noblock=1, numinputs=0) casadi::Opti::ConstraintType &OUTPUT (casadi::Opti::ConstraintType m) {
 $1 = &m;
}


#ifdef SWIGPYTHON

%define make_property(class, name)
  %rename(_ ## name) class ## :: ## name;
  %extend class {
    %pythoncode %{
      @property
      def name(self):
        return self._ ## name()
    %}
  }
%enddef


make_property(casadi::Opti, debug);
make_property(casadi::Opti, advanced);
make_property(casadi::OptiSol, opti);

%define make_property_opti(name)
  make_property(casadi::Opti, name);
%enddef

make_property(casadi::OptiSol, debug);
make_property_opti(f)
make_property_opti(g)
make_property_opti(x)
make_property_opti(p)
make_property_opti(lam_g)
make_property_opti(lbg)
make_property_opti(ubg)
make_property_opti(nx)
make_property_opti(np)
make_property_opti(ng)

make_property(casadi::Opti, casadi_solver);
%define opti_metadata_modifiers(class)
  %rename(_variable) class ## :: variable;
  %rename(_parameter) class ## :: parameter;
  %rename(_subject_to) class ## :: subject_to;
  %extend class {
    %pythoncode %{
      def parameter(self,*args):
        import sys
        import os
        frame = sys._getframe(1)
        meta = {"stacktrace": {"file":os.path.abspath(frame.f_code.co_filename),"line":frame.f_lineno,"name":frame.f_code.co_name}}
        ret = self._parameter(*args)
        self.update_user_dict(ret, meta)
        return ret

      def variable(self,*args):
        import sys
        import os
        frame = sys._getframe(1)
        meta = {"stacktrace": {"file":os.path.abspath(frame.f_code.co_filename),"line":frame.f_lineno,"name":frame.f_code.co_name}}
        ret = self._variable(*args)
        self.update_user_dict(ret, meta)
        return ret

      def subject_to(self,*args):
        if len(args)==0:
          return self._subject_to()
        import sys
        import os
        frame = sys._getframe(1)
        meta = {"stacktrace": {"file":os.path.abspath(frame.f_code.co_filename),"line":frame.f_lineno,"name":frame.f_code.co_name}}
        ret = self._subject_to(*args)
        self.update_user_dict(args[0], meta)
        return ret
    %}
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
          meta = struct('stacktrace', st(1));
        else
          meta = struct;
        end
        out = self.internal_variable(varargin{:});
        self.update_user_dict(out, meta);
      end
      function out = parameter(self, varargin)
        st = dbstack('-completenames',1);
        if length(st)>0
          meta = struct('stacktrace', st(1));
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
          meta = struct('stacktrace', st(1));
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

// Cleanup for dependent modules
%exception {
  $action
}
