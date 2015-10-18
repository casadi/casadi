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

 // Include all public CasADi C++
%{
#include <casadi/casadi.hpp>
#include <casadi/core/casadi_interrupt.hpp>
%}

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
      return PyErr_CheckSignals();
    }


  }

%}
%init %{
  // Set logger functions
  casadi::Logger::writeWarn = casadi::pythonlogger;
  casadi::Logger::writeProg = casadi::pythonlogger;
  casadi::Logger::writeDebug = casadi::pythonlogger;
  casadi::Logger::writeAll = casadi::pythonlogger;

  // @jgillis: please document
  casadi::InterruptHandler::checkInterrupted = casadi::pythoncheckinterrupted;
%}
#elif defined(SWIGMATLAB)
%{
  namespace casadi {
    // Redirect printout to mexPrintf
    static void mexlogger(const char* s, std::streamsize num, bool error) {
      mexPrintf("%.*s", static_cast<int>(num), s);
    }

    // Flush the command window buffer (needed in gui mode)
    static void mexflush(bool error) {
      mexEvalString("drawnow('update');");
    }

    // Undocumented matlab feature
    extern "C" bool utIsInterruptPending();

    static bool mexcheckinterrupted() {
      return utIsInterruptPending();
    }
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

  // Set library path
  casadi::CasadiOptions::setCasadiPath(path);

  // @jgillis: please document
  mxArray *warning_rhs[] = {mxCreateString("error"),

                            mxCreateString("SWIG:OverloadError")};
  mexCallMATLAB(0, 0, 2, warning_rhs, "warning");
  mxDestroyArray(warning_rhs[0]);
  mxDestroyArray(warning_rhs[1]);

  
  // Set logger functions
  casadi::Logger::writeWarn = casadi::mexlogger;
  casadi::Logger::writeProg = casadi::mexlogger;
  casadi::Logger::writeDebug = casadi::mexlogger;
  casadi::Logger::writeAll = casadi::mexlogger;
  casadi::Logger::flush = casadi::mexflush;

  // @jgillis: please document
  casadi::InterruptHandler::checkInterrupted = casadi::mexcheckinterrupted;
%}
#endif

// Turn off the warnings that certain methods are effectively ignored, this seams to be a false warning,
// for example vertcat(SXVector), vertcat(DMatrixVector) and vertcat(MXVector) appears to work fine
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

#ifdef SWIGMATLAB
%rename(disp) repr;
#else
%ignore print;
%ignore repr;
#endif

%begin %{
#define SWIG_PYTHON_OUTPUT_TUPLE
%}

// Print representation
#ifdef SWIGMATLAB
#define SWIG_REPR disp
#else
#define SWIG_REPR __repr__
#endif

// Print description
#ifdef SWIGMATLAB
#define SWIG_STR print
#else
#define SWIG_STR __str__
#endif


//#endif // SWIGPYTHON


#ifdef SWIGPYTHON
%pythoncode %{

import contextlib

@contextlib.contextmanager
def internalAPI():
    backup = CasadiOptions.getAllowedInternalAPI()
    CasadiOptions.setAllowedInternalAPI(True)
    try:
      yield
    finally:
      CasadiOptions.setAllowedInternalAPI(backup)

class _copyableObject(_object):
  def __copy__(self):
    return self.__class__(self)

  def __deepcopy__(self,dummy=None):
    return self.__class__(self)

_object = _copyableObject

_swig_repr_default = _swig_repr
def _swig_repr(self):
  if hasattr(self,'getRepresentation'):
    return self.getRepresentation()
  else:
    return _swig_repr_default(self)

%}
#endif // WITH_SWIGPYTHON

#if defined(SWIGPYTHON) || defined(SWIGMATLAB)
%include "doc_merged.i"
#else
%include "doc.i"
#endif

%feature("autodoc", "1");

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

%define DEPRECATED_MSG(MSG)
if (deprecated("$decl",MSG)) SWIG_fail;
%enddef

%define INTERNAL_MSG()
if (internal("$decl")) SWIG_fail;
%enddef

#ifdef SWIGPYTHON
%wrapper %{
int deprecated(const std::string & c,const std::string & a) {
  std::string msg = "This CasADi function (" + c + ") is deprecated. " + a;
  return PyErr_WarnEx(PyExc_DeprecationWarning,msg.c_str(),3);
}
int internal(const std::string & c) {
  if (casadi::CasadiOptions::allowed_internal_api) return 0;
  std::string msg = "This CasADi function (" + c + ") is not part of the public API. Use at your own risk.";
  return PyErr_WarnEx(PyExc_SyntaxWarning,msg.c_str(),3);
}
%}
#endif // SWIGPYTHON

#ifdef SWIGMATLAB
%wrapper %{
int deprecated(const std::string & c,const std::string & a) {
  std::string msg = "This CasADi function (" + c + ") is deprecated. " + a;
  mexWarnMsgIdAndTxt("SWIG:DeprecationWarning",msg.c_str());
  return 0;
}
int internal(const std::string & c) {
  if (casadi::CasadiOptions::allowed_internal_api) return 0;
  std::string msg = "This CasADi function (" + c + ") is not part of the public API. Use at your own risk.";
  mexWarnMsgIdAndTxt("SWIG:SyntaxWarning",msg.c_str());
  return 0;
}
%}
#endif // SWIGMATLAB

#ifndef SWIGXML
%{
#define CATCH_OR_NOT(...) \
if (casadi::CasadiOptions::catch_errors_swig) { \
  try { \
    __VA_ARGS__ \
  } catch(const std::exception& e) { \
    SWIG_exception(SWIG_RuntimeError, e.what()); \
  } \
} else { \
  __VA_ARGS__ \
}

%}
#endif

// Exceptions handling
%include "exception.i"
%exception {
  CATCH_OR_NOT($action)
}

// Python sometimes takes an approach to not check, but just try.
// It expects a python error to be thrown.
%exception __int__ {
 try {
    $action
  } catch (const std::exception& e) { \
  SWIG_exception(SWIG_RuntimeError, e.what()); \
  }
}

// See https://github.com/casadi/casadi/issues/701
// Recent numpys will only catch TypeError or ValueError in printing logic
%exception __nonzero__ {
 try {
    $action
    // foobar
  } catch (const std::exception& e) { \
  SWIG_exception(SWIG_TypeError, e.what()); \
  }
}

#ifdef SWIGPYTHON
%feature("director:except") {
	if ($error != NULL) {
    SWIG_PYTHON_THREAD_BEGIN_BLOCK;
    PyErr_Print();
    SWIG_PYTHON_THREAD_END_BLOCK; 
		throw Swig::DirectorMethodException();
	}
}
#endif //SWIGPYTHON

%include "internal.i"
%include "deprecated.i"

#ifdef SWIGPYTHON

%{
#define SWIG_FILE_WITH_INIT
#include "numpy.hpp"
#define SWIG_PYTHON_CAST_MODE 1
%}

%init %{
import_array();
%}

#endif // SWIGPYTHON

%{
#define SWIG_Error_return(code, msg)  { std::cerr << "Error occured in CasADi SWIG interface code:" << std::endl << "  "<< msg << std::endl;SWIG_Error(code, msg); return 0; }
%}

#ifndef SWIGXML

%fragment("casadi_decl", "header") {
  namespace casadi {
    /* Check if Null or None */
    bool is_null(GUESTOBJECT *p);

    /* Convert a pointer in interfaced language to C++
     * Input: GUESTOBJECT pointer p
     * Output: Pointer to pointer: At input, pointer to pointer to temporary
     * The routine will either:
     *   - Do nothing, if 0
     *   - Change the pointer
     *   - Change the temporary object
     * Returns true upon success, else false
     */
    bool to_ptr(GUESTOBJECT *p, bool** m);
    bool to_ptr(GUESTOBJECT *p, int** m);
    bool to_ptr(GUESTOBJECT *p, double** m);
    bool to_ptr(GUESTOBJECT *p, std::string** m);
    bool to_ptr(GUESTOBJECT *p, casadi::Slice** m);
    bool to_ptr(GUESTOBJECT *p, casadi::Sparsity** m);
    bool to_ptr(GUESTOBJECT *p, casadi::DMatrix** m);
    bool to_ptr(GUESTOBJECT *p, casadi::IMatrix** m);
    bool to_ptr(GUESTOBJECT *p, casadi::SX** m);
    bool to_ptr(GUESTOBJECT *p, casadi::MX** m);
    bool to_ptr(GUESTOBJECT *p, casadi::Function** m);
    bool to_ptr(GUESTOBJECT *p, casadi::GenericType** m);
#ifdef SWIGMATLAB
    bool to_ptr(GUESTOBJECT *p, std::pair<int, int>** m);
#endif // SWIGMATLAB
    template<typename M1, typename M2> bool to_ptr(GUESTOBJECT *p, std::pair<M1, M2>** m);
#ifdef SWIGMATLAB
    bool to_ptr(GUESTOBJECT *p, std::vector<std::string>** m);
#endif // SWIGMATLAB
    template<typename M> bool to_ptr(GUESTOBJECT *p, std::vector<M>** m);
    template<typename M> bool to_ptr(GUESTOBJECT *p, std::map<std::string, M>** m);

    // Same as the above, but with pointer instead of pointer to pointer
    template<typename M> bool to_val(GUESTOBJECT *p, M* m);

    // Check if conversion is possible
    template<typename M> bool can_convert(GUESTOBJECT *p) { return to_ptr(p, static_cast<M**>(0));}

    // Assign to a vector, if conversion is allowed
    template<typename E, typename M> bool assign_vector(E* d, int sz, std::vector<M>** m);

    /* Convert result from CasADi to interfaced language */
    GUESTOBJECT* from_ptr(const casadi::GenericType *a);
    GUESTOBJECT* from_ptr(const bool *a);
    GUESTOBJECT* from_ptr(const int *a);
    GUESTOBJECT* from_ptr(const double *a);
    GUESTOBJECT* from_ptr(const std::string *a);
    GUESTOBJECT* from_ptr(const casadi::Slice *a);
    GUESTOBJECT* from_ptr(const casadi::Sparsity *a);
    GUESTOBJECT* from_ptr(const casadi::DMatrix *a);
    GUESTOBJECT* from_ptr(const casadi::IMatrix *a);
    GUESTOBJECT* from_ptr(const casadi::SX *a);
    GUESTOBJECT* from_ptr(const casadi::MX *a);
    GUESTOBJECT* from_ptr(const casadi::Function *a);
#ifdef SWIGMATLAB
    GUESTOBJECT* from_ptr(const std::pair<int, int>* a);
#endif // SWIGMATLAB
    template<typename M1, typename M2> GUESTOBJECT* from_ptr(const std::pair<M1, M2>* a);
#ifdef SWIGMATLAB
    GUESTOBJECT* from_ptr(const std::vector<std::string> *a);
#endif // SWIGMATLAB
    template<typename M> GUESTOBJECT* from_ptr(const std::vector<M> *a);
    template<typename M> GUESTOBJECT* from_ptr(const std::map<std::string, M> *a);

    // Same as the above, but with reference instead of pointer
    template<typename M> GUESTOBJECT* from_ref(const M& m) { return from_ptr(&m);}

    // Same as the above, but with a temporary object
    template<typename M> GUESTOBJECT* from_tmp(M m) { return from_ptr(&m);}
#ifdef SWIGMATLAB
    // Get sparsity pattern
    Sparsity getSparsity(const mxArray* p);

    // Number of nonzeros
    size_t getNNZ(const mxArray* p);
#endif // SWIGMATLAB

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

    // Check if int
    template<typename T> struct is_int {
      static inline bool check() {return false;}
    };

    template<> struct is_int<int> {
      static inline bool check() {return true;}
    };

    // Traits for assign vector
    template<typename E, typename M> struct traits_assign_vector {
      inline static bool assign(E* d, int sz, std::vector<M>** m) {
        // Not allowed by default
        return false;
      }
    };

    // int-to-int
    template<> struct traits_assign_vector<int, int> {
      inline static bool assign(int* d, int sz, std::vector<int>** m) {
        if (m) **m = std::vector<int>(d, d+sz);
        return true;
      }
    };

    // long-to-int
    template<> struct traits_assign_vector<long, int> {
      inline static bool assign(long* d, int sz, std::vector<int>** m) {
        if (m) **m = std::vector<int>(d, d+sz);
        return true;
      }
    };

    // long-to-double
    template<> struct traits_assign_vector<long, double> {
      inline static bool assign(long* d, int sz, std::vector<double>** m) {
        if (m) **m = std::vector<double>(d, d+sz);
        return true;
      }
    };

    // int-to-double
    template<> struct traits_assign_vector<int, double> {
      inline static bool assign(int* d, int sz, std::vector<double>** m) {
        if (m) **m = std::vector<double>(d, d+sz);
        return true;
      }
    };

    // double-to-double
    template<> struct traits_assign_vector<double, double> {
      inline static bool assign(double* d, int sz, std::vector<double>** m) {
        if (m) **m = std::vector<double>(d, d+sz);
        return true;
      }
    };

    // Assign to a vector, if conversion is allowed
    template<typename E, typename M> bool assign_vector(E* d, int sz, std::vector<M>** m) {
      return traits_assign_vector<E, M>::assign(d, sz, m);
    }

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
    Sparsity getSparsity(const mxArray* p) {
      // Get sparsity pattern
      size_t nrow = mxGetM(p);
      size_t ncol = mxGetN(p);

      if (mxIsSparse(p)) {
        // Sparse storage in MATLAB
        mwIndex *Jc = mxGetJc(p);
        mwIndex *Ir = mxGetIr(p);

        // Store in vectors
        std::vector<int> colind(ncol+1);
        std::copy(Jc, Jc+colind.size(), colind.begin());
        std::vector<int> row(colind.back());
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

%fragment("casadi_int", "header", fragment="casadi_aux", fragment=SWIG_AsVal_frag(int), fragment=SWIG_AsVal_frag(long)) {
  namespace casadi {
    bool to_ptr(GUESTOBJECT *p, int** m) {
      // Treat Null
      if (is_null(p)) return false;

      // Standard typemaps
      if (SWIG_IsOK(SWIG_AsVal(int)(p, m ? *m : 0))) return true;

#ifdef SWIGPYTHON
      // Numpy integer
      if (PyArray_IsScalar(p, Integer)) {
        int tmp = PyArray_PyIntAsInt(p);
        if (!PyErr_Occurred()) {
          if (m) **m = tmp;
          return true;
        }
        PyErr_Clear();
      }
#endif // SWIGPYTHON

      // long within int bounds
      {
        long tmp;
        if (SWIG_IsOK(SWIG_AsVal(long)(p, &tmp))) {
          // Check if within bounds
          if (tmp>=std::numeric_limits<int>::min() && tmp<=std::numeric_limits<int>::max()) {
            if (m) **m = static_cast<int>(tmp);
            return true;
          }
        }
      }

      // Scalar IMatrix
      {
        IMatrix *m2;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2), $descriptor(casadi::Matrix<int>*), 0))
            && m2->isscalar()) {
          if (m) **m = m2->getIntValue();
          return true;
        }
      }

      // No match
      return false;
    }

    GUESTOBJECT * from_ptr(const int *a) {
#ifdef SWIGPYTHON
      return PyInt_FromLong(*a);
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

      // Scalar DMatrix
      {
        DMatrix *m2;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2), $descriptor(casadi::Matrix<double>*), 0))
            && m2->isscalar()) {
          if (m) **m = m2->getValue();
          return true;
        }
      }

      // Scalar IMatrix
      {
        IMatrix *m2;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2), $descriptor(casadi::Matrix<int>*), 0))
            && m2->isscalar()) {
          if (m) **m = m2->getValue();
          return true;
        }
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
    // MATLAB n-by-m char array mapped to vector of length m

#ifdef SWIGMATLAB
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
      } else {
        return false;
      }
    }
#endif // SWIGMATLAB

    template<typename M> bool to_ptr(GUESTOBJECT *p, std::vector<M>** m) {
      // Treat Null
      if (is_null(p)) return false;
#ifdef SWIGPYTHON
      // 1D numpy array
      if (is_array(p) && array_numdims(p)==1 && array_type(p)!=NPY_OBJECT && array_is_native(p)) {
        int sz = array_size(p,0);

        // Make sure we have a contigous array with int datatype
        int array_is_new_object;
        PyArrayObject* array;

        // Trying NPY_INT
        if (assign_vector<int, M>(0, 0, 0)) {
          array = obj_to_array_contiguous_allow_conversion(p, NPY_INT, &array_is_new_object);
          if (array) {
            int *d = reinterpret_cast<int*>(array_data(array));
            int flag = assign_vector(d, sz, m);
            if (array_is_new_object) Py_DECREF(array);
            return flag;
          }
        }

        // Trying NPY_LONG
        if (assign_vector<long, M>(0, 0, 0)) {
          array = obj_to_array_contiguous_allow_conversion(p, NPY_LONG, &array_is_new_object);
          if (array) {
            long* d= reinterpret_cast<long*>(array_data(array));
            int flag = assign_vector(d, sz, m);
            if (array_is_new_object) Py_DECREF(array);
            return flag;
          }
        }

        // Trying NPY_DOUBLE
        if (assign_vector<double, M>(0, 0, 0)) {
          array = obj_to_array_contiguous_allow_conversion(p, NPY_DOUBLE, &array_is_new_object);
          if (array) {
            double* d= reinterpret_cast<double*>(array_data(array));
            int flag = assign_vector(d, sz, m);
            if (array_is_new_object) Py_DECREF(array);
            return flag;
          }
        }

        // No match
        return false;
      }
      // Python sequence
      if (PyList_Check(p) || PyTuple_Check(p)) {

        // Iterator to the sequence
        PyObject *it = PyObject_GetIter(p);
        if (!it) {
          PyErr_Clear();
          return false;
        }

        // Get size
        Py_ssize_t sz = PySequence_Size(p);
        if (sz==-1) {
          PyErr_Clear();
          return false;
        }

        // Allocate elements
        if (m) {
          (**m).clear();
          (**m).reserve(sz);
        }

        // Temporary
        M tmp;

        // Iterate over sequence
        for (Py_ssize_t i=0; i!=sz; ++i) {
          PyObject *pe=PyIter_Next(it);
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
      }
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
      // Cell arrays (only row vectors)
      if (mxGetClassID(p)==mxCELL_CLASS) {
        int nrow = mxGetM(p), ncol = mxGetN(p);
        if (nrow==1 || (nrow==0 && ncol==0)) {
          // Allocate elements
          if (m) {
            (**m).clear();
            (**m).reserve(ncol);
          }

          // Temporary
          M tmp;

          // Loop over elements
          for (int i=0; i<ncol; ++i) {
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
#endif // SWIGMATLAB
      // No match
      return false;
    }

#ifdef SWIGMATLAB
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
      for (int k=0; k<a->size(); ++k) {
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
      for (int k=0; k<a->size(); ++k) {
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

%fragment("casadi_function", "header", fragment="casadi_aux") {
  namespace casadi {
    bool to_ptr(GUESTOBJECT *p, Function** m) {
      // Treat Null
      if (is_null(p)) return false;

      // GenericType already?
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
      if (to_generic<int>(p, m)
          || to_generic<double>(p, m)
          || to_generic<std::string>(p, m)
          || to_generic<std::vector<int> >(p, m)
          || to_generic<std::vector<double> >(p, m)
          || to_generic<std::vector<std::string> >(p, m)
          || to_generic<std::vector<std::vector<int> > >(p, m)
          || to_generic<casadi::Function>(p, m)
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
      case OT_BOOLEAN: return from_tmp(a->asBool());
      case OT_INTEGER: return from_tmp(a->asInt());
      case OT_REAL: return from_tmp(a->asDouble());
      case OT_STRING: return from_tmp(a->asString());
      case OT_INTEGERVECTOR: return from_tmp(a->asIntVector());
      case OT_INTEGERVECTORVECTOR: return from_tmp(a->asIntVectorVector());
      case OT_REALVECTOR: return from_tmp(a->asDoubleVector());
      case OT_STRINGVECTOR: return from_tmp(a->asStringVector());
      case OT_DICT: return from_tmp(a->asDict());
      case OT_FUNCTION: return from_tmp(a->asFunction());
#ifdef SWIGPYTHON
      case OT_NULL: return Py_None;
#endif // SWIGPYTHON
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
      if (PyString_Check(p)) {
        if (m) (*m)->clear();
        if (m) (*m)->append(PyString_AsString(p));
        return true;
      }
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
      if (mxIsChar(p) && mxGetM(p)==1) {
	if (m) {
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
      // Python int
      if (PyInt_Check(p)) {
        if (m) {
          (**m).start_ = PyInt_AsLong(p);
          (**m).stop_ = (**m).start_+1;
          if ((**m).stop_==0) (**m).stop_ = std::numeric_limits<int>::max();
        }
        return true;
      }
      // Python slice
      if (PySlice_Check(p)) {
        PySliceObject *r = (PySliceObject*)(p);
        if (m) {
          (**m).start_ = (r->start == Py_None || PyInt_AsLong(r->start) < std::numeric_limits<int>::min())
            ? std::numeric_limits<int>::min() : PyInt_AsLong(r->start);
          (**m).stop_  = (r->stop ==Py_None || PyInt_AsLong(r->stop)> std::numeric_limits<int>::max())
            ? std::numeric_limits<int>::max() : PyInt_AsLong(r->stop);
          if(r->step !=Py_None) (**m).step_  = PyInt_AsLong(r->step);
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
          if (!PyString_Check(key)) return false;
          if (m) {
            M *v=&(**m)[std::string(PyString_AsString(key))], *v2=v;
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
	int len = mxGetNumberOfFields(p);
	for (int k=0; k<len; ++k) {
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
	  for (int k=0; k<fields.size(); ++k) mxDestroyArray(fields[k]);
	  return 0;
	}
	fields.push_back(f);	
      }
      
      // Create return object
      mxArray *p = mxCreateStructMatrix(1, 1, fields.size(),
					fieldnames.empty() ? 0 : &fieldnames[0]);
      for (int k=0; k<fields.size(); ++k) mxSetFieldByNumber(p, 0, k, fields[k]);
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
    bool to_ptr(GUESTOBJECT *p, std::pair<int, int>** m) {
      // (int,int) mapped to 2-by-1 double matrix
      if (mxIsDouble(p) && mxGetNumberOfDimensions(p)==2 && !mxIsSparse(p)
          && mxGetM(p)==1 && mxGetN(p)==2) {
        double* data = static_cast<double*>(mxGetData(p));
        int first = static_cast<int>(data[0]);
        int second = static_cast<int>(data[1]);
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
    GUESTOBJECT* from_ptr(const std::pair<int, int>* a) {
      // (int,int) mapped to 2-by-1 double matrix
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


%fragment("casadi_dvector", "header", fragment="casadi_aux") {
  namespace casadi {
    int to_ptr(GUESTOBJECT *p, std::vector<double> **m) {
      // Treat Null
      if (is_null(p)) return false;

      // Convert to DMatrix
      DMatrix tmp, *tmp_ptr=&tmp;
      if (to_ptr(p, &tmp_ptr) && tmp_ptr->iscolumn()) {
        if (m) tmp_ptr->get(**m);
        return true;
      }

      // No match
      return false;
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
                                    $descriptor(casadi::Matrix<casadi::SXElement>*), 0))) {
        return true;
      }

      // Object is an DMatrix
      {
        // Pointer to object
        DMatrix *m2;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2),
                                      $descriptor(casadi::Matrix<double>*), 0))) {
          if (m) **m=*m2;
          return true;
        }
      }

      // Object is an IMatrix
      {
        // Pointer to object
        IMatrix *m2;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2),
                                      $descriptor(casadi::Matrix<int>*), 0))) {
          if (m) **m=*m2;
          return true;
        }
      }

      // Object is a sparsity pattern
      {
        Sparsity *m2;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2),
                                      $descriptor(casadi::Sparsity*), 0))) {
          if (m) **m=SX::ones(*m2);
          return true;
        }
      }

      // Double scalar
      {
        double tmp;
        if (to_val(p, &tmp)) {
          if (m) **m=tmp;
          return true;
        }
      }

      // Integer scalar
      {
        int tmp;
        if (to_val(p, &tmp)) {
          if (m) **m=tmp;
          return true;
        }
      }

      // Try first converting to a temporary DMatrix
      {
        DMatrix tmp, *mt=&tmp;
        if(casadi::to_ptr(p, m ? &mt : 0)) {
          if (m) **m = *mt;
          return true;
        }
      }

#ifdef SWIGPYTHON
      // Numpy arrays will be cast to dense SX
      if (is_array(p)) {
        if (array_type(p) != NPY_OBJECT) return false;
        if (array_numdims(p)>2 || array_numdims(p)<1) return false;
        int nrows = array_size(p,0); // 1D array is cast into column vector
        int ncols  = array_numdims(p)==2 ? array_size(p,1) : 1;
        PyArrayIterObject* it = (PyArrayIterObject*)PyArray_IterNew(p);
        casadi::SX mT;
        if (m) mT = casadi::SX::zeros(ncols, nrows);
        int k=0;
        casadi::SX tmp, *tmp2;
        PyObject *pe;
        while (it->index < it->size) {
          pe = *((PyObject**) PyArray_ITER_DATA(it));
          tmp2=&tmp;
          if (!to_ptr(pe, &tmp2) || !tmp2->isscalar()) {
            Py_DECREF(it);
            return false;
          }
          if (m) mT(k++) = *tmp2;
          PyArray_ITER_NEXT(it);
        }
        Py_DECREF(it);
        if (m) **m = mT.T();
        return true;
      }
      // Object has __SX__ method
      if (PyObject_HasAttrString(p,"__SX__")) {
        char cmd[] = "__SX__";
        PyObject *cr = PyObject_CallMethod(p, cmd, 0);
        if (!cr) return false;
        int flag = to_ptr(cr, m);
        Py_DECREF(cr);
        return flag;
      }
#endif // SWIGPYTHON

      // No match
      return false;
    }

    GUESTOBJECT* from_ptr(const SX *a) {
      return SWIG_NewPointerObj(new SX(*a), $descriptor(casadi::Matrix<casadi::SXElement> *), SWIG_POINTER_OWN);
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

      // Object is an DMatrix
      {
        // Pointer to object
        DMatrix *m2;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2),
                                      $descriptor(casadi::Matrix<double>*), 0))) {
          if (m) **m=*m2;
          return true;
        }
      }

      // Object is a sparsity pattern
      {
        Sparsity *m2;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2),
                                      $descriptor(casadi::Sparsity*), 0))) {
          if (m) **m=MX::ones(*m2);
          return true;
        }
      }

      // Try first converting to a temporary DMatrix
      {
        DMatrix tmp, *mt=&tmp;
        if(casadi::to_ptr(p, m ? &mt : 0)) {
          if (m) **m = *mt;
          return true;
        }
      }

#ifdef SWIGPYTHON
      if (PyObject_HasAttrString(p,"__MX__")) {
        char cmd[] = "__MX__";
        PyObject *cr = PyObject_CallMethod(p, cmd, 0);
        if (!cr) return false;
        int flag = to_ptr(cr, m);
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

      bool ret = strcmp(PyString_AsString(classname),name)==0;
      Py_DECREF(classo);Py_DECREF(classname);
      return ret;
    }
#endif // SWIGPYTHON

    bool to_ptr(GUESTOBJECT *p, DMatrix** m) {
      // Treat Null
      if (is_null(p)) return false;

      // DMatrix already?
      if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                    $descriptor(casadi::Matrix<double>*), 0))) {
        return true;
      }

      // Object is an IMatrix
      {
        // Pointer to object
        IMatrix *m2;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2),
                                      $descriptor(casadi::Matrix<int>*), 0))) {
          if (m) **m=*m2;
          return true;
        }
      }

      // Object is a sparsity pattern
      {
        Sparsity *m2;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2),
                                      $descriptor(casadi::Sparsity*), 0))) {
          if (m) **m=DMatrix::ones(*m2);
          return true;
        }
      }

      // Double scalar
      {
        double tmp;
        if (to_val(p, &tmp)) {
          if (m) **m=tmp;
          return true;
        }
      }

      // Integer scalar
      {
        int tmp;
        if (to_val(p, &tmp)) {
          if (m) **m=tmp;
          return true;
        }
      }

#ifdef SWIGPYTHON
      // Object has __DMatrix__ method
      if (PyObject_HasAttrString(p,"__DMatrix__")) {
        char name[] = "__DMatrix__";
        PyObject *cr = PyObject_CallMethod(p, name, 0);
        if (!cr) return false;
        int result = to_val(cr, m ? *m : 0);
        Py_DECREF(cr);
        return result;
      }
      // Numpy arrays will be cast to dense Matrix<double>
      if (is_array(p)) {
        int array_is_new_object;
        PyArrayObject* array = obj_to_array_contiguous_allow_conversion(p, NPY_DOUBLE, &array_is_new_object);
        if (!array) return false;
        int nrow, ncol;
        switch (array_numdims(p)) {
        case 0:
          // Scalar
          nrow=ncol=1;
          break;
        case 1:
          // Vector
          nrow=array_size(p, 0);
          ncol=1;
          break;
        case 2:
          // Matrix
          nrow=array_size(p, 0);
          ncol=array_size(p, 1);
          break;
        default:
          // More than two dimension unsupported
          if (array_is_new_object) Py_DECREF(array);
          return false;
        }
        if (m) {
          **m = casadi::Matrix<double>::zeros(nrow, ncol);
          casadi::Matrix<double>::iterator it=(*m)->begin();
          double* d = reinterpret_cast<double*>(array_data(array));
          for (int cc=0; cc<ncol; ++cc) {
            for (int rr=0; rr<nrow; ++rr) {
              *it++ = d[cc+rr*ncol];
            }
          }
        }

        // Free memory
        if (array_is_new_object) Py_DECREF(array);
        return true;
      }

      // scipy's csc_matrix will be cast to sparse DMatrix
      if(PyObjectHasClassName(p, "csc_matrix")) {

        // Get the dimensions of the csc_matrix
        PyObject * shape = PyObject_GetAttrString( p, "shape"); // need's to be decref'ed
        if (!shape) return false;
        if(!PyTuple_Check(shape) || PyTuple_Size(shape)!=2) {
          Py_DECREF(shape);
          return false;
        }
        int nrows=PyInt_AsLong(PyTuple_GetItem(shape,0));
        int ncols=PyInt_AsLong(PyTuple_GetItem(shape,1));
        Py_DECREF(shape);

        bool ret= false;

        PyObject * narray=0;
        PyObject * row=0;
        PyObject * colind=0;
        PyArrayObject* array=0;
        PyArrayObject* array_row=0;
        PyArrayObject* array_colind=0;

        int array_is_new_object=0;
        int row_is_new_object=0;
        int colind_is_new_object=0;

        // Fetch data
        narray=PyObject_GetAttrString( p, "data"); // need's to be decref'ed
        if (!narray || !is_array(narray) || array_numdims(narray)!=1) goto cleanup;
        array = obj_to_array_contiguous_allow_conversion(narray,NPY_DOUBLE,&array_is_new_object);
        if (!array) goto cleanup;

        // Construct the 'row' vector needed for initialising the correct sparsity
        row = PyObject_GetAttrString(p,"indices"); // need's to be decref'ed
        if (!row || !is_array(row) || array_numdims(row)!=1) goto cleanup;
        array_row = obj_to_array_contiguous_allow_conversion(row,NPY_INT,&row_is_new_object);
        if (!array_row) goto cleanup;

        // Construct the 'colind' vector needed for initialising the correct sparsity
        colind = PyObject_GetAttrString(p,"indptr"); // need's to be decref'ed
        if (!colind || !is_array(colind) || array_numdims(colind)!=1) goto cleanup;
        array_colind = obj_to_array_contiguous_allow_conversion(colind,NPY_INT,&colind_is_new_object);
        if (!array_colind) goto cleanup;
        {
          int size=array_size(array,0); // number on non-zeros
          double* d=(double*) array_data(array);
          std::vector<double> v(d,d+size);

          int* rowd=(int*) array_data(array_row);
          std::vector<int> rowv(rowd,rowd+size);

          int* colindd=(int*) array_data(array_colind);
          std::vector<int> colindv(colindd,colindd+(ncols+1));

          if (m) **m = casadi::Matrix<double>(casadi::Sparsity(nrows,ncols,colindv,rowv), v, false);

          ret = true;
        }

      cleanup: // yes that's right; goto.
        // Rather that than a pyramid of conditional memory-deallocation
        // TODO(jaeandersson): Create a helper struct and put the below in the destructor
        if (array_is_new_object && array) Py_DECREF(array);
        if (narray) Py_DECREF(narray);
        if (row_is_new_object && array_row) Py_DECREF(array_row);
        if (row) Py_DECREF(row);
        if (colind_is_new_object && array_colind) Py_DECREF(array_colind);
        if (colind) Py_DECREF(colind);
        return ret;
      }
      if(PyObject_HasAttrString(p,"tocsc")) {
        char name[] = "tocsc";
        PyObject *cr = PyObject_CallMethod(p, name,0);
        if (!cr) return false;
        int result = to_val(cr, m ? *m : 0);
        Py_DECREF(cr);
        return result;
      }

      {
        std::vector <double> t;
        int res = to_val(p, &t);
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
          **m = casadi::DMatrix(getSparsity(p));
          double* data = static_cast<double*>(mxGetData(p));
          (*m)->setNZ(data);
        }
        return true;
      }
#endif // SWIGMATLAB

      // First convert to IMatrix
      if (can_convert<IMatrix>(p)) {
        IMatrix tmp;
        if (to_val(p, &tmp)) {
          if (m) **m=tmp;
          return true;
        }
      }

      // No match
      return false;
    }

    GUESTOBJECT* from_ptr(const DMatrix *a) {
      return SWIG_NewPointerObj(new DMatrix(*a), $descriptor(casadi::Matrix<double>*), SWIG_POINTER_OWN);
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
    bool to_ptr(GUESTOBJECT *p, IMatrix** m) {
      // Treat Null
      if (is_null(p)) return false;

      // IMatrix already?
      if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                    $descriptor(casadi::Matrix<int>*), 0))) {
        return true;
      }

      // Object is a sparsity pattern
      {
        Sparsity *m2;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2),
                                      $descriptor(casadi::Sparsity*), 0))) {
          if (m) **m=IMatrix::ones(*m2);
          return true;
        }
      }

      // First convert to integer
      {
        int tmp;
        if (to_val(p, &tmp)) {
          if (m) **m=tmp;
          return true;
        }
      }

#ifdef SWIGPYTHON
      // Numpy arrays will be cast to dense Matrix<int>
      if (is_array(p)) {
        int array_is_new_object;
        bool is_long=false;
        PyArrayObject* array = obj_to_array_contiguous_allow_conversion(p, NPY_INT, &array_is_new_object);
        if (!array) {
          // Trying NPY_LONG
          is_long=true;
          PyErr_Clear();
          array = obj_to_array_contiguous_allow_conversion(p, NPY_LONG, &array_is_new_object);
        }
        if (!array) return false;
        int nrow, ncol;
        switch (array_numdims(p)) {
        case 0:
          // Scalar
          nrow=ncol=1;
          break;
        case 1:
          // Vector
          nrow=array_size(p, 0);
          ncol=1;
          break;
        case 2:
          // Matrix
          nrow=array_size(p, 0);
          ncol=array_size(p, 1);
          break;
        default:
          // More than two dimension unsupported
          if (array_is_new_object) Py_DECREF(array);
          return false;
        }
        if (m) {
          **m = casadi::Matrix<int>::zeros(nrow, ncol);
          casadi::Matrix<int>::iterator it=(*m)->begin();
          if (is_long) {
            long* d = reinterpret_cast<long*>(array_data(array));
            for (int cc=0; cc<ncol; ++cc) {
              for (int rr=0; rr<nrow; ++rr) {
                *it++ = d[cc+rr*ncol];
              }
            }
          } else {
            int* d = reinterpret_cast<int*>(array_data(array));
            for (int cc=0; cc<ncol; ++cc) {
              for (int rr=0; rr<nrow; ++rr) {
                *it++ = d[cc+rr*ncol];
              }
            }
          }
        }

        // Free memory
        if (array_is_new_object) Py_DECREF(array);
        return true;
      }

      if (PyObject_HasAttrString(p,"__IMatrix__")) {
        char cmd[] = "__IMatrix__";
        PyObject *cr = PyObject_CallMethod(p, cmd, 0);
        if (!cr) return false;
        int result = to_val(cr, m ? *m : 0);
        Py_DECREF(cr);
        return result;
      }

      {
        std::vector <int> t;
        int res = to_val(p, &t);
        if (m) **m = casadi::Matrix<int>(t);
        return res;
      }
      return true;
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
      // In MATLAB, it is common to use floating point values to represent integers
      if (mxIsDouble(p) && mxGetNumberOfDimensions(p)==2) {
        double* data = static_cast<double*>(mxGetData(p));

        // Check if all integers
        bool all_integers=true;
        size_t sz = getNNZ(p);
        for (size_t i=0; i<sz; ++i) {
          if (data[i] != int(data[i])) {
            all_integers = false;
            break;
          }
        }

        // If successful
        if (all_integers) {
          if (m) {
            **m = casadi::IMatrix(getSparsity(p));
            for (size_t i=0; i<sz; ++i) {
              (*m)->at(i) = int(data[i]);
            }
          }
          return true;
        }
      }
#endif // SWIGMATLAB

      // No match
      return false;
    }

    GUESTOBJECT* from_ptr(const IMatrix *a) {
      return SWIG_NewPointerObj(new IMatrix(*a), $descriptor(casadi::Matrix<int>*), SWIG_POINTER_OWN);
    }
  } // namespace casadi
 }

// Collect all fragments
%fragment("casadi_all", "header", fragment="casadi_aux,casadi_bool,casadi_int,casadi_double,casadi_vector,casadi_function,casadi_generictype,casadi_string,casadi_slice,casadi_map,casadi_pair,casadi_dvector,casadi_sx,casadi_mx,casadi_dmatrix,casadi_sparsity,casadi_imatrix") { }

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
%typemap(in, noblock=1, fragment="casadi_all") xType {
  if (!casadi::to_val($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Cannot convert input to " xName ".");
 }

 // Pass input by value, cleanup
%typemap(freearg, noblock=1) xType {}

 // Pass input by reference, check if matches
%typemap(typecheck, noblock=1, precedence=xPrec, fragment="casadi_all") const xType& {
  $1 = casadi::to_ptr($input, static_cast< xType **>(0));
 }

 // Pass input by reference, convert argument
%typemap(in, noblock=1, fragment="casadi_all") const xType & (xType m) {
  $1 = &m;
  if (!casadi::to_ptr($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input to " xName ".");
 }

 // Pass input by reference, cleanup
%typemap(freearg, noblock=1) const xType & {}
%enddef

 // Define all output typemaps
%define %casadi_output_typemaps(xName, xType...)

 // Return-by-value
%typemap(out, noblock=1, fragment="casadi_all") xType, const xType {
  if(!($result = casadi::from_ref(static_cast<const xType &>($1)))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to " xName ".");
}

// Return a const-ref behaves like return-by-value
%typemap(out, noblock=1, fragment="casadi_all") const xType& {
  if(!($result = casadi::from_ptr(static_cast<const xType *>($1)))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to " xName ".");
}

// Inputs marked OUTPUT are also returned by the function, ...
%typemap(argout,noblock=1,fragment="casadi_all") xType &OUTPUT {
  %append_output(casadi::from_ptr(static_cast<xType *>($1)));
 }

// ... and the corresponding inputs are ignored
%typemap(in, numinputs=0) xType &OUTPUT (xType m) {
 $1 = &m;
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
  %append_output(casadi::from_ptr(static_cast<xType *>($1)));
 }

// ... but kept as inputs
%typemap(in, noblock=1, fragment="casadi_all") xType &INOUT (xType m) {
  $1 = &m;
  if (!casadi::to_ptr($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input to " xName ".");
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
%define PREC_GENERICTYPE 22 %enddef
%define PREC_DICT 21 %enddef
%define PREC_SPARSITY 90 %enddef
%define PREC_IVector 92 %enddef
%define PREC_IVectorVector 92 %enddef
%define PREC_VECTOR 92 %enddef
%define PREC_PAIR_SLICE_SLICE 93 %enddef
%define PREC_SLICE 94 %enddef
%define PREC_PAIR_IVector_IVector 96 %enddef
%define PREC_IMatrix 97 %enddef
%define PREC_IMatrixVector 98 %enddef
%define PREC_IMatrixVectorVector 98 %enddef
%define PREC_DVector 99 %enddef
%define PREC_DMatrix 100 %enddef
%define PREC_DMatrixVector 101 %enddef
%define PREC_DMatrixVectorVector 101 %enddef
%define PREC_SX 103 %enddef
%define PREC_SXVector 103 %enddef
%define PREC_SXVectorVector 103 %enddef
%define PREC_MX 104 %enddef
%define PREC_MXVector 105 %enddef
%define PREC_MXVectorVector 106 %enddef
%define PREC_CREATOR 150 %enddef
%define PREC_STRING 180 %enddef
%define PREC_FUNCTION 200 %enddef

#ifndef SWIGXML

 // std::ostream & is not typemapped to anything useful and should be ignored
 // (or possibly turned into a string output) 
%typemap(in, numinputs=0) std::ostream &stream ""

%casadi_typemaps("str", PREC_STRING, std::string)
%casadi_template("[str]", PREC_STRING, std::vector<std::string>)
%casadi_typemaps("Sparsity", PREC_SPARSITY, casadi::Sparsity)
%casadi_template("[Sparsity]", PREC_SPARSITY, std::vector< casadi::Sparsity>)
%casadi_template("[[Sparsity]]", PREC_SPARSITY, std::vector<std::vector< casadi::Sparsity> >)
%casadi_template("str:Sparsity", PREC_SPARSITY, std::map<std::string, casadi::Sparsity >)
%casadi_template("str:[Sparsity]", PREC_SPARSITY, std::map<std::string, std::vector<casadi::Sparsity > >)
%casadi_template("(str:Sparsity,[str])", PREC_SPARSITY, std::pair<std::map<std::string, casadi::Sparsity >, std::vector<std::string> >)
%casadi_typemaps("bool", SWIG_TYPECHECK_BOOL, bool)
%casadi_template("[bool]", SWIG_TYPECHECK_BOOL, std::vector<bool>)
%casadi_template("[[bool]]", SWIG_TYPECHECK_BOOL, std::vector<std::vector<bool> >)
%casadi_typemaps("int", SWIG_TYPECHECK_INTEGER, int)
%casadi_template("(int,int)", SWIG_TYPECHECK_INTEGER, std::pair<int,int>)
%casadi_template("[int]", PREC_IVector, std::vector<int>)
%casadi_template("[[int]]", PREC_IVectorVector, std::vector<std::vector<int> >)
%casadi_typemaps("double", SWIG_TYPECHECK_DOUBLE, double)
%casadi_template("[double]", SWIG_TYPECHECK_DOUBLE, std::vector<double>)
%casadi_template("[[double]]", SWIG_TYPECHECK_DOUBLE, std::vector<std::vector<double> >)
%casadi_typemaps("SX", PREC_SX, casadi::Matrix<casadi::SXElement>)
%casadi_template("[SX]", PREC_SXVector, std::vector< casadi::Matrix<casadi::SXElement> >)
%casadi_template("[[SX]]", PREC_SXVectorVector, std::vector<std::vector< casadi::Matrix<casadi::SXElement> > >)
%casadi_template("str:SX", PREC_SX, std::map<std::string, casadi::Matrix<casadi::SXElement> >)
%casadi_template("(str:SX,[str])", PREC_SX, std::pair<std::map<std::string, casadi::Matrix<casadi::SXElement> >, std::vector<std::string> >)
%casadi_typemaps("MX", PREC_MX, casadi::MX)
%casadi_template("[MX]", PREC_MXVector, std::vector<casadi::MX>)
%casadi_template("[[MX]]", PREC_MXVectorVector, std::vector<std::vector<casadi::MX> >)
%casadi_template("str:MX", PREC_MX, std::map<std::string, casadi::MX>)
%casadi_template("(str:MX,[str])", PREC_MX, std::pair<std::map<std::string, casadi::MX >, std::vector<std::string> >)
%casadi_typemaps("DMatrix", PREC_DMatrix, casadi::Matrix<double>)
%casadi_template("[DMatrix]", PREC_DMatrixVector, std::vector< casadi::Matrix<double> >)
%casadi_template("[[DMatrix]]", PREC_DMatrixVectorVector, std::vector<std::vector< casadi::Matrix<double> > >)
%casadi_template("str:DMatrix", PREC_DMatrix, std::map<std::string, casadi::Matrix<double> >)
%casadi_template("(str:DMatrix,[str])", PREC_DMatrix, std::pair<std::map<std::string, casadi::Matrix<double> >, std::vector<std::string> >)
%casadi_typemaps("IMatrix", PREC_IMatrix, casadi::Matrix<int>)
%casadi_template("[IMatrix]", PREC_IMatrixVector, std::vector< casadi::Matrix<int> >)
%casadi_template("[[IMatrix]]", PREC_IMatrixVectorVector, std::vector<std::vector< casadi::Matrix<int> > >)
%casadi_typemaps("GenericType", PREC_GENERICTYPE, casadi::GenericType)
%casadi_template("[GenericType]", PREC_GENERICTYPE, std::vector<casadi::GenericType>)
%casadi_typemaps("Slice", PREC_SLICE, casadi::Slice)
%casadi_typemaps("Function", PREC_FUNCTION, casadi::Function)
%casadi_template("[Function]", PREC_FUNCTION, std::vector<casadi::Function>)
%casadi_template("(Function,Function)", PREC_FUNCTION, std::pair<casadi::Function, casadi::Function>)
%casadi_template("Dict", PREC_DICT, std::map<std::string, casadi::GenericType>)

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
import _casadi
%}
#endif // SWIGPYTHON

%include <casadi/core/function/schemes_metadata.hpp>

// Init hooks
#ifdef SWIGPYTHON
#ifdef WITH_PYTHON_INTERRUPTS
%{
#include <pythonrun.h>
void SigIntHandler(int) {
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

try:
  from numpy import sin, cos, tan, sqrt, log, exp, floor, ceil, fmod, fmin, fmax, sinh, cosh, tanh, arcsin, arccos, arctan, arctan2, fabs, sign, arctanh, arcsinh, arccosh, copysign
except:
  sin = lambda x: x.sin()
  cos = lambda x: x.cos()
  tan = lambda x: x.tan()
  arcsin = lambda x: x.arcsin()
  arccos = lambda x: x.arccos()
  arctan = lambda x: x.arctan()
  sqrt = lambda x: x.sqrt()
  log = lambda x: x.log()
  exp = lambda x: x.exp()
  floor = lambda x: x.floor()
  ceil = lambda x: x.ceil()
  fmin = lambda x,y: x.fmin(y)
  fmax = lambda x,y: x.fmax(y)
  sinh = lambda x: x.sinh()
  cosh = lambda x: x.cosh()
  tanh = lambda x: x.tanh()
  fabs = lambda x: x.fabs()
  sign = lambda x: x.sign()
  arctan2 = lambda x,y: x.arctan2(y)
  arctanh = lambda x: x.arctanh()
  arcsinh = lambda x: x.arcsinh()
  arccosh = lambda x: x.arccosh()
  copysign = lambda x,y: x.copysign(y)
%}
#endif // SWIGPYTHON

%rename("%(regex:/friendwrap_(?!ML)(.*)/\\1/)s") ""; // Strip leading friendwrap_ unless followed by ML
#ifndef SWIGMATLAB
%rename("%(regex:/zz_(?!ML)(.*)/\\1/)s") ""; // Strip leading zz_ unless followed by ML
#endif

%rename(row) get_row;
%rename(colind) get_colind;
%rename(sparsity) getSparsity;

#ifdef SWIGPYTHON
%ignore T;
%rename(__float__) getValue;
%rename(__int__) getIntValue;

%rename(logic_and) friendwrap_and;
%rename(logic_or) friendwrap_or;
%rename(logic_not) friendwrap_not;
%rename(logic_all) friendwrap_all;
%rename(logic_any) friendwrap_any;

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
%rename(mtimes) friendwrap_mul;
%feature("varargin","1") friendwrap_vertcat;
%feature("varargin","1") friendwrap_horzcat;
%feature("varargin","1") friendwrap_veccat;
%feature("optionalunpack","1") size;

// Raise an error if "this" not correct
%typemap(check) SWIGTYPE *self %{
if (!$1) {
  SWIG_Error(SWIG_RuntimeError, "Invalid 'self' object");
  SWIG_fail;
 }
%}

// Workarounds, pending proper fix
%rename(nonzero) __nonzero__;
%rename(hash) __hash__;
#endif // SWIGMATLAB

#ifdef SWIGPYTHON
%pythoncode %{
def prod(self,*args):
    raise Exception("'prod' is not supported anymore in CasADi. Use 'mul' to do matrix multiplication.")
def dot(self,*args):
    raise Exception("'dot' is not supported anymore in CasADi. Use 'mul' to do matrix multiplication.")

class NZproxy:
  def __init__(self,matrix):
    self.matrix = matrix

  def __getitem__(self,s):
    return self.matrix.getNZ(False, s)

  def __setitem__(self,s,val):
    return self.matrix.setNZ(val, False, s)

  def __len__(self):
    return self.matrix.nnz()

  def __iter__(self):
    for k in range(len(self)):
      yield self[k]
%}

%define %matrix_convertors
%pythoncode %{

    def toMatrix(self):
        import numpy as n
        return n.matrix(self.toArray())

    def __iter__(self):
      for k in self.nz:
        yield k

%}
%enddef
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
        with internalAPI():
          if isinstance(s, tuple) and len(s)==2:
            return self.get(False, s[0], s[1])
          return self.get(False, s)

    def __setitem__(self,s,val):
        with internalAPI():
          if isinstance(s,tuple) and len(s)==2:
            return self.set(val, False, s[0], s[1])
          return self.set(val, False, s)

    @property
    def nz(self):
      return NZproxy(self)

    def prod(self,*args):
        raise Exception("'prod' is not supported anymore in CasADi. Use 'mul' to do matrix multiplication.")

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
        return self.toArray()

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
      casadi_assert(ch==':');
      return Slice();
    }
  } // namespace casadi
%}

%define %matrix_helpers(Type)
    // Get a submatrix (index-1)
    const Type paren(char rr) const {
      casadi_assert(rr==':');
      return vec(*$self);
    }
    const Type paren(const Matrix<int>& rr) const {
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
    const Type paren(char rr, const Matrix<int>& cc) const {
      Type m;
      $self->get(m, true, casadi::char2Slice(rr), cc);
      return m;
    }
    const Type paren(const Matrix<int>& rr, char cc) const {
      Type m;
      $self->get(m, true, rr, casadi::char2Slice(cc));
      return m;
    }
    const Type paren(const Matrix<int>& rr, const Matrix<int>& cc) const {
      Type m;
      $self->get(m, true, rr, cc);
      return m;
    }

    // Set a submatrix (index-1)
    void paren_asgn(const Type& m, char rr) { $self->set(m, true, casadi::char2Slice(rr));}
    void paren_asgn(const Type& m, const Matrix<int>& rr) { $self->set(m, true, rr);}
    void paren_asgn(const Type& m, const Sparsity& sp) { $self->set(m, true, sp);}
    void paren_asgn(const Type& m, char rr, char cc) { $self->set(m, true, casadi::char2Slice(rr), casadi::char2Slice(cc));}
    void paren_asgn(const Type& m, char rr, const Matrix<int>& cc) { $self->set(m, true, casadi::char2Slice(rr), cc);}
    void paren_asgn(const Type& m, const Matrix<int>& rr, char cc) { $self->set(m, true, rr, casadi::char2Slice(cc));}
    void paren_asgn(const Type& m, const Matrix<int>& rr, const Matrix<int>& cc) { $self->set(m, true, rr, cc);}

    // Get nonzeros (index-1)
    const Type brace(char rr) const { Type m; $self->getNZ(m, true, casadi::char2Slice(rr)); return m;}
    const Type brace(const Matrix<int>& rr) const { Type m; $self->getNZ(m, true, rr); return m;}

    // Set nonzeros (index-1)
    void setbrace(const Type& m, char rr) { $self->setNZ(m, true, casadi::char2Slice(rr));}
    void setbrace(const Type& m, const Matrix<int>& rr) { $self->setNZ(m, true, rr);}

    // 'end' function (needed for end syntax in MATLAB)
    inline int end(int i, int n) const {
      return n==1 ? $self->numel() : i==1 ? $self->size1() : $self->size2();
    }

    // Transpose using the A' syntax in addition to A.'
    Type ctranspose() const { return $self->T();}

%enddef
#endif

#ifndef SWIGPYTHON
%define %matrix_convertors
%enddef
#endif

%include <casadi/core/printable_object.hpp>

#ifdef SWIGPYTHON
%rename(SWIG_STR) getDescription;
#endif // SWIGPYTHON

%template(PrintSharedObject) casadi::PrintableObject<casadi::SharedObject>;
%template(PrintSlice)        casadi::PrintableObject<casadi::Slice>;
%template(PrintIMatrix)      casadi::PrintableObject<casadi::Matrix<int> >;
%template(PrintDMatrix)      casadi::PrintableObject<casadi::Matrix<double> >;
//%template(PrintSX)           casadi::PrintableObject<casadi::Matrix<casadi::SXElement> >;
%template(PrintNlpBuilder)     casadi::PrintableObject<casadi::NlpBuilder>;
%template(PrintVariable)        casadi::PrintableObject<casadi::Variable>;
%template(PrintDaeBuilder)     casadi::PrintableObject<casadi::DaeBuilder>;

%include <casadi/core/shared_object.hpp>
%include <casadi/core/std_vector_tools.hpp>
%include <casadi/core/weak_ref.hpp>
%include <casadi/core/casadi_types.hpp>
%include <casadi/core/generic_type.hpp>
%include <casadi/core/options_functionality.hpp>

namespace casadi {
  %extend OptionsFunctionality {
    void setOption(const std::string &name, const std::string& val){$self->setOption(name,val);}
    void setOption(const std::string &name, const std::vector<int>& val){$self->setOption(name,val);}
    void setOption(const std::string &name, const std::vector<double>& val){$self->setOption(name,val);}
    void setOption(const std::string &name, double val){$self->setOption(name,val);}
    void setOption(const std::string &name, int val){$self->setOption(name,val);}
    void setOption(const std::string &name, bool val){$self->setOption(name,val);}
    void setOption(const std::string &name, const std::vector< std::vector<int> >& val){$self->setOption(name,val);}
  }
} // namespace casadi

%include <casadi/core/matrix/sparsity_interface.hpp>

%template(SpSparsity) casadi::SparsityInterface<casadi::Sparsity>;
%include <casadi/core/matrix/sparsity.hpp>

// Logic for pickling
#ifdef SWIGPYTHON
namespace casadi{
%extend Sparsity {
  %pythoncode %{
    def __setstate__(self, state):
        if state:
          self.__init__(state["nrow"],state["ncol"],state["colind"],state["row"])
        else:
          self.__init__()

    def __getstate__(self):
        if self.isNull(): return {}
        return {"nrow": self.size1(), "ncol": self.size2(), "colind": numpy.array(self.colind(),dtype=int), "row": numpy.array(self.row(),dtype=int)}
  %}
}

} // namespace casadi
#endif // SWIGPYTHON

/* There is no reason to expose the Slice class to e.g. Python or MATLAB. Only if an interfaced language
   lacks a slice type, the type should be exposed here */
// #if !(defined(SWIGPYTHON) || defined(SWIGMATLAB))
%include <casadi/core/matrix/slice.hpp>
 //#endif

%template(SpIMatrix)        casadi::SparsityInterface<casadi::Matrix<int> >;
%template(SpDMatrix)        casadi::SparsityInterface<casadi::Matrix<double> >;
%template(SpSX)             casadi::SparsityInterface<casadi::Matrix<casadi::SXElement> >;
%template(SpMX)             casadi::SparsityInterface<casadi::MX>;

%include <casadi/core/matrix/generic_matrix.hpp>

%template(GenIMatrix)        casadi::GenericMatrix<casadi::Matrix<int> >;
%template(GenDMatrix)        casadi::GenericMatrix<casadi::Matrix<double> >;
%template(GenSX)             casadi::GenericMatrix<casadi::Matrix<casadi::SXElement> >;
%template(GenMX)             casadi::GenericMatrix<casadi::MX>;

%include <casadi/core/matrix/generic_expression.hpp>

%template(ExpIMatrix)        casadi::GenericExpression<casadi::Matrix<int> >;
%template(ExpDMatrix)        casadi::GenericExpression<casadi::Matrix<double> >;
%template(ExpSX)             casadi::GenericExpression<casadi::Matrix<casadi::SXElement> >;
%template(ExpMX)             casadi::GenericExpression<casadi::MX>;

// Prefix symbols
#if defined(SWIGMATLAB) || defined(SWIGXML)
%define %HIDE(SYM) friendwrap_ ## SYM %enddef
#else
%define %HIDE(SYM) casadi_ ## SYM %enddef
#endif
%define %SHOW(SYM) friendwrap_ ## SYM %enddef

// Flags to allow differentiating the wrapping by type
#define IS_GLOBAL   0x1
#define IS_MEMBER   0x10
#define IS_SPARSITY 0x100
#define IS_DMATRIX  0x1000
#define IS_IMATRIX  0x10000
#define IS_SX       0x100000
#define IS_MX       0x1000000

%define SPARSITY_INTERFACE_FUN_BASE(DECL, FLAG, M)
#if FLAG & IS_MEMBER

 DECL M %SHOW(horzcat)(const std::vector< M > &v) {
  return horzcat(v);
 }
 DECL M %SHOW(vertcat)(const std::vector< M > &v) {
 return vertcat(v);
 }
 DECL std::vector< M >
 %SHOW(horzsplit)(const M& v, const std::vector<int>& offset) {
 return horzsplit(v, offset);
 }
 DECL std::vector< M > %SHOW(horzsplit)(const M& v, int incr=1) {
 return horzsplit(v, incr);
 }
 DECL std::vector< M >
 %SHOW(vertsplit)(const M& v, const std::vector<int>& offset) {
 return vertsplit(v, offset);
 }
 DECL std::vector<int >
 %SHOW(offset)(const std::vector< M > &v, bool vert=true) {
 return offset(v, vert);
 }
 DECL std::vector< M >
 %SHOW(vertsplit)(const M& v, int incr=1) {
 return vertsplit(v, incr);
 }
 DECL M %SHOW(blockcat)(const std::vector< std::vector< M > > &v) {
 return blockcat(v);
 }
 DECL M %SHOW(blockcat)(const M& A, const M& B, const M& C, const M& D) {
 return vertcat(horzcat(A, B), horzcat(C, D));
 }
 DECL std::vector< std::vector< M > >
 %SHOW(blocksplit)(const M& x, const std::vector<int>& vert_offset,
 const std::vector<int>& horz_offset) {
 return blocksplit(x, vert_offset, horz_offset);
 }
 DECL std::vector< std::vector< M > >
 %SHOW(blocksplit)(const M& x, int vert_incr=1, int horz_incr=1) {
 return blocksplit(x, vert_incr, horz_incr);
 }
 DECL M %SHOW(diagcat)(const std::vector< M > &A) {
 return diagcat(A);
 }
 DECL std::vector< M >
 %SHOW(diagsplit)(const M& x, const std::vector<int>& output_offset1,
 const std::vector<int>& output_offset2) {
 return diagsplit(x, output_offset1, output_offset2);
 }
 DECL std::vector< M >
 %SHOW(diagsplit)(const M& x, const std::vector<int>& output_offset) {
 return diagsplit(x, output_offset);
 }
 DECL std::vector< M > %SHOW(diagsplit)(const M& x, int incr=1) {
 return diagsplit(x, incr);
 }
 DECL std::vector< M >
 %SHOW(diagsplit)(const M& x, int incr1, int incr2) {
 return diagsplit(x, incr1, incr2);
 }
 DECL M %SHOW(veccat)(const std::vector< M >& x) {
 return veccat(x);
 }
 DECL M %SHOW(mul)(const M& X, const M& Y) {
 return mul(X, Y);
 }
 DECL M %SHOW(mul)(const std::vector< M > &args) {
 return mul(args);
 }
 DECL M %SHOW(mac)(const M& X, const M& Y, const M& Z) {
 return mac(X, Y, Z);
 }
 DECL M %SHOW(transpose)(const M& X) {
 return X.T();
 }
 DECL M %SHOW(vec)(const M& a) {
 return vec(a);
 }
 DECL M %SHOW(vecNZ)(const M& a) {
 return vecNZ(a);
 }
 DECL M %SHOW(reshape)(const M& a, int nrow, int ncol) {
 return reshape(a, nrow, ncol);
 }
 DECL M %SHOW(reshape)(const M& a, std::pair<int, int> rc) {
 return reshape(a, rc.first, rc.second);
 }
 DECL M %SHOW(reshape)(const M& a, const Sparsity& sp) {
 return reshape(a, sp);
 }
 DECL int %SHOW(sprank)(const M& A) {
 return sprank(A);
 }
 DECL int %SHOW(norm_0_mul)(const M& x, const M& y) {
 return norm_0_mul(x, y);
 }
 DECL M %SHOW(triu)(const M& a, bool includeDiagonal=true) {
 return triu(a, includeDiagonal);
 }
 DECL M %SHOW(tril)(const M& a, bool includeDiagonal=true) {
 return tril(a, includeDiagonal);
 }
 DECL M %SHOW(kron)(const M& a, const M& b) {
 return kron(a, b);
 }
 DECL M %SHOW(repmat)(const M& A, int n, int m=1) {
 return repmat(A, n, m);
 }
 DECL M %SHOW(repmat)(const M& A, const std::pair<int, int>& rc) {
 return repmat(A, rc.first, rc.second);
 }
#endif
%enddef

%define SPARSITY_INTERFACE_ALL(DECL, FLAG)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_SPARSITY), Sparsity)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_MX), MX)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_IMATRIX), Matrix<int>)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_DMATRIX), Matrix<double>)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_SX), Matrix<SXElement>)
%enddef

#ifdef SWIGMATLAB
  %define SPARSITY_INTERFACE_FUN(DECL, FLAG, M)
    SPARSITY_INTERFACE_FUN_BASE(DECL, FLAG, M)
    #if FLAG & IS_MEMBER
     DECL int %SHOW(length)(const M &v) {
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
DECL M %SHOW(mpower)(const M& x, const M& n) {
  return mpower(x, n);
}

DECL M %HIDE(mrdivide)(const M& x, const M& y) {
  return mrdivide(x, y);
}

DECL M %HIDE(mldivide)(const M& x, const M& y) {
  return mldivide(x, y);
}

DECL std::vector< M > %SHOW(symvar)(const M& x) {
  return symvar(x);
}

DECL M %SHOW(quad_form)(const M& X, const M& A) {
  return quad_form(X, A);
}

DECL M %SHOW(quad_form)(const M& X) {
  return quad_form(X);          
}

DECL M %SHOW(sum_square)(const M& X) {
  return sum_square(X);         
}

DECL M %SHOW(linspace)(const M& a, const M& b, int nsteps) {
  return linspace(a, b, nsteps);
}

DECL M %SHOW(cross)(const M& a, const M& b, int dim = -1) {
  return cross(a, b, dim);      
}

DECL M %SHOW(det)(const M& A) {
  return det(A);                
}

DECL M %SHOW(inv)(const M& A) {
  return inv(A);                
}

DECL M %SHOW(trace)(const M& a) {
  return trace(a);              
}

DECL M %SHOW(tril2symm)(const M& a) {
  return tril2symm(a);          
}

DECL M %SHOW(triu2symm)(const M& a) {
  return triu2symm(a);          
}

DECL M %SHOW(norm_F)(const M& x) {
  return norm_F(x);             
}

DECL M %SHOW(norm_2)(const M& x) {
  return norm_2(x);             
}

DECL M %SHOW(norm_1)(const M& x) {
  return norm_1(x);             
}

DECL M %SHOW(norm_inf)(const M& x) {
  return norm_inf(x);           
}

DECL M %SHOW(sumCols)(const M& x) {
  return sumCols(x);            
}

DECL M %SHOW(sumRows)(const M& x) {
  return sumRows(x);            
}

DECL M %SHOW(inner_prod)(const M& x, const M& y) {
  return inner_prod(x, y);      
}

DECL M %SHOW(outer_prod)(const M& x, const M& y) {
  return outer_prod(x, y);      
}

DECL M %SHOW(nullspace)(const M& A) {
  return nullspace(A);          
}

DECL M %SHOW(polyval)(const M& p, const M& x) {
  return polyval(p, x);         
}

DECL M %SHOW(diag)(const M& A) {
  return diag(A);               
}

DECL M %SHOW(unite)(const M& A, const M& B) {
  return unite(A, B);           
}

DECL M %SHOW(densify)(const M& x) {
  return densify(x);            
}

DECL M %SHOW(project)(const M& A, const Sparsity& sp, bool intersect=false) {
  return project(A, sp, intersect);    
}

DECL M %SHOW(if_else)(const M& cond, const M& if_true, 
                    const M& if_false, bool short_circuit=true) {
  return if_else(cond, if_true, if_false, short_circuit);   
}

DECL M %SHOW(conditional)(const M& ind, const std::vector< M > &x,
                        const M& x_default, bool short_circuit=true) {
  return conditional(ind, x, x_default, short_circuit);     
}

DECL bool %SHOW(dependsOn)(const M& f, const M& arg) {
  return dependsOn(f, arg);     
}

DECL M %SHOW(solve)(const M& A, const M& b) {
  return solve(A, b);
}

DECL M %SHOW(solve)(const M& A, const M& b,
                       const std::string& lsolver,
                       const casadi::Dict& opts = casadi::Dict()) {
  return solve(A, b, lsolver, opts);
}

DECL M %SHOW(pinv)(const M& A) {
  return pinv(A);
}

DECL M %SHOW(pinv)(const M& A, const std::string& lsolver,
                      const casadi::Dict& opts = casadi::Dict()) {
  return pinv(A, lsolver, opts);
}

DECL M %SHOW(jacobian)(const M &ex, const M &arg) {
  return jacobian(ex, arg);
}

DECL M %SHOW(gradient)(const M &ex, const M &arg) {
  return gradient(ex, arg);
}

DECL M %SHOW(tangent)(const M &ex, const M &arg) {
  return tangent(ex, arg);
}

DECL M %SHOW(hessian)(const M& ex, const M& arg, M& OUTPUT1) {
  return hessian(ex, arg, OUTPUT1);
}

DECL int %SHOW(countNodes)(const M& A) {
  return countNodes(A);
}

DECL std::string %SHOW(getOperatorRepresentation)(const M& xb,
                                                  const std::vector<std::string>& args) {
  return getOperatorRepresentation(xb, args);
}
DECL M %SHOW(repsum)(const M& A, int n, int m=1) {
  return repsum(A, n, m);
}

#endif // FLAG & IS_MEMBER

#if FLAG & IS_GLOBAL
DECL M %SHOW(substitute)(const M& ex, const M& v, const M& vdef) {
  return substitute(ex, v, vdef);
}

DECL std::vector< M > %SHOW(substitute)(const std::vector< M >& ex,
                                         const std::vector< M >& v,
                                         const std::vector< M >& vdef) {
  return substitute(ex, v, vdef);
}

DECL void %SHOW(substituteInPlace)(const std::vector< M >& v,
                                      std::vector< M >& INOUT1,
                                      std::vector< M >& INOUT2,
                                      bool reverse=false) {
  return substituteInPlace(v, INOUT1, INOUT2, reverse);
}

DECL void %SHOW(extractShared)(const std::vector< M >& ex,
                               std::vector< M >& OUTPUT1,
                               std::vector< M >& OUTPUT2,
                               std::vector< M >& OUTPUT3,
                               const std::string& v_prefix="v_",
                               const std::string& v_suffix="") {
  extractShared(ex, OUTPUT1, OUTPUT2, OUTPUT3, v_prefix, v_suffix);
}

#endif // FLAG & IS_GLOBAL
%enddef

%define GENERIC_MATRIX_ALL(DECL, FLAG)
GENERIC_MATRIX_FUN(DECL, (FLAG | IS_MX), MX)
GENERIC_MATRIX_FUN(DECL, (FLAG | IS_IMATRIX), Matrix<int>)
GENERIC_MATRIX_FUN(DECL, (FLAG | IS_DMATRIX), Matrix<double>)
GENERIC_MATRIX_FUN(DECL, (FLAG | IS_SX), Matrix<SXElement>)
%enddef

%define GENERIC_EXPRESSION_FUN(DECL, FLAG, M) 
#if FLAG & IS_MEMBER
DECL M %HIDE(plus)(const M& x, const M& y) { return x+y; }
DECL M %HIDE(minus)(const M& x, const M& y) { return x-y; }
DECL M %HIDE(times)(const M& x, const M& y) { return x*y; }
DECL M %HIDE(rdivide)(const M& x, const M& y) { return x/y; }
DECL M %HIDE(ldivide)(const M& x, const M& y) { return y/x; }
DECL M %HIDE(lt)(const M& x, const M& y) { return x<y; }
DECL M %HIDE(le)(const M& x, const M& y) { return x<=y; }
DECL M %HIDE(gt)(const M& x, const M& y) { return x>y; }
DECL M %HIDE(ge)(const M& x, const M& y) { return x>=y; }
DECL M %HIDE(eq)(const M& x, const M& y) { return x==y; }
DECL M %HIDE(ne)(const M& x, const M& y) { return x!=y; }
DECL M %SHOW(and)(const M& x, const M& y) { return x&&y; }
DECL M %SHOW(or)(const M& x, const M& y) { return x||y; }
DECL M %SHOW(not)(const M& x) { return !x; }
DECL M %HIDE(abs)(const M& x) { return fabs(x); }
DECL M %HIDE(sqrt)(const M& x) { return sqrt(x); }
DECL M %HIDE(sin)(const M& x) { return sin(x); }
DECL M %HIDE(cos)(const M& x) { return cos(x); }
DECL M %HIDE(tan)(const M& x) { return tan(x); }
DECL M %HIDE(atan)(const M& x) { return atan(x); }
DECL M %HIDE(asin)(const M& x) { return asin(x); }
DECL M %HIDE(acos)(const M& x) { return acos(x); }
DECL M %HIDE(tanh)(const M& x) { return tanh(x); }
DECL M %HIDE(sinh)(const M& x) { return sinh(x); }
DECL M %HIDE(cosh)(const M& x) { return cosh(x); }
DECL M %HIDE(atanh)(const M& x) { return atanh(x); }
DECL M %HIDE(asinh)(const M& x) { return asinh(x); }
DECL M %HIDE(acosh)(const M& x) { return acosh(x); }
DECL M %HIDE(exp)(const M& x) { return exp(x); }
DECL M %HIDE(log)(const M& x) { return log(x); }
DECL M %HIDE(log10)(const M& x) { return log10(x); }
DECL M %HIDE(floor)(const M& x) { return floor(x); }
DECL M %HIDE(ceil)(const M& x) { return ceil(x); }
DECL M %HIDE(erf)(const M& x) { return erf(x); }
DECL M %SHOW(erfinv)(const M& x) { return erfinv(x); }
DECL M %HIDE(sign)(const M& x) { return sign(x); }
DECL M %HIDE(power)(const M& x, const M& n) { return pow(x, n); }
DECL M %HIDE(mod)(const M& x, const M& y) { return fmod(x, y); }
DECL M %HIDE(atan2)(const M& x, const M& y) { return atan2(x, y); }
DECL M %HIDE(min)(const M& x, const M& y) { return fmin(x, y); }
DECL M %HIDE(max)(const M& x, const M& y) { return fmax(x, y); }
DECL M %SHOW(simplify)(const M& x) { return simplify(x); }
DECL bool %SHOW(is_equal)(const M& x, const M& y, int depth=0) { return is_equal(x, y, depth); }
DECL bool %SHOW(iszero)(const M& x) { return iszero(x); }
DECL M %HIDE(copysign)(const M& x, const M& y) { return copysign(x, y); }
DECL M %HIDE(constpow)(const M& x, const M& y) { return constpow(x, y); }
#endif // FLAG & IS_MEMBER
%enddef

%define GENERIC_EXPRESSION_ALL(DECL, FLAG) 
GENERIC_EXPRESSION_FUN(DECL, (FLAG | IS_MX), MX)
GENERIC_EXPRESSION_FUN(DECL, (FLAG | IS_IMATRIX), Matrix<int>)
GENERIC_EXPRESSION_FUN(DECL, (FLAG | IS_DMATRIX), Matrix<double>)
GENERIC_EXPRESSION_FUN(DECL, (FLAG | IS_SX), Matrix<SXElement>)
%enddef

%define MATRIX_FUN(DECL, FLAG, M)
#if FLAG & IS_MEMBER
DECL M %SHOW(all)(const M& x) {
  return all(x);
}

DECL M %SHOW(any)(const M& x) {
  return any(x);
}

DECL M %SHOW(adj)(const M& A) {
  return adj(A);
}

DECL M %SHOW(getMinor)(const M& x, int i, int j) {
  return getMinor(x, i, j);
}

DECL M %SHOW(cofactor)(const M& x, int i, int j) {
  return cofactor(x, i, j);
}

DECL void %SHOW(qr)(const M& A, M& OUTPUT1, M& OUTPUT2) {
  return qr(A, OUTPUT1, OUTPUT2);
}

DECL M %SHOW(chol)(const M& A) {
  return chol(A);
}

DECL M %SHOW(norm_inf_mul)(const M& x, const M& y) {
  return norm_inf_mul(x, y);
}

DECL M %SHOW(sparsify)(const M& A, double tol=0) {
  return sparsify(A, tol);
}

DECL void %SHOW(expand)(const M& ex, M& OUTPUT1, M& OUTPUT2) {
  expand(ex, OUTPUT1, OUTPUT2);
}

DECL M %SHOW(pw_const)(const M &t, const M& tval, const M& val) {
  return pw_const(t, tval, val);
}

DECL M %SHOW(pw_lin)(const M& t, const M& tval, const M& val) {
  return pw_lin(t, tval, val);
}

DECL M %SHOW(heaviside)(const M& x) {
  return heaviside(x);
}

DECL M %SHOW(rectangle)(const M& x) {
  return rectangle(x);
}

DECL M %SHOW(triangle)(const M& x) {
  return triangle(x);
}

DECL M %SHOW(ramp)(const M& x) {
  return ramp(x);
}

DECL M %SHOW(gauss_quadrature)(const M& f, const M& x,
                               const M& a, const M& b,
                               int order=5) {
  return gauss_quadrature(f, x, a, b, order);
}

DECL M %SHOW(gauss_quadrature)(const M& f, const M& x,
                               const M& a, const M& b,
                               int order, const M& w) {
  return gauss_quadrature(f, x, a, b, order, w);
}

DECL M %SHOW(jacobianTimesVector)(const M& ex, const M& arg, const M& v,
                                  bool transpose_jacobian=false) {
  return jacobianTimesVector(ex, arg, v, transpose_jacobian);
}

DECL M %SHOW(taylor)(const M& ex, const M& x, const M& a=0, int order=1) {
  return taylor(ex, x, a, order);
}

DECL M %SHOW(mtaylor)(const M& ex, const M& x, const M& a, int order=1) {
  return mtaylor(ex, x, a, order);
}

DECL M %SHOW(mtaylor)(const M& ex, const M& x, const M& a, int order,
                      const std::vector<int>& order_contributions) {
  return mtaylor(ex, x, a, order, order_contributions);
}

DECL M %SHOW(poly_coeff)(const M& ex,
                         const M&x) {
  return poly_coeff(ex, x);
}

DECL M %SHOW(poly_roots)(const M& p) {
  return poly_roots(p);
}

DECL M %SHOW(eig_symbolic)(const M& m) {
  return eig_symbolic(m);
}

#endif
%enddef

%define MATRIX_ALL(DECL, FLAG)
MATRIX_FUN(DECL, (FLAG | IS_IMATRIX), Matrix<int>)
MATRIX_FUN(DECL, (FLAG | IS_DMATRIX), Matrix<double>)
MATRIX_FUN(DECL, (FLAG | IS_SX), Matrix<SXElement>)
%enddef

%define MX_FUN(DECL, FLAG, M)
#if FLAG & IS_MEMBER
DECL M %SHOW(find)(const M& x) {
  return find(x);
}
#endif // FLAG & IS_MEMBER

#if FLAG & IS_GLOBAL
DECL std::vector< M >
%SHOW(matrix_expand)(const std::vector< M >& e,
                     const std::vector< M > &boundary = std::vector< M >(),
                     const Dict& options = Dict()) {
  return matrix_expand(e, boundary, options);
}

DECL M %SHOW(matrix_expand)(const M& e,
                            const std::vector< M > &boundary = std::vector< M >(),
                            const Dict& options = Dict()) {
  return matrix_expand(e, boundary, options);
}

DECL M %SHOW(graph_substitute)(const M& ex, const std::vector< M >& v,
                         const std::vector< M > &vdef) {
  return graph_substitute(ex, v, vdef);
}

DECL std::vector< M >
%SHOW(graph_substitute)(const std::vector< M > &ex,
                 const std::vector< M > &v,
                 const std::vector< M > &vdef) {
  return graph_substitute(ex, v, vdef);
}

#endif
%enddef

%define MX_ALL(DECL, FLAG)
MX_FUN(DECL, (FLAG | IS_MX), MX)
%enddef

%template(PrintSX)           casadi::PrintableObject<casadi::Matrix<casadi::SXElement> >;

%include <casadi/core/matrix/matrix.hpp>

%template(IMatrix)           casadi::Matrix<int>;
%template(DMatrix)           casadi::Matrix<double>;

%extend casadi::Matrix<double> {
   %template(DMatrix) Matrix<int>;
};

namespace casadi{
  %extend Matrix<double> {

    void assign(const casadi::Matrix<double>&rhs) { (*$self)=rhs; }
    %matrix_convertors
    %matrix_helpers(casadi::Matrix<double>)

  }
  %extend Matrix<int> {

    void assign(const casadi::Matrix<int>&rhs) { (*$self)=rhs; }
    %matrix_convertors
    %matrix_helpers(casadi::Matrix<int>)

  }
}

// Extend DMatrix with SWIG unique features
namespace casadi{
  %extend Matrix<double> {
    // Convert to a dense matrix
    GUESTOBJECT* full() const {
#ifdef SWIGPYTHON
      npy_intp dims[2] = {$self->size1(), $self->size2()};
      PyObject* ret = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
      double* d = static_cast<double*>(array_data(ret));
      $self->get(d, true); // Row-major
      return ret;
#elif defined(SWIGMATLAB)
      mxArray *p  = mxCreateDoubleMatrix($self->size1(), $self->size2(), mxREAL);
      double* d = static_cast<double*>(mxGetData(p));
      $self->get(d); // Column-major
      return p;
#else
      return 0;
#endif
    }

#ifdef SWIGMATLAB
    // Convert to a sparse matrix
    GUESTOBJECT* sparse() const {
      mxArray *p  = mxCreateSparse($self->size1(), $self->size2(), $self->nnz(), mxREAL);
      $self->getNZ(static_cast<double*>(mxGetData(p)));
      std::copy($self->colind(), $self->colind()+$self->size2()+1, mxGetJc(p));
      std::copy($self->row(), $self->row()+$self->nnz(), mxGetIr(p));
      return p;
    }
#endif
  }
} // namespace casadi


#ifdef SWIGPYTHON
namespace casadi{
%extend Matrix<double> {
%pythoncode %{
  def toArray(self):
    import numpy as n
    if isinstance(self,IMatrix):
      return n.array(self.get(),n.int).reshape((self.shape[1],self.shape[0])).T
    else:
      return n.array(self.get()).reshape((self.shape[1],self.shape[0])).T
%}

%python_array_wrappers(999.0)

// The following code has some trickery to fool numpy ufunc.
// Normally, because of the presence of __array__, an ufunctor like nump.sqrt
// will unleash its activity on the output of __array__
// However, we wish DMatrix to remain a DMatrix
// So when we receive a call from a functor, we return a dummy empty array
// and return the real result during the postprocessing (__array_wrap__) of the functor.
%pythoncode %{
  def __array_custom__(self,*args,**kwargs):
    if "dtype" in kwargs and not(isinstance(kwargs["dtype"],n.double)):
      return n.array(self.toArray(),dtype=kwargs["dtype"])
    else:
      return self.toArray()
%}

%pythoncode %{
  def toCsc_matrix(self):
    import numpy as n
    import warnings
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      from scipy.sparse import csc_matrix
    return csc_matrix( (self.nonzeros(),self.row(),self.colind()), shape = self.shape, dtype=n.double )

  def tocsc(self):
    return self.toCsc_matrix()

%}

%pythoncode %{
  def __nonzero__(self):
    if self.numel()!=1:
      raise Exception("Only a scalar can be cast to a float")
    if self.nnz()==0:
      return 0
    return float(self)!=0
%}

%pythoncode %{
  def __abs__(self):
    return abs(float(self))
%}

}; // extend Matrix<double>

%extend Matrix<int> {

  %python_array_wrappers(998.0)

  %pythoncode %{
    def toArray(self):
      import numpy as n
      r = n.zeros((self.size1(),self.size2()))
      d = self.nonzeros_int()
      for j in range(self.size2()):
        for k in range(self.colind(j),self.colind(j+1)):
          i = self.row(k)
          r[i,j] = d[k]
      return r
  %}

  %pythoncode %{
    def __abs__(self):
      return abs(int(self))
  %}
} // extend Matrix<int>


// Logic for pickling

%extend Matrix<int> {

  %pythoncode %{
    def __setstate__(self, state):
        sp = Sparsity.__new__(Sparsity)
        sp.__setstate__(state["sparsity"])
        self.__init__(sp,state["data"])

    def __getstate__(self):
        return {"sparsity" : self.sparsity().__getstate__(), "data": numpy.array(self.nonzeros_int(),dtype=int)}
  %}
}

%extend Matrix<double> {

  %pythoncode %{
    def __setstate__(self, state):
        sp = Sparsity.__new__(Sparsity)
        sp.__setstate__(state["sparsity"])
        self.__init__(sp,state["data"])

    def __getstate__(self):
        return {"sparsity" : self.sparsity().__getstate__(), "data": numpy.array(self.nonzeros(),dtype=float)}
  %}

}


} // namespace casadi
#endif // SWIGPYTHON

%include <casadi/core/sx/sx_element.hpp>

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
            return DMatrix.ones(self).toArray()
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

  fmin_backup = fmin
  fmax_backup = fmax

  def fmin(x,y):
    pass

  def fmax(x,y):
    pass

  _min_ufunc = numpy.frompyfunc(fmin,2,1)
  _max_ufunc = numpy.frompyfunc(fmax,2,1)

  fmin = fmin_backup
  fmax = fmax_backup

  _defaultmin = min
  def min(*args,**kwargs):
    if len(args)==2 and len(kwargs)==0 and (hasattr(args[0],'fmin') or hasattr(args[1],'fmin')):
      return _min_ufunc(*args)
    else:
      return _defaultmin(*args,**kwargs)

  _defaultmax = max
  def max(*args,**kwargs):
    if len(args)==2 and len(kwargs)==0 and (hasattr(args[0],'fmax') or hasattr(args[1],'fmax')):
      return _max_ufunc(*args)
    else:
      return _defaultmax(*args,**kwargs)
except:
  pass
%}
#endif // SWIGPYTHON

namespace casadi {
%extend Matrix<SXElement>{

    %matrix_convertors
    %matrix_helpers(casadi::Matrix<casadi::SXElement>)

    #ifdef SWIGPYTHON
    %pythoncode %{
    def toArray(self):
      import numpy as n
      r = n.array((),dtype=object)
      r.resize(self.size1(),self.size2())
      for j in range(self.size2()):
        for el in range(self.colind(j),self.colind(j+1)):
          i=self.row(el)
          r[i,j] = self.nz[el]
      return r
    %}

  %python_array_wrappers(1001.0)
  #endif // SWIGPYTHON

};

} // namespace casadi

#ifdef SWIGPYTHON
#include <arrayobject.h>
%template()    std::vector<PyObject*>;
#endif // SWIGPYTHON

%template(SX) casadi::Matrix<casadi::SXElement>;
%extend casadi::Matrix<casadi::SXElement> {
   %template(SX) Matrix<int>;
   %template(SX) Matrix<double>;
};

%include <casadi/core/mx/mx.hpp>

%extend casadi::MX{
  %matrix_helpers(casadi::MX)

  #ifdef SWIGPYTHON
  %python_array_wrappers(1002.0)

  %pythoncode %{
  def __array_custom__(self,*args,**kwargs):
    import numpy as np
    if np.__version__=="1.8.1": #1083
      return np.array(np.nan)
    raise Exception("MX cannot be converted to an array. MX.__array__ purely exists to allow ufunc/numpy goodies")

  def __iter__(self):
    return self.nz.__iter__()

  %}
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

%include <casadi/core/function/io_interface.hpp>

%template(IOInterfaceFunction) casadi::IOInterface<casadi::Function>;

%extend casadi::IOInterface<casadi::Function> {
  casadi::Matrix<double> getInput(int iind=0) const             { static_cast<const casadi::Function*>(return $self->input(iind);}
  casadi::Matrix<double> getInput(const std::string &iname) const             { return $self->input($self->inputIndex(iname)); }
  casadi::Matrix<double> getOutput(int oind=0) const            { static_cast<const casadi::Function*>(return $self->output(oind);}
}

%include <casadi/core/function/io_scheme.hpp>
%include <casadi/core/function/function.hpp>
%feature("copyctor", "0") casadi::CodeGenerator;
%include <casadi/core/function/code_generator.hpp>

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
namespace casadi {
  %extend GenericExpressionCommon {
    %pythoncode %{
      def __add__(x, y): return casadi_plus(x, y)
      def __radd__(x, y): return casadi_plus(y, x)
      def __sub__(x, y): return casadi_minus(x, y)
      def __rsub__(x, y): return casadi_minus(y, x)
      def __mul__(x, y): return casadi_times(x, y)
      def __rmul__(x, y): return casadi_times(y, x)
      def __div__(x, y): return casadi_rdivide(x, y)
      def __rdiv__(x, y): return casadi_rdivide(y, x)
      def __truediv__(x, y): return casadi_rdivide(x, y)
      def __rtruediv__(x, y): return casadi_rdivide(y, x)
      def __lt__(x, y): return casadi_lt(x, y)
      def __rlt__(x, y): return casadi_lt(y, x)
      def __le__(x, y): return casadi_le(x, y)
      def __rle__(x, y): return casadi_le(y, x)
      def __gt__(x, y): return casadi_lt(y, x)
      def __rgt__(x, y): return casadi_lt(x, y)
      def __ge__(x, y): return casadi_le(y, x)
      def __rge__(x, y): return casadi_le(x, y)
      def __eq__(x, y): return casadi_eq(x, y)
      def __req__(x, y): return casadi_eq(y, x)
      def __ne__(x, y): return casadi_ne(x, y)
      def __rne__(x, y): return casadi_ne(y, x)
      def __pow__(x, n): return casadi_power(x, n)
      def __rpow__(n, x): return casadi_power(x, n)
      def __arctan2__(x, y): return casadi_atan2(x, y)
      def __rarctan2__(y, x): return casadi_atan2(x, y)
      def fmin(x, y): return casadi_min(x, y)
      def fmax(x, y): return casadi_max(x, y)
      def __fmin__(x, y): return casadi_min(x, y)
      def __rfmin__(y, x): return casadi_min(x, y)
      def __fmax__(x, y): return casadi_max(x, y)
      def __rfmax__(y, x): return casadi_max(x, y)
      def logic_and(x, y): return casadi_and(x, y)
      def logic_or(x, y): return casadi_or(x, y)
      def fabs(x): return casadi_abs(x)
      def sqrt(x): return casadi_sqrt(x)
      def sin(x): return casadi_sin(x)
      def cos(x): return casadi_cos(x)
      def tan(x): return casadi_tan(x)
      def arcsin(x): return casadi_asin(x)
      def arccos(x): return casadi_acos(x)
      def arctan(x): return casadi_atan(x)
      def sinh(x): return casadi_sinh(x)
      def cosh(x): return casadi_cosh(x)
      def tanh(x): return casadi_tanh(x)
      def arcsinh(x): return casadi_asinh(x)
      def arccosh(x): return casadi_acosh(x)
      def arctanh(x): return casadi_atanh(x)
      def exp(x): return casadi_exp(x)
      def log(x): return casadi_log(x)
      def log10(x): return casadi_log10(x)
      def floor(x): return casadi_floor(x)
      def ceil(x): return casadi_ceil(x)
      def erf(x): return casadi_erf(x)
      def sign(x): return casadi_sign(x)
      def fmod(x, y): return casadi_mod(x, y)
      def __copysign__(x, y): return casadi_copysign(x, y)
      def __rcopysign__(y, x): return casadi_copysign(x, y)
      def copysign(x, y): return casadi_copysign(x, y)
      def rcopysign(y, x): return casadi_copysign(x, y)
      def __constpow__(x, y): return casadi_constpow(x, y)
      def __rconstpow__(y, x): return casadi_constpow(x, y)
      def constpow(x, y): return casadi_constpow(x, y)
      def rconstpow(y, x): return casadi_constpow(x, y)
    %}
  }

  %extend GenericMatrixCommon {
    %pythoncode %{
      def __mldivide__(x, y): return casadi_mldivide(x, y)
      def __rmldivide__(y, x): return casadi_mldivide(x, y)
      def __mrdivide__(x, y): return casadi_mrdivide(x, y)
      def __rmrdivide__(y, x): return casadi_mrdivide(x, y)
      def __mpower__(x, y): return casadi_mpower(x, y)
      def __rmpower__(y, x): return casadi_mpower(x, y)
    %}
  }

} // namespace casadi
#endif // SWIGPYTHON

%feature("director") casadi::Callback;

%include <casadi/core/function/compiler.hpp>
%include <casadi/core/function/linear_solver.hpp>
%include <casadi/core/function/implicit_function.hpp>
%include <casadi/core/function/integrator.hpp>
%include <casadi/core/function/simulator.hpp>
%include <casadi/core/function/nlp_solver.hpp>
%include <casadi/core/function/qp_solver.hpp>
%include <casadi/core/function/stabilized_qp_solver.hpp>
%include <casadi/core/function/callback.hpp>
%include <casadi/core/function/nullspace.hpp>

%include "autogenerated.i"

%include <casadi/core/casadi_options.hpp>
%include <casadi/core/casadi_meta.hpp>
%include <casadi/core/misc/integration_tools.hpp>
%include <casadi/core/misc/nlp_builder.hpp>
%include <casadi/core/misc/variable.hpp>
%include <casadi/core/misc/dae_builder.hpp>
%include <casadi/core/misc/xml_file.hpp>
#ifdef SWIGPYTHON
%pythoncode %{


def swig_typename_convertor_cpp2python(s):
  import re
  s = s.replace("C/C++ prototypes","Python usages")
  s = s.replace("casadi::","")
  s = s.replace("MXDict","str:MX")
  s = s.replace("SXDict","str:SX")
  s = s.replace("std::string","str")
  s = s.replace(" const &","")
  s = s.replace("friendwrap_","")
  s = re.sub(r"\b((\w+)(< \w+ >)?)::\2\b",r"\1",s)
  s = re.sub("(const )?Matrix< ?SXElement *>( &)?",r"SX",s)
  s = re.sub("(const )?GenericMatrix< ?(\w+) *>( ?&)?",r"\2 ",s)
  s = re.sub("(const )?Matrix< ?int *>( ?&)?",r"IMatrix ",s)
  s = re.sub("(const )?Matrix< ?double *>( ?&)?",r"DMatrix ",s)
  s = re.sub("(const )?Matrix< ?(\w+) *>( ?&)?",r"array(\2) ",s)
  s = re.sub("(const )?GenericMatrix< ?([\w\(\)]+) *>( ?&)?",r"\2 ",s)
  s = re.sub(r"const (\w+) &",r"\1 ",s)
  s = re.sub(r"< [\w\(\)]+ +>\(",r"(",s)
  for i in range(5):
    s = re.sub(r"(const )? ?std::pair< ?([\w\(\)\]\[: ]+?) ?, ?([\w\(\)\]\[: ]+?) ?> ?&?",r"(\2,\3) ",s)
    s = re.sub(r"(const )? ?std::vector< ?([\w\(\)\[\] ]+) ?(, ?std::allocator< ?\2 ?>)? ?> ?&?",r"[\2] ",s)
  s = re.sub(r"\b(\w+)(< \w+ >)?::\1",r"\1",s)
  s = s.replace("casadi::","")
  s = s.replace("IOInterface< Function >","Function")
  s = s.replace("::",".")
  s = s.replace(".operator ()","")
  s = re.sub(r"([A-Z]\w+)Vector",r"[\1]",s)
  return s
  
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
    return "|".join(set([swig_typename_convertor_python2cpp(i) for i in a.keys()])) +":"+ "|".join(set([swig_typename_convertor_python2cpp(i) for i in a.values()]))

  return type(a).__name__


def swig_monkeypatch(v,cl=True):
  import re
  if hasattr(v,"__monkeypatched__"):
    return v
  def foo(*args,**kwargs):
    try:
      return v(*args,**kwargs)
    except NotImplementedError as e:
      import sys
      exc_info = sys.exc_info()
      if e.message.startswith("Wrong number or type of arguments for overloaded function"):

        s = e.args[0]
        s = s.replace("'new_","'")
        #s = re.sub(r"overloaded function '(\w+?)_(\w+)'",r"overloaded function '\1.\2'",s)
        m = re.search("overloaded function '([\w\.]+)'",s)
        if m:
          name = m.group(1)
          name = name.replace(".__call__","")
        else:
          name = "method"
        ne = NotImplementedError(swig_typename_convertor_cpp2python(s)+"You have: %s(%s)\n" % (name,", ".join(map(swig_typename_convertor_python2cpp,args[1:] if cl else args)+ ["%s=%s" % (k,swig_typename_convertor_python2cpp(vv)) for k,vv in kwargs.items()])))
        raise ne.__class__, ne, exc_info[2].tb_next
      else:
        raise exc_info[1], None, exc_info[2].tb_next
    except TypeError as e:
      import sys
      exc_info = sys.exc_info()
      
      methodname = "method"
      try:
        methodname = exc_info[2].tb_next.tb_frame.f_code.co_name
      except:
        pass

      if e.message.startswith("in method '"):
        s = e.args[0]
        s = re.sub(r"method '(\w+?)_(\w+)'",r"method '\1.\2'",s)
        m = re.search("method '([\w\.]+)'",s)
        if m:
          name = m.group(1)
          name = name.replace(".__call__","")
        else:
          name = "method"
        ne = TypeError(swig_typename_convertor_cpp2python(s)+" expected.\nYou have: %s(%s)\n" % (name,", ".join(map(swig_typename_convertor_python2cpp,args[1:] if cl else args))))
        raise ne.__class__, ne, exc_info[2].tb_next
      elif e.message.startswith("Expecting one of"):
        s = e.args[0]
        conversion = {"mul": "*", "div": "/", "add": "+", "sub": "-","le":"<=","ge":">=","lt":"<","gt":">","eq":"==","pow":"**"}
        if methodname.startswith("__") and methodname[2:-2] in conversion:
          ne = TypeError(swig_typename_convertor_cpp2python(s)+"\nYou try to do: %s %s %s.\n" % (  swig_typename_convertor_python2cpp(args[0]),conversion[methodname[2:-2]] ,swig_typename_convertor_python2cpp(args[1]) ))
        elif methodname.startswith("__r") and methodname[3:-2] in conversion:
          ne = TypeError(swig_typename_convertor_cpp2python(s)+"\nYou try to do: %s %s %s.\n" % ( swig_typename_convertor_python2cpp(args[1]),  conversion[methodname[3:-2]], swig_typename_convertor_python2cpp(args[0]) ))
        else:
          ne = TypeError(swig_typename_convertor_cpp2python(s)+"\nYou have: (%s)\n" % (", ".join(map(swig_typename_convertor_python2cpp,args[1:] if cl else args))))
        raise ne.__class__, ne, exc_info[2].tb_next
      else:
        s = e.args[0]
        ne = TypeError(s+"\nYou have: (%s)\n" % (", ".join(map(swig_typename_convertor_python2cpp,args[1:] if cl else args) + ["%s=%s" % (k,swig_typename_convertor_python2cpp(vv)) for k,vv in kwargs.items()]  )))
        raise ne.__class__, ne, exc_info[2].tb_next
    except AttributeError as e:
      import sys
      exc_info = sys.exc_info()
      if e.message=="type object 'object' has no attribute '__getattr__'":
        # swig 3.0 bug
        ne = AttributeError("Unkown attribute: %s has no attribute '%s'." % (str(args[1]),args[2]))
        raise ne.__class__, ne, exc_info[2].tb_next
      else:
        raise exc_info[1], None, exc_info[2].tb_next
    except Exception as e:
      import sys
      exc_info = sys.exc_info()
      raise exc_info[1], None, exc_info[2].tb_next
      
  if v.__doc__ is not None:
    foo.__doc__ = swig_typename_convertor_cpp2python(v.__doc__)
  foo.__name__ = v.__name__
  foo.__monkeypatched__ = True
  return foo

import inspect

def swig_improvedcall(v):
  def newcall(self,*args,**kwargs):
    if len(args)>0 and len(kwargs)>0:
      raise Exception("You cannot mix positional and keyword arguments in __call__")
    if len(kwargs)>0:
      return v(self,kwargs)
    else:
      return v(self,*args)

  newcall.__name__ = v.__name__
  newcall.__doc__ = v.__doc__ + "\nYou can also call with keyword arguments if the Function has a known scheme\nExample: nlp(x=x)\n"
  return newcall

for name,cl in locals().items():
  if not inspect.isclass(cl): continue
  for k,v in inspect.getmembers(cl, inspect.ismethod):
    if k == "__del__" or v.__name__ == "<lambda>": continue
    vv = v
    if k=="__call__" and issubclass(cl,Function):
      vv = swig_improvedcall(v)
    setattr(cl,k,swig_monkeypatch(vv))
  for k,v in inspect.getmembers(cl, inspect.isfunction):
    setattr(cl,k,staticmethod(swig_monkeypatch(v,cl=False)))
  
for name,v in locals().items():
  if not inspect.isfunction(v): continue
  if name.startswith("swig") : continue
  p = swig_monkeypatch(v,cl=False)
  #setattr(casadi,name,p)
  import sys
  setattr(sys.modules[__name__], name, p)


%}

#endif

// Cleanup for dependent modules
%exception {
  $action
}
