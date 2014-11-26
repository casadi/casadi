#include "stdio.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <arrayobject.h>

/* The following code originally appeared in
 * enthought/kiva/agg/src/numeric.i written by Eric Jones.  It was
 * translated from C++ to C by John Hunter.  Bill Spotz has modified
 * it slightly to fix some minor bugs, upgrade to numpy (all
 * versions), add some comments and some functionality.
 */

#ifndef NDARRAY_VERSION
#define NDARRAY_VERSION NPY_VERSION
#endif

#ifndef NPY_ALIGNED
#define NPY_ALIGNED NPY_ARRAY_ALIGNED
#endif


#ifndef NPY_ARRAY_F_CONTIGUOUS
#define NPY_ARRAY_F_CONTIGUOUS NPY_F_CONTIGUOUS
#endif

#ifndef NPY_ARRAY_FARRAY
#define NPY_ARRAY_FARRAY NPY_FARRAY
#endif

#ifndef NPY_ARRAY_CARRAY
#define NPY_ARRAY_CARRAY NPY_CARRAY
#endif

#ifndef NPY_ARRAY_DEFAULT
#define NPY_ARRAY_DEFAULT NPY_DEFAULT
#endif

#ifndef NPY_ARRAY_ENSUREARRAY
#define NPY_ARRAY_ENSUREARRAY NPY_ENSUREARRAY
#endif
/* Macros to extract array attributes.
 */
#define is_array(a)            ((a) && PyArray_Check((PyArrayObject *)a))
#define array_type(a)          (int)(PyArray_TYPE((PyArrayObject *)a))
#define array_is_contiguous(a) (PyArray_ISCONTIGUOUS((PyArrayObject *)a))
#define array_is_native(a)     (PyArray_ISNOTSWAPPED((PyArrayObject *)a))

#if NDARRAY_VERSION < 0x01000009
#define array_numdims(a)       (((PyArrayObject *)a)->nd)
#define array_dimensions(a)    (((PyArrayObject *)a)->dimensions)
#define array_size(a,i)        (((PyArrayObject *)a)->dimensions[i])
#define array_data(a)          (((PyArrayObject *)a)->data)
#ifndef PyArray_CHKFLAGS
#define PyArray_CHKFLAGS(m, FLAGS) \
	((((PyArrayObject *)(m))->flags & (FLAGS)) == (FLAGS))
#endif
#else
#define array_numdims(a)    (PyArray_NDIM(((PyArrayObject *)a)))
#define array_dimensions(a) (PyArray_DIMS(((PyArrayObject *)a)))
#define array_size(a,i)     (PyArray_DIM(((PyArrayObject *)a),i))
#define array_data(a)       (PyArray_DATA(((PyArrayObject *)a)))
#endif
	
/* Support older NumPy data type names
*/
#if NDARRAY_VERSION < 0x01000000
#define NPY_BOOL        PyArray_BOOL
#define NPY_BYTE        PyArray_BYTE
#define NPY_UBYTE       PyArray_UBYTE
#define NPY_SHORT       PyArray_SHORT
#define NPY_USHORT      PyArray_USHORT
#define NPY_INT         PyArray_INT
#define NPY_UINT        PyArray_UINT
#define NPY_LONG        PyArray_LONG
#define NPY_ULONG       PyArray_ULONG
#define NPY_LONGLONG    PyArray_LONGLONG
#define NPY_ULONGLONG   PyArray_ULONGLONG
#define NPY_FLOAT       PyArray_FLOAT
#define NPY_DOUBLE      PyArray_DOUBLE
#define NPY_LONGDOUBLE  PyArray_LONGDOUBLE
#define NPY_CFLOAT      PyArray_CFLOAT
#define NPY_CDOUBLE     PyArray_CDOUBLE
#define NPY_CLONGDOUBLE PyArray_CLONGDOUBLE
#define NPY_OBJECT      PyArray_OBJECT
#define NPY_STRING      PyArray_STRING
#define NPY_UNICODE     PyArray_UNICODE
#define NPY_VOID        PyArray_VOID
#define NPY_NTYPES      PyArray_NTYPES
#define NPY_NOTYPE      PyArray_NOTYPE
#define NPY_CHAR        PyArray_CHAR
#define NPY_USERDEF     PyArray_USERDEF
#define npy_intp        intp
#endif

/* Given a PyObject, return a string describing its type.
 */
const char* pytype_string(PyObject* py_obj) {
  if (py_obj == NULL          ) return "C NULL value";
  if (py_obj == Py_None       ) return "Python None" ;
  if (PyCallable_Check(py_obj)) return "callable"    ;
  if (PyString_Check(  py_obj)) return "string"      ;
  if (PyInt_Check(     py_obj)) return "int"         ;
  if (PyFloat_Check(   py_obj)) return "float"       ;
  if (PyDict_Check(    py_obj)) return "dict"        ;
  if (PyList_Check(    py_obj)) return "list"        ;
  if (PyTuple_Check(   py_obj)) return "tuple"       ;
  if (PyFile_Check(    py_obj)) return "file"        ;
  if (PyModule_Check(  py_obj)) return "module"      ;
  if (PyInstance_Check(py_obj)) return "instance"    ;

  return "unknown type";
}

/* Given a NumPy typecode, return a string describing the type.
 */
const char* typecode_string(int typecode) {
  static const char* type_names[25] = {"bool", "byte", "unsigned byte",
				 "short", "unsigned short", "int",
				 "unsigned int", "long", "unsigned long",
				 "long long", "unsigned long long",
				 "float", "double", "long double",
				 "complex float", "complex double",
				 "complex long double", "object",
				 "string", "unicode", "void", "ntypes",
				 "notype", "char", "unknown"};
  return typecode < 24 ? type_names[typecode] : type_names[24];
}

/* Make sure input has correct numpy type.  Allow character and byte
 * to match.  Also allow int and long to match.  This is deprecated.
 * You should use PyArray_EquivTypenums() instead.
 */
int type_match(int actual_type, int desired_type) {
  return PyArray_EquivTypenums(actual_type, desired_type);
}

/* Given a PyObject pointer, cast it to a PyArrayObject pointer if
 * legal.  If not, set the python error string appropriately and
 * return NULL.
 */
PyArrayObject* obj_to_array_no_conversion(PyObject* input, int typecode) {
  PyArrayObject* ary = NULL;
  if (is_array(input) && (typecode == NPY_NOTYPE ||
			  PyArray_EquivTypenums(array_type(input), typecode))) {
    ary = (PyArrayObject*) input;
  }
  else if is_array(input) {
    const char* desired_type = typecode_string(typecode);
    const char* actual_type  = typecode_string(array_type(input));
    PyErr_Format(PyExc_TypeError, 
		 "Array of type '%s' required.  Array of type '%s' given", 
		 desired_type, actual_type);
    ary = NULL;
  }
  else {
    const char * desired_type = typecode_string(typecode);
    const char * actual_type  = pytype_string(input);
    PyErr_Format(PyExc_TypeError, 
		 "Array of type '%s' required.  A '%s' was given", 
		 desired_type, actual_type);
    ary = NULL;
  }
  return ary;
}

/* Convert the given PyObject to a NumPy array with the given
 * typecode.  On success, return a valid PyArrayObject* with the
 * correct type.  On failure, the python error string will be set and
 * the routine returns NULL.
 */
PyArrayObject* obj_to_array_allow_conversion(PyObject* input, int typecode,
                                             int* is_new_object) {
  PyArrayObject* ary = NULL;
  PyObject* py_obj;
  if (is_array(input) && (typecode == NPY_NOTYPE ||
			  PyArray_EquivTypenums(array_type(input),typecode))) {
    ary = (PyArrayObject*) input;
    *is_new_object = 0;
  }
  else {
    py_obj = PyArray_FromObject(input, typecode, 0, 0);
    /* If NULL, PyArray_FromObject will have set python error value.*/
    ary = (PyArrayObject*) py_obj;
    *is_new_object = 1;
  }
  return ary;
}

/* Given a PyArrayObject, check to see if it is contiguous.  If so,
 * return the input pointer and flag it as not a new object.  If it is
 * not contiguous, create a new PyArrayObject using the original data,
 * flag it as a new object and return the pointer.
 */
PyArrayObject* make_contiguous(PyArrayObject* ary, int* is_new_object,
                               int min_dims, int max_dims) {
  PyArrayObject* result;
  if (array_is_contiguous(ary)) {
    result = ary;
    *is_new_object = 0;
  }
  else {
    result = (PyArrayObject*) PyArray_ContiguousFromObject((PyObject*)ary, 
							   array_type(ary), 
							   min_dims,
							   max_dims);
    *is_new_object = 1;
  }
  return result;
}

/* Convert a given PyObject to a contiguous PyArrayObject of the
 * specified type.  If the input object is not a contiguous
 * PyArrayObject, a new one will be created and the new object flag
 * will be set.
 */
PyArrayObject* obj_to_array_contiguous_allow_conversion(PyObject* input,
                                                        int typecode,
                                                        int* is_new_object) {
  int is_new1 = 0;
  int is_new2 = 0;
  PyArrayObject* ary2;
  PyArrayObject* ary1 = obj_to_array_allow_conversion(input, typecode, 
						      &is_new1);
  if (ary1) {
    ary2 = make_contiguous(ary1, &is_new2, 0, 0);
    if ( is_new1 && is_new2) {
      Py_DECREF(ary1);
    }
    ary1 = ary2;    
  }
  *is_new_object = is_new1 || is_new2;
  return ary1;
}

/* Test whether a python object is contiguous.  If array is
 * contiguous, return 1.  Otherwise, set the python error string and
 * return 0.
 */
int require_contiguous(PyArrayObject* ary) {
  int contiguous = 1;
  if (!array_is_contiguous(ary)) {
    PyErr_SetString(PyExc_TypeError,
		    "Array must be contiguous.  A non-contiguous array was given");
    contiguous = 0;
  }
  return contiguous;
}

/* Require that a numpy array is not byte-swapped.  If the array is
 * not byte-swapped, return 1.  Otherwise, set the python error string
 * and return 0.
 */
int require_native(PyArrayObject* ary) {
  int native = 1;
  if (!array_is_native(ary)) {
    PyErr_SetString(PyExc_TypeError,
		    "Array must have native byteorder.  A byte-swapped array was given");
    native = 0;
  }
  return native;
}

/* Require the given PyArrayObject to have a specified number of
 * dimensions.  If the array has the specified number of dimensions,
 * return 1.  Otherwise, set the python error string and return 0.
 */
int require_dimensions(PyArrayObject* ary, int exact_dimensions) {
  int success = 1;
  if (array_numdims(ary) != exact_dimensions) {
    PyErr_Format(PyExc_TypeError, 
		 "Array must have %d dimensions.  Given array has %d dimensions", 
		 exact_dimensions, array_numdims(ary));
    success = 0;
  }
  return success;
}

/* Require the given PyArrayObject to have one of a list of specified
 * number of dimensions.  If the array has one of the specified number
 * of dimensions, return 1.  Otherwise, set the python error string
 * and return 0.
 */
int require_dimensions_n(PyArrayObject* ary, int* exact_dimensions, int n) {
  int success = 0;
  int i;
  char dims_str[255] = "";
  char s[255];
  for (i = 0; i < n && !success; i++) {
    if (array_numdims(ary) == exact_dimensions[i]) {
      success = 1;
    }
  }
  if (!success) {
    for (i = 0; i < n-1; i++) {
      sprintf(s, "%d, ", exact_dimensions[i]);                
      strcat(dims_str,s);
    }
    sprintf(s, " or %d", exact_dimensions[n-1]);            
    strcat(dims_str,s);
    PyErr_Format(PyExc_TypeError, 
		 "Array must be have %s dimensions.  Given array has %d dimensions",
		 dims_str, array_numdims(ary));
  }
  return success;
}    

/* Require the given PyArrayObject to have a specified shape.  If the
 * array has the specified shape, return 1.  Otherwise, set the python
 * error string and return 0.
 */
int require_size(PyArrayObject* ary, npy_intp* size, int n) {
  int i;
  int success = 1;
  int len;
  char desired_dims[255] = "[";
  char s[255];
  char actual_dims[255] = "[";
  for(i=0; i < n;i++) {
    if (size[i] != -1 &&  size[i] != array_size(ary,i)) {
      success = 0;    
    }
  }
  if (!success) {
    for (i = 0; i < n; i++) {
      if (size[i] == -1) {
	sprintf(s, "*,");                
      }
      else
      {
	sprintf(s, "%ld,", (long int) size[i]);                
      }    
      strcat(desired_dims,s);
    }
    len = strlen(desired_dims);
    desired_dims[len-1] = ']';
    for (i = 0; i < n; i++) {
      sprintf(s, "%ld,", (long int) array_size(ary,i));                            
      strcat(actual_dims,s);
    }
    len = strlen(actual_dims);
    actual_dims[len-1] = ']';
    PyErr_Format(PyExc_TypeError, 
		 "Array must be have shape of %s.  Given array has shape of %s",
		 desired_dims, actual_dims);
  }
  return success;
}

/* End John Hunter translation (with modifications by Bill Spotz)
 */
