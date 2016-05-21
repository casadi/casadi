#include "stdio.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <arrayobject.h>

#ifndef NDARRAY_VERSION
#define NDARRAY_VERSION NPY_VERSION
#endif

#ifndef NPY_ALIGNED
#define NPY_ALIGNED NPY_ARRAY_ALIGNED
#endif

#define is_array(a)            ((a) && PyArray_Check((PyArrayObject *)a))
#define array_type(a)          (int)(PyArray_TYPE((PyArrayObject *)a))
#define array_is_contiguous(a) (PyArray_ISCONTIGUOUS((PyArrayObject *)a))
#define array_is_native(a)     (PyArray_ISNOTSWAPPED((PyArrayObject *)a))
#define array_numdims(a)    (PyArray_NDIM(((PyArrayObject *)a)))
#define array_dimensions(a) (PyArray_DIMS(((PyArrayObject *)a)))
#define array_size(a,i)     (PyArray_DIM(((PyArrayObject *)a),i))
#define array_data(a)       (PyArray_DATA(((PyArrayObject *)a)))

#if PY_MAJOR_VERSION >= 3
#define PyInt_Type PyLong_Type
#endif 

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
  PyErr_Clear();
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
  PyErr_Clear();
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
