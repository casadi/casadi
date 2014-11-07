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

namespace casadi{

%{
#ifndef SWIGXML

/// casadi::Slice
template<> char meta< casadi::Slice >::expected_message[] = "Expecting Slice or number";
template <>
int meta< casadi::Slice >::as(GUESTOBJECT *p, casadi::Slice &m) {
  NATIVERETURN(casadi::Slice,m)
#ifdef SWIGPYTHON
  if (PyInt_Check(p)) {
    m.start_ = PyInt_AsLong(p);
    m.stop_ = m.start_+1;
    if (m.stop_==0) m.stop_ = std::numeric_limits<int>::max();
    return true;
  } else if (PySlice_Check(p)) {
    PySliceObject *r = (PySliceObject*)(p);
    m.start_ = (r->start == Py_None || PyInt_AsLong(r->start) < std::numeric_limits<int>::min()) ? std::numeric_limits<int>::min() : PyInt_AsLong(r->start);
    m.stop_  = (r->stop ==Py_None || PyInt_AsLong(r->stop)> std::numeric_limits<int>::max()) ? std::numeric_limits<int>::max() : PyInt_AsLong(r->stop) ;
    if(r->step !=Py_None) m.step_  = PyInt_AsLong(r->step);
    return true;
  } else {
    return false;
  }
#else
  return false;
#endif
}

/// casadi::IndexList
template<> char meta< casadi::IndexList >::expected_message[] = "Expecting Slice or number or list of ints";
template <>
int meta< casadi::IndexList >::as(GUESTOBJECT *p, casadi::IndexList &m) {
  NATIVERETURN(casadi::IndexList,m)
  if (meta< int >::couldbe(p)) {
    m.type = casadi::IndexList::INT;
    meta< int >::as(p,m.i);
  } else if (meta< std::vector<int> >::couldbe(p)) {
    m.type = casadi::IndexList::IVECTOR;
    return meta< std::vector<int> >::as(p,m.iv);
  } else if (meta< casadi::Slice>::couldbe(p)) {
    m.type = casadi::IndexList::SLICE;
    return meta< casadi::Slice >::as(p,m.slice);
  } else {
    return false;
  }
  return true;
}

#endif // SWIGXML



%}

%{
template<> swig_type_info** meta< casadi::Slice >::name = &SWIGTYPE_p_casadi__Slice;
template<> swig_type_info** meta< casadi::IndexList >::name = &SWIGTYPE_p_casadi__IndexList;
%}

%my_generic_const_typemap(PRECEDENCE_SLICE,casadi::Slice);
%my_generic_const_typemap(PRECEDENCE_IndexVector,casadi::IndexList);


#ifdef SWIGPYTHON

#ifdef WITH_NUMPY
/**

Accepts: 2D numpy.ndarray, numpy.matrix (contiguous, native byte order, datatype double)   - DENSE
         1D numpy.ndarray, numpy.matrix (contiguous, native byte order, datatype double)   - SPARSE
         2D scipy.csc_matrix
*/

%typemap(in,numinputs=1) (double * val,int len,int stride1, int stride2,SparsityType sp)  {
	PyObject* p = $input;
	$3 = 0;
	$4 = 0;
	if (is_array(p)) {
			if (!(array_is_native(p) && array_type(p)==NPY_DOUBLE))
			  SWIG_exception_fail(SWIG_TypeError, "Array should be native & of datatype double");
			  
	    if (!(array_is_contiguous(p))) {
	      if (PyArray_CHKFLAGS((PyArrayObject *) p,NPY_ALIGNED)) {
	        $3 = PyArray_STRIDE((PyArrayObject *) p,0)/sizeof(double);
	        $4 = PyArray_STRIDE((PyArrayObject *) p,1)/sizeof(double);
	      } else {
			   SWIG_exception_fail(SWIG_TypeError, "Array should be contiguous or aligned");
	      }
	    }
	    
			if (array_numdims(p)==2) {
				if (!(array_size(p,0)==arg1->size1() && array_size(p,1)==arg1->size2()) ) {
				  std::stringstream s;
				  s << "SWIG::typemap(in) (double *val,int len,SparsityType sp) " << std::endl;
				  s << "Array is not of correct shape.";
				  s << "Expecting shape (" << arg1->size1() << "," << arg1->size2() << ")" << ", but got shape (" << array_size(p,0) << "," << array_size(p,1) <<") instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			  }
			  $5 = casadi::DENSETRANS;
			  $2 = array_size(p,0)*array_size(p,1);
			  $1 = (double*) array_data(p);
			} else if (array_numdims(p)==1) {
				if (!(array_size(p,0)==arg1->size()) ) {
				  std::stringstream s;
				  s << "SWIG::typemap(in) (double *val,int len,SparsityType sp) " << std::endl;
				  s << "Array is not of correct size. Should match number of non-zero elements.";
				  s << "Expecting " << array_size(p,0) << " non-zeros, but got " << arg1->size() <<" instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			  }
			  $5 = casadi::SPARSE;
			  $2 = array_size(p,0);
			  $1 = (double*) array_data(p);
			} else {
			  SWIG_exception_fail(SWIG_TypeError, "Expecting 1D or 2D numpy.ndarray");
			}
	} else if (PyObjectHasClassName(p,"csc_matrix")) {
			$5 = casadi::SPARSE;
			PyObject * narray=PyObject_GetAttrString( p, "data"); // narray needs to be decref'ed
			if (!(array_is_contiguous(narray) && array_is_native(narray) && array_type(narray)==NPY_DOUBLE))
			  SWIG_exception_fail(SWIG_TypeError, "csc_matrix should be contiguous, native & of datatype double");
			$2 = array_size(narray,0);
			if (!(array_size(narray,0)==arg1->size() ) ) {
					std::stringstream s;
				  s << "SWIG::typemap(in) (double *val,int len,SparsityType sp) " << std::endl;
				  s << "csc_matrix does not have correct number of non-zero elements.";
				  s << "Expecting " << arg1->size() << " non-zeros, but got " << array_size(narray,0) << " instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
		      Py_DECREF(narray);
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			}
			$1 = (double*) array_data(narray);
			Py_DECREF(narray);
	} else {
			SWIG_exception_fail(SWIG_TypeError, "Unrecognised object");
	}
	
}

/**
Accepts: 2D numpy.ndarray, numpy.matrix (any setting of contiguous, native byte order, datatype)  - DENSE
         1D numpy.ndarray, numpy.matrix (any setting of contiguous, native byte order, datatype double) - SPARSE
         2D scipy.csc_matrix (any setting of contiguous, native byte order, datatype double) 
*/
%typemap(in,numinputs=1) (const double *val,int len,SparsityType sp) (PyArrayObject* array=NULL, int array_is_new_object=0)  {
	PyObject* p = $input;
	if (is_array(p)) {
			array = obj_to_array_contiguous_allow_conversion(p,NPY_DOUBLE,&array_is_new_object);
			if (array_numdims(array)==2) {
				if (!(array_size(array,0)==arg1->size1() && array_size(array,1)==arg1->size2()) ) {
				  std::stringstream s;
				  s << "SWIG::typemap(in) (const double *val,int len,SparsityType sp) " << std::endl;
				  s << "Array is not of correct shape.";
				  s << "Expecting shape (" << arg1->size1() << "," << arg1->size2() << ")" << ", but got shape (" << array_size(array,0) << "," << array_size(array,1) <<") instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			  }
			  $3 = casadi::DENSETRANS;
			  $2 = array_size(array,0)*array_size(array,1);
			  $1 = (double*) array_data(array);
			} else if (array_numdims(array)==1) {
				if (!(array_size(array,0)==arg1->size()) ) {
				  std::stringstream s;
				  s << "SWIG::typemap(in) (const double *val,int len,SparsityType sp) " << std::endl;
				  s << "Array is not of correct size. Should match number of non-zero elements.";
				  s << "Expecting " << arg1->size() << " non-zeros, but got " << array_size(array,0) << " instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			  }
			  $3 = casadi::SPARSE;
			  $2 = array_size(array,0);
			  $1 = (double*) array_data(array);
			} else {
			  SWIG_exception_fail(SWIG_TypeError, "Expecting 1D or 2D numpy.ndarray");
			}
	} else if (PyObjectHasClassName(p,"csc_matrix")) {
			$3 = casadi::SPARSE;
			PyObject * narray=PyObject_GetAttrString( p, "data"); // narray needs to be decref'ed
			$2 = array_size(narray,0);
			if (!(array_size(narray,0)==arg1->size() ) ) {
					std::stringstream s;
				  s << "SWIG::typemap(in) (const double *val,int len,SparsityType sp) " << std::endl;
				  s << "csc_matrix does not have correct number of non-zero elements.";
				  s << "Expecting " << arg1->size() << " non-zeros, but got " << array_size(narray,0) << " instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
          Py_DECREF(narray);
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			}
			array = obj_to_array_contiguous_allow_conversion(narray,NPY_DOUBLE,&array_is_new_object);
			$1 = (double*) array_data(array);
			Py_DECREF(narray);
	} else {
			SWIG_exception_fail(SWIG_TypeError, "Unrecognised object");
	}
	
}

%typemap(freearg) (const double *val,int len,SparsityType sp) {
    if (array_is_new_object$argnum && array$argnum) { Py_DECREF(array$argnum); }
}


%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) (double * val,int len,int stride1, int stride2,SparsityType sp) {
  PyObject* p = $input;
  if (((is_array(p) && array_numdims(p) < 3)  && array_type(p)!=NPY_OBJECT)|| PyObjectHasClassName(p,"csc_matrix")) {
    $1=1;
  } else {
    $1=0;
  }
}

%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) (const double * val,int len,SparsityType sp) {
  PyObject* p = $input;
  if (((is_array(p) && array_numdims(p) < 3)  && array_type(p)!=NPY_OBJECT)|| PyObjectHasClassName(p,"csc_matrix")) {
    $1=1;
  } else {
    $1=0;
  }
}
#endif // WITH_NUMPY
#endif // SWIGPYTHON


} // namespace casadi
