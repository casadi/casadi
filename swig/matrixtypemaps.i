%{
#include "casadi/matrix/crs_sparsity.hpp"
#include "casadi/matrix/matrix.hpp"
%}

%include "typemaps.i"

namespace CasADi{
%extend Matrix<double> {
/// Create a 2D contiguous NP_DOUBLE numpy.ndarray
%pythoncode %{
  def toArray(self):
    import numpy as n
    r = n.zeros((self.size1(),self.size2()))
    self.get(r)
    return r
%}
// Can this be done in C?
%pythoncode %{
  def toMatrix(self):
    import numpy as n
    return n.matrix(self.getArray())
%}

%pythoncode %{
  def toList(self):
    return list(self)
%}

%pythoncode %{
  def toTuple(self):
    return tuple(self)
%}

// Can this be done in C?
%pythoncode %{
  def toCsr_matrix(self):
    import numpy as n
    from scipy.sparse import csr_matrix
    return csr_matrix( (self,self.sparsity().col(),self.sparsity().rowind()), shape = (self.size1(),self.size2()), dtype=n.double )
%}
};
}

/** If the array is of type double, contiguous and in native byte order, this function is efficient.
* Other types of numpy array will trigger conversion, requiring temporary allocation of memory.
*/
#ifdef WITH_NUMPY
%inline %{
bool PyObjectHasClassName(PyObject* p, const char * name) {
	return strcmp(PyString_AsString(PyObject_GetAttrString(PyObject_GetAttrString( p, "__class__"),"__name__")),name)==0;
}

CasADi::Matrix<double> * typemapDMatrixHelper(PyObject* p, PyArrayObject* array, int& array_is_new_object, bool& freearg) {
if (is_array(p)) { // Numpy arrays will be cast to dense Matrix<double>
		if (array_numdims(p)>2 || array_numdims(p)<1)
			SWIG_Error(SWIG_TypeError, "asMatrixDouble: Number of dimensions must be 1 or 2.");
		int nrows = array_size(p,0); // 1D array is cast into column vector
		int ncols  = 1;
		if (array_numdims(p)==2)
			ncols=array_size(p,1); 
		int size=nrows*ncols; // number of elements in the dense matrix
		if (!array_is_native(p)) 
			SWIG_Error(SWIG_TypeError, "asMatrixDouble: array byte order should be native.");
		// Make sure we have a contigous array with double datatype
		array = obj_to_array_contiguous_allow_conversion(p,NPY_DOUBLE,&array_is_new_object);
		double* d=(double*) array->data;
		std::vector<double> v(d,d+size);

		// Construct a dense Matrix<double>
		// Memory will be freed in typemap(freearg)
		return new CasADi::Matrix<double>(v, nrows, ncols);

	} else if(PyObjectHasClassName(p,"csr_matrix")) { // scipy's csr_matrix will be cast to sparse Matrix<double>
		PyObject * narray=PyObject_GetAttrString( p, "data");
		if (!(is_array(narray) && array_numdims(narray)==1))
			SWIG_Error(SWIG_TypeError, "asMatrixDouble: data should be numpy array");
		array = obj_to_array_contiguous_allow_conversion(narray,NPY_DOUBLE,&array_is_new_object);
		int size=array_size(narray,0); // number on non-zeros

		double* d=(double*) array->data;
		std::vector<double> v(d,d+size);

		// Get the dimensions of the csr_matrix
		PyObject * shape = PyObject_GetAttrString( p, "shape");
		int nrows=PyInt_AsLong(PyTuple_GetItem(shape,0));
		int ncols=PyInt_AsLong(PyTuple_GetItem(shape,1));
		
		// Construct the 'col' vector needed for initialising the correct sparsity
		PyObject * col = PyObject_GetAttrString(p,"indices");
		if (!(is_array(col) && array_numdims(col)==1 && array_type(col)==NPY_INT))
			SWIG_Error(SWIG_TypeError, "asMatrixDouble: data.indices should be numpy array");
		int* cold=(int*) array_data(col);
		std::vector<int> colv(cold,cold+size);
			
		// Construct the 'rowind' vector needed for initialising the correct sparsity
		PyObject * rowind = PyObject_GetAttrString(p,"indptr");
		if (!(is_array(rowind) && array_numdims(rowind)==1 && array_type(rowind)==NPY_INT))
			SWIG_Error(SWIG_TypeError, "asMatrixDouble: data.indptr should be numpy array");
		int* rowindd=(int*) array_data(rowind);
		std::vector<int> rowindv(rowindd,rowindd+(nrows+1));
			
		return new CasADi::Matrix<double>(nrows,ncols,colv,rowindv, v);
	} else if (PyObjectHasClassName(p,"DMatrix")) { // Dmatrix object get passed on as-is.
		void *pd = 0 ;
		int res = SWIG_ConvertPtr(p, &pd,SWIGTYPE_p_CasADi__MatrixT_double_t, 0 |  0 );
		if (!SWIG_IsOK(res)) {
			SWIG_Error(SWIG_TypeError, "asMatrixDouble: DMatrix cast problem");
		}
		freearg = false;
		return reinterpret_cast< CasADi::Matrix< double > * >(pd);
	} else {
		std::cerr << "type:" << PyString_AsString(PyObject_GetAttrString(PyObject_GetAttrString( p, "__class__"),"__name__")) << std::endl;
		SWIG_Error(SWIG_TypeError, "asMatrixDouble: unknown type");
	}
}




%}
#endif


#ifdef WITH_NUMPY
namespace CasADi{
%typemap(in) const Matrix<double> &  (PyArrayObject* array=NULL, int array_is_new_object, bool freearg = true){
	$1 = typemapDMatrixHelper($input,array,array_is_new_object, freearg);
}


%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) const Matrix<double> & {
    PyObject* p = $input;
    if ((is_array(p) && array_numdims(p)==2) || PyObjectHasClassName(p,"csr_matrix") || PyObjectHasClassName(p,"DMatrix")) {
	$1=1;
    } else {
	$1=0;
    }
}


%typemap(freearg) const Matrix<double> & {
    if (freearg$argnum) {
	    if ($1) 
		delete $1;
    }
    if (array_is_new_object$argnum && array$argnum) { Py_DECREF(array$argnum); }
}




/**
Accepts: 2D numpy.ndarray, numpy.matrix (contiguous, native byte order, datatype double)   - DENSE
         1D numpy.ndarray, numpy.matrix (contiguous, native byte order, datatype double)   - SPARSE
         2D scipy.csr_matrix
*/

%typemap(in,numinputs=1) (double *val,int len,Sparsity sp)  {
	PyObject* p = $input;
	if (is_array(p)) {
			if (!(array_is_contiguous(p) && array_is_native(p) && array_type(p)==NPY_DOUBLE))
			  SWIG_exception_fail(SWIG_TypeError, "Array should be contiguous, native & of datatype double");
			if (array_numdims(p)==2) {
				if (!(array_size(p,0)==arg1->size1() && array_size(p,1)==arg1->size2()) )
			    SWIG_exception_fail(SWIG_TypeError, "Array is not of correct shape.");
			  $3 = CasADi::DENSE;
			  $2 = array_size(p,0)*array_size(p,1);
			  $1 = (double*) array_data(p);
			} else if (array_numdims(p)==1) {
				if (!(array_size(p,0)==arg1->size()) )
			    SWIG_exception_fail(SWIG_TypeError, "Array is not of correct size. Should match number of non-zero elements.");
			  $3 = CasADi::SPARSE;
			  $2 = array_size(p,0);
			  $1 = (double*) array_data(p);
			} else {
			  SWIG_exception_fail(SWIG_TypeError, "Expecting 1D or 2D numpy.ndarray");
			}
	} else if (PyObjectHasClassName(p,"csr_matrix")) {
			$3 = CasADi::SPARSE;
			PyObject * narray=PyObject_GetAttrString( p, "data");
			if (!(array_is_contiguous(narray) && array_is_native(narray) && array_type(narray)==NPY_DOUBLE))
			  SWIG_exception_fail(SWIG_TypeError, "csr_matrix should be contiguous, native & of datatype double");
			$2 = array_size(narray,0);
			if (!(array_size(narray,0)==arg1->size() ) )
			  SWIG_exception_fail(SWIG_TypeError, "csr_matrix does not have correct number of non-zero elements");
			$1 = (double*) array_data(narray);
	} else {
			SWIG_exception_fail(SWIG_TypeError, "Unrecognised object");
	}
	
}

/**
Accepts: 2D numpy.ndarray, numpy.matrix (any setting of contiguous, native byte order, datatype)  - DENSE
         1D numpy.ndarray, numpy.matrix (any setting of contiguous, native byte order, datatype double) - SPARSE
         2D scipy.csr_matrix (any setting of contiguous, native byte order, datatype double) 
*/
%typemap(in,numinputs=1) (const double *val,int len,Sparsity sp) (PyArrayObject* array, int array_is_new_object=0)  {
	PyObject* p = $input;
	if (is_array(p)) {
			array = obj_to_array_contiguous_allow_conversion(p,NPY_DOUBLE,&array_is_new_object);
			if (array_numdims(p)==2) {
				if (!(array_size(p,0)==arg1->size1() && array_size(p,1)==arg1->size2()) )
			    SWIG_exception_fail(SWIG_TypeError, "Array is not of correct shape.");
			  $3 = CasADi::DENSE;
			  $2 = array_size(p,0)*array_size(p,1);
			  $1 = (double*) array_data(p);
			} else if (array_numdims(p)==1) {
				if (!(array_size(p,0)==arg1->size()) )
			    SWIG_exception_fail(SWIG_TypeError, "Array is not of correct size. Should match number of non-zero elements.");
			  $3 = CasADi::SPARSE;
			  $2 = array_size(p,0);
			  $1 = (double*) array_data(p);
			} else {
			  SWIG_exception_fail(SWIG_TypeError, "Expecting 1D or 2D numpy.ndarray");
			}
	} else if (PyObjectHasClassName(p,"csr_matrix")) {
			$3 = CasADi::SPARSE;
			PyObject * narray=PyObject_GetAttrString( p, "data");
			$2 = array_size(narray,0);
			if (!(array_size(narray,0)==arg1->size() ) )
			  SWIG_exception_fail(SWIG_TypeError, "csr_matrix does not have correct number of non-zero elements");
			array = obj_to_array_contiguous_allow_conversion(narray,NPY_DOUBLE,&array_is_new_object);
			$1 = (double*) array_data(array);
	} else {
			SWIG_exception_fail(SWIG_TypeError, "Unrecognised object");
	}
	
}

%typemap(freearg) (const double *val,int len,Sparsity sp) {
    if (array_is_new_object$argnum && array$argnum) { Py_DECREF(array$argnum); }
}


%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) (double * val,int len,Sparsity sp) {
    PyObject* p = $input;
    if ((is_array(p) && array_numdims(p) < 3 )|| PyObjectHasClassName(p,"csr_matrix")) {
	$1=1;
    } else {
	$1=0;
    }
}

%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) (const double * val,int len,Sparsity sp) {
    PyObject* p = $input;
    if ((is_array(p) && array_numdims(p) < 3 ) || PyObjectHasClassName(p,"csr_matrix")) {
	$1=1;
    } else {
	$1=0;
    }
}

}
#endif



