#ifdef SWIGOCTAVE
%rename(__el_mul__) __mul__;
%rename(__el_div__) __div__;
%rename(__rel_mul__) __rmul__;
%rename(__rel_div__) __rdiv__;
%rename(__el_pow__) __pow__;
%rename(__rel_pow__) __rpow__;
%rename(__mul__) prod;
%rename(__rmul__) rprod;
%rename(__transpose__) trans;
%rename(__div__) __mrdivide__;
%rename(__rdiv__) __rmrdivide__;
%rename(__pow__) __mpower__;
%rename(__rpow__) __rmpower__;
#endif // SWIGOCTAVE

%{
#include "casadi/matrix/crs_sparsity.hpp"
#include "casadi/matrix/matrix.hpp"
#include <sstream>
#include "casadi/casadi_exception.hpp"

// to allow for typechecking
#include "casadi/sx/sx.hpp"

// to typecheck for MX
#include "casadi/mx/mx.hpp"
// to have prod available
#include "casadi/mx/mx_tools.hpp"
%}

#ifdef SWIGPYTHON
%outputRefOwn(CasADi::CRSSparsity)
%outputRefOwn(std::vector< CasADi::SX >)
%outputRefOwn(std::vector< int >)
%outputRefOwn(std::vector< double >)
%outputRefOwn(CasADi::Matrix< double >)
%outputRefOwn(CasADi::Matrix< CasADi::SX >)
#endif // SWIGPYTHON

#ifdef SWIGPYTHON
%define %python_matrix_convertors
%pythoncode %{
        
    def toList(self):
        return list(self.data())
        
    def toMatrix(self):
        import numpy as n
        return n.matrix(self.toArray())
        
%}
%enddef 
%define %python_matrix_helpers(Type)
%pythoncode %{
    @property
    def shape(self):
        return (self.size1(),self.size2())
        
    def reshape(self,arg):
        return reshape(self,arg)
        
    @property
    def T(self):
        return trans(self)
        
    def __getitem__(self,s):
        if isinstance(s,tuple) and len(s)==2:
          return self.__Cgetitem__(s[0],s[1])  
        return self.__Cgetitem__(s)

    def __setitem__(self,s,val):
        if isinstance(s,tuple) and len(s)==2:
          return self.__Csetitem__(s[0],s[1],val)  
        return self.__Csetitem__(s,val)
     
%}
%enddef 
    
%define %python_array_wrappers(arraypriority)
%pythoncode %{

  __array_priority__ = arraypriority

  def __array_wrap__(self,out_arr,context=None):
    name = context[0].__name__
    args = list(context[1])
    
    if len(context[1])==3:
      raise Exception("Error with %s. Looks like you are using an assignment operator, such as 'a+=b' where 'a' is a numpy type. This is not supported, and cannot be supported without changing numpy." % name)

    if "vectorized" in name:
        name = name[:-len(" (vectorized)")]
    
    selfM = self
    
    if isinstance(self,SX):  # SX get's promoted to SXMatrix first
      selfM = SXMatrix(self)
    
    conversion = {"multiply": "mul", "divide": "div", "subtract":"sub","power":"pow"}
    if name in conversion:
      name = conversion[name]
    if len(context[1])==2 and context[1][1] is self and not(context[1][0] is self):
      name = 'r' + name
      args.reverse()
    if not(hasattr(selfM,name)):
      name = '__' + name + '__'
    fun=getattr(selfM, name)
    return fun(*args[1:])
     
     
  def __array__(self,*args,**kwargs):
    import numpy as n
    if len(args) > 1 and isinstance(args[1],tuple) and isinstance(args[1][0],n.ufunc):
      if len(args[1][1])==3:
        raise Exception("Error with %s. Looks like you are using an assignment operator, such as 'a+=b'. This is not supported when 'a' is a numpy type, and cannot be supported without changing numpy itself. Either upgrade a to a CasADi type first, or use 'a = a + b'. " % args[1][0].__name__)
      return n.array([1])
    else:
      if hasattr(self,'__array_custom__'):
        return self.__array_custom__(*args,**kwargs)
      else:
        return self.toArray()
      
%}
%enddef
#endif // SWIGPYTHON


#ifdef SWIGOCTAVE
%define %python_matrix_convertors
%enddef 
%define %python_matrix_helpers(Type)

  Type __hermitian__() const { return trans((*$self)); }
  
  std::vector<int> __dims__() const {
    std::vector<int> ret(2);
    ret[0] = $self->size1();
    ret[1] = $self->size2();
    return ret;
  }
  
%enddef 
#endif // SWIGOCTAVE


#ifdef SWIGOCTAVE
%rename(__paren__) indexed_one_based;
%rename(__paren__) indexed;
%rename(__paren_asgn__) indexed_one_based_assignment;
%rename(__paren_asgn__) indexed_assignment;
%rename(__vertcat__) vertcat;
%rename(__horzcat__) horzcat;
#endif
#ifdef SWIGPYTHON
%rename(__Cgetitem__) indexed_zero_based;
%rename(__Cgetitem__) indexed;
%rename(__Csetitem__) indexed_zero_based_assignment;
%rename(__Csetitem__) indexed_assignment;
#endif


#ifdef SWIGOCTAVE
namespace CasADi{
%extend Matrix<double> {
/// Create a 2D contiguous NP_DOUBLE numpy.ndarray

octave_value toSparse() {
  int nz = (*$self).size(), nr = (*$self).size1(), nc = (*$self).size2();
  
  Array<int> Ar(nz);
  Array<int> Ac(nz);
  
  std::vector<int> vr = (*$self).sparsity().getRow();
  Array<double> mydata(nz);
  const std::vector<double> &cdata = (*$self).data();
  
  for (int k=0;k<nz;k++) {
    Ar(k)=vr[k];
    Ac(k)=(*$self).sparsity().col()[k];
    mydata(k)=cdata[k];
  }
  
  return octave_value(SparseMatrix(mydata,Ar,Ac,nr,nc));
}

binopsrFull(CasADi::Matrix<double>)
binopsFull(const CasADi::Matrix<CasADi::SX> & b,,CasADi::Matrix<CasADi::SX>,CasADi::Matrix<CasADi::SX>)
binopsFull(const CasADi::MX & b,,CasADi::MX,CasADi::MX)

}; // extend Matrix<double>
} // namespace CasADi
#endif // SWIGOCTAVE


namespace CasADi{
%extend Matrix<double> {

void assign(const CasADi::Matrix<double>&rhs) { (*$self)=rhs; }
%python_matrix_convertors
%python_matrix_helpers(CasADi::Matrix<double>)

}
}


#ifdef SWIGPYTHON
namespace CasADi{
%extend Matrix<double> {
/// Create a 2D contiguous NP_DOUBLE numpy.ndarray

#ifdef WITH_NUMPY
PyObject* arrayView() {
  if ($self->size()!=$self->numel()) 
    throw  CasADi::CasadiException("Matrix<double>::arrayview() can only construct arrayviews for dense DMatrices.");
  npy_intp dims[2];
  dims[0] = $self->size1();
  dims[1] = $self->size2();
  std::vector<double> &v = $self->data();
  return PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, &v[0]);
}
#endif // WITH_NUMPY

%pythoncode %{
  def __eq__(self,other):
    return _casadi.__eq__(self,other)
%}
    
%pythoncode %{
  def toArray(self,shared=False):
    import numpy as n
    if shared:
      if self.size()!=self.numel():
        raise Expection("toArray(shared=True) only possible for dense arrays.")
      return self.arrayView()
    else:
      r = n.zeros((self.size1(),self.size2()))
      self.get(r)
    return r
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
  def toCsr_matrix(self):
    import numpy as n
    from scipy.sparse import csr_matrix
    return csr_matrix( (list(self.data()),self.sparsity().col(),self.sparsity().rowind()), shape = (self.size1(),self.size2()), dtype=n.double )
%}

%pythoncode %{
  def __float__(self):
    if self.numel()!=1:
      raise Exception("Only a scalar can be cast to a float")
    if self.size()==0:
      return 0.0
    return self.toScalar()
%}

%pythoncode %{
  def __int__(self):
    if self.numel()!=1:
      raise Exception("Only a scalar can be cast to a float")
    if self.size()==0:
      return 0
    return int(self.toScalar())
%}

%pythoncode %{
  def __nonzero__(self):
    if self.numel()!=1:
      raise Exception("Only a scalar can be cast to a float")
    if self.size()==0:
      return 0
    return self.toScalar()!=0
%}

%pythoncode %{
  def __abs__(self):
    return abs(self.__float__())
%}

binopsrFull(CasADi::Matrix<double>)
binopsFull(const CasADi::Matrix<CasADi::SX> & b,,CasADi::Matrix<CasADi::SX>,CasADi::Matrix<CasADi::SX>)
binopsFull(const CasADi::MX & b,,CasADi::MX,CasADi::MX)

}; // extend Matrix<double>
} // namespace CasADi
#endif // SWIGPYTHON

  


namespace CasADi{

%my_generic_const_typemap(PRECEDENCE_SLICE,CasADi::Slice);
%my_generic_const_typemap(PRECEDENCE_IndexVector,CasADi::IndexList);
%my_generic_const_typemap(PRECEDENCE_DMatrix,CasADi::Matrix<double>);

#ifdef SWIGPYTHON

#ifdef WITH_NUMPY
/**

Accepts: 2D numpy.ndarray, numpy.matrix (contiguous, native byte order, datatype double)   - DENSE
         1D numpy.ndarray, numpy.matrix (contiguous, native byte order, datatype double)   - SPARSE
         2D scipy.csr_matrix
*/

%typemap(in,numinputs=1) (double * val,int len,int stride1, int stride2,Sparsity sp)  {
	PyObject* p = $input;
	$3 = 0;
	$4 = 0;
	if (is_array(p)) {
			if (!(array_is_native(p) && array_type(p)==NPY_DOUBLE))
			  SWIG_exception_fail(SWIG_TypeError, "Array should be native & of datatype double");
			  
	    if (!(array_is_contiguous(p))) {
	      if (PyArray_CHKFLAGS(p,NPY_ALIGNED)) {
	        $3 = PyArray_STRIDE(p,0)/sizeof(double);
	        $4 = PyArray_STRIDE(p,1)/sizeof(double);
	      } else {
			   SWIG_exception_fail(SWIG_TypeError, "Array should be contiguous or aligned");
	      }
	    }
	    
			if (array_numdims(p)==2) {
				if (!(array_size(p,0)==arg1->size1() && array_size(p,1)==arg1->size2()) ) {
				  std::stringstream s;
				  s << "SWIG::typemap(in) (double *val,int len,Sparsity sp) " << std::endl;
				  s << "Array is not of correct shape.";
				  s << "Expecting shape (" << arg1->size1() << "," << arg1->size2() << ")" << ", but got shape (" << array_size(p,0) << "," << array_size(p,1) <<") instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			  }
			  $5 = CasADi::DENSE;
			  $2 = array_size(p,0)*array_size(p,1);
			  $1 = (double*) array_data(p);
			} else if (array_numdims(p)==1) {
				if (!(array_size(p,0)==arg1->size()) ) {
				  std::stringstream s;
				  s << "SWIG::typemap(in) (double *val,int len,Sparsity sp) " << std::endl;
				  s << "Array is not of correct size. Should match number of non-zero elements.";
				  s << "Expecting " << array_size(p,0) << " non-zeros, but got " << arg1->size() <<" instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			  }
			  $5 = CasADi::SPARSE;
			  $2 = array_size(p,0);
			  $1 = (double*) array_data(p);
			} else {
			  SWIG_exception_fail(SWIG_TypeError, "Expecting 1D or 2D numpy.ndarray");
			}
	} else if (PyObjectHasClassName(p,"csr_matrix")) {
			$5 = CasADi::SPARSE;
			PyObject * narray=PyObject_GetAttrString( p, "data"); // narray needs to be decref'ed
			if (!(array_is_contiguous(narray) && array_is_native(narray) && array_type(narray)==NPY_DOUBLE))
			  SWIG_exception_fail(SWIG_TypeError, "csr_matrix should be contiguous, native & of datatype double");
			$2 = array_size(narray,0);
			if (!(array_size(narray,0)==arg1->size() ) ) {
					std::stringstream s;
				  s << "SWIG::typemap(in) (double *val,int len,Sparsity sp) " << std::endl;
				  s << "csr_matrix does not have correct number of non-zero elements.";
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
         2D scipy.csr_matrix (any setting of contiguous, native byte order, datatype double) 
*/
%typemap(in,numinputs=1) (const double *val,int len,Sparsity sp) (PyArrayObject* array, int array_is_new_object=0)  {
	PyObject* p = $input;
	if (is_array(p)) {
			array = obj_to_array_contiguous_allow_conversion(p,NPY_DOUBLE,&array_is_new_object);
			if (array_numdims(array)==2) {
				if (!(array_size(array,0)==arg1->size1() && array_size(array,1)==arg1->size2()) ) {
				  std::stringstream s;
				  s << "SWIG::typemap(in) (const double *val,int len,Sparsity sp) " << std::endl;
				  s << "Array is not of correct shape.";
				  s << "Expecting shape (" << arg1->size1() << "," << arg1->size2() << ")" << ", but got shape (" << array_size(array,0) << "," << array_size(array,1) <<") instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			  }
			  $3 = CasADi::DENSE;
			  $2 = array_size(array,0)*array_size(array,1);
			  $1 = (double*) array_data(array);
			} else if (array_numdims(array)==1) {
				if (!(array_size(array,0)==arg1->size()) ) {
				  std::stringstream s;
				  s << "SWIG::typemap(in) (const double *val,int len,Sparsity sp) " << std::endl;
				  s << "Array is not of correct size. Should match number of non-zero elements.";
				  s << "Expecting " << arg1->size() << " non-zeros, but got " << array_size(array,0) << " instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			  }
			  $3 = CasADi::SPARSE;
			  $2 = array_size(array,0);
			  $1 = (double*) array_data(array);
			} else {
			  SWIG_exception_fail(SWIG_TypeError, "Expecting 1D or 2D numpy.ndarray");
			}
	} else if (PyObjectHasClassName(p,"csr_matrix")) {
			$3 = CasADi::SPARSE;
			PyObject * narray=PyObject_GetAttrString( p, "data"); // narray needs to be decref'ed
			$2 = array_size(narray,0);
			if (!(array_size(narray,0)==arg1->size() ) ) {
					std::stringstream s;
				  s << "SWIG::typemap(in) (const double *val,int len,Sparsity sp) " << std::endl;
				  s << "csr_matrix does not have correct number of non-zero elements.";
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

%typemap(freearg) (const double *val,int len,Sparsity sp) {
    if (array_is_new_object$argnum && array$argnum) { Py_DECREF(array$argnum); }
}


%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) (double * val,int len,int stride1, int stride2,Sparsity sp) {
  PyObject* p = $input;
  if ((is_array(p) && array_numdims(p) < 3)  && array_type(p)!=NPY_OBJECT|| PyObjectHasClassName(p,"csr_matrix")) {
    $1=1;
  } else {
    $1=0;
  }
}

%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) (const double * val,int len,Sparsity sp) {
  PyObject* p = $input;
  if ((is_array(p) && array_numdims(p) < 3)  && array_type(p)!=NPY_OBJECT|| PyObjectHasClassName(p,"csr_matrix")) {
    $1=1;
  } else {
    $1=0;
  }
}
#endif // WITH_NUMPY
#endif // SWIGPYTHON


} // namespace CasADi
