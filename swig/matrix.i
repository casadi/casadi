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

%inline %{
template<> swig_type_info** meta< CasADi::IndexList >::name = &SWIGTYPE_p_CasADi__IndexList;
template<> swig_type_info** meta< CasADi::Matrix<double> >::name = &SWIGTYPE_p_CasADi__MatrixT_double_t;
template<> swig_type_info** meta< std::vector< CasADi::Matrix<double> > >::name = &SWIGTYPE_p_std__vectorT_CasADi__MatrixT_double_t_std__allocatorT_CasADi__MatrixT_double_t_t_t;
template<> swig_type_info** meta< CasADi::Slice >::name = &SWIGTYPE_p_CasADi__Slice;
template<> swig_type_info** meta< CasADi::CRSSparsity >::name = &SWIGTYPE_p_CasADi__CRSSparsity;
%}

%outputConstRefCopy(CasADi::CRSSparsity)


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
#endif // SWIGPYTHON


#ifdef SWIGOCTAVE
%define concat_operator(ReturnType,LType,RType)
  // __concat__ should fill in the rhs into the lhs with offset (offsetx,offsety)
  // lhs is assumed to be equal in size or larger than rhs
  ReturnType __concat__(const LType &lhs,const RType &rhs, int offsetx, int offsety) const {
    ReturnType ret = lhs;
    ret(CasADi::range(offsetx,offsetx+rhs.size1()),CasADi::range(offsety,offsety+rhs.size2())) = rhs;
    return ret;
  }
%enddef
%define %python_matrix_convertors
%enddef 
%define %python_matrix_helpers(Type)

  int ndims() const {
    return 2;
  }
  
  int length() const {
    return $self->size1()>$self->size2() ? $self->size1(): $self->size2() ;
  }
  
  Type __hermitian__() const { return trans((*$self)); }
  
  
  Type __resize__(int nrows, int ncols) const {
    return Type(nrows,ncols);
  }
  
  concat_operator(Type,Type,Type)
  
  std::vector<int> __dims__() const {
    std::vector<int> ret(2);
    ret[0] = $self->size1();
    ret[1] = $self->size2();
    return ret;
  }
  
%enddef 
#endif // SWIGOCTAVE

/// CasADi::Matrix<double>
#ifdef SWIGPYTHON
%inline %{
template<> char meta< CasADi::Matrix<double> >::expected_message[] = "Expecting numpy.array2D, numpy.matrix, csr_matrix, DMatrix";

template <>
int meta< CasADi::Matrix<double> >::as(PyObject * p,CasADi::Matrix<double> &m) {
  NATIVERETURN(CasADi::Matrix<double>,m)
  if (is_array(p)) { // Numpy arrays will be cast to dense Matrix<double>
    if (array_numdims(p)>2 || array_numdims(p)<1) {
      SWIG_Error(SWIG_TypeError, "asMatrixDouble: Number of dimensions must be 1 or 2.");
      std::stringstream s;
      s << "SWIG::typemapDMatrixHelper:";
      s << "Number of dimensions must be 1 or 2.";
      s << "Got " << array_numdims(p) << " instead.";
      const std::string tmp(s.str());
      const char* cstr = tmp.c_str();
      SWIG_Error(SWIG_TypeError,  cstr);
    }
    int nrows = array_size(p,0); // 1D array is cast into column vector
    int ncols  = 1;
    if (array_numdims(p)==2)
      ncols=array_size(p,1); 
    int size=nrows*ncols; // number of elements in the dense matrix
    if (!array_is_native(p)) 
      SWIG_Error(SWIG_TypeError, "asMatrixDouble: array byte order should be native.");
    // Make sure we have a contigous array with double datatype
    int array_is_new_object;
    PyArrayObject* array = obj_to_array_contiguous_allow_conversion(p,NPY_DOUBLE,&array_is_new_object);

    double* d=(double*) array->data;
    std::vector<double> v(d,d+size);
    
    m = CasADi::Matrix<double>(v, nrows, ncols);
                  
    // Free memory
    if (array_is_new_object)
      Py_DECREF(array); 
  } else if(PyObjectHasClassName(p,"csr_matrix")) { // scipy's csr_matrix will be cast to sparse Matrix<double>
    PyObject * narray=PyObject_GetAttrString( p, "data"); // need's to be decref'ed
    if (!(is_array(narray) && array_numdims(narray)==1))
      SWIG_Error(SWIG_TypeError, "asMatrixDouble: data should be numpy array");
    int array_is_new_object;
    PyArrayObject* array = obj_to_array_contiguous_allow_conversion(narray,NPY_DOUBLE,&array_is_new_object);
    int size=array_size(array,0); // number on non-zeros
    double* d=(double*) array->data;
    std::vector<double> v(d,d+size);

    // Get the dimensions of the csr_matrix
    PyObject * shape = PyObject_GetAttrString( p, "shape"); // need's to be decref'ed
    int nrows=PyInt_AsLong(PyTuple_GetItem(shape,0));
    int ncols=PyInt_AsLong(PyTuple_GetItem(shape,1));
		
    // Construct the 'col' vector needed for initialising the correct sparsity
    PyObject * col = PyObject_GetAttrString(p,"indices"); // need's to be decref'ed
    if (!(is_array(col) && array_numdims(col)==1 && array_type(col)==NPY_INT))
      SWIG_Error(SWIG_TypeError, "asMatrixDouble: data.indices should be numpy array");
    int* cold=(int*) array_data(col);
    std::vector<int> colv(cold,cold+size);
    
    // Construct the 'rowind' vector needed for initialising the correct sparsity
    PyObject * rowind = PyObject_GetAttrString(p,"indptr"); // need's to be decref'ed
    if (!(is_array(rowind) && array_numdims(rowind)==1 && array_type(rowind)==NPY_INT))
      SWIG_Error(SWIG_TypeError, "asMatrixDouble: data.indptr should be numpy array");
    int* rowindd=(int*) array_data(rowind);
    std::vector<int> rowindv(rowindd,rowindd+(nrows+1));
    
    m = CasADi::Matrix<double>(nrows,ncols,colv,rowindv, v);
    
    Py_DECREF(narray);Py_DECREF(shape);Py_DECREF(col);Py_DECREF(rowind);
    
    if (array_is_new_object)
      Py_DECREF(array);
  } else if (meta< double >::couldbe(p)) {
    double t;
    int res = meta< double >::as(p,t);
    m = t;
    return res;
  } else if ( meta< double >::couldbe_sequence(p)) {
    std::vector <double> t;
    int res = meta< double >::as_vector(p,t);
    m = CasADi::Matrix<double>(t,t.size(),1);
    return res;
  } else {
    SWIG_Error(SWIG_TypeError, "asDMatrix: unrecognised type. Should have been caught by typemap(typecheck)");
    return false;
  }
  return true;
}

// Disallow 1D numpy arrays. Allowing them may introduce conflicts with other typemaps or overloaded methods
template <>
bool meta< CasADi::Matrix<double> >::couldbe(PyObject * p) {
  return meta< double >::couldbe(p) || ((is_array(p) && array_numdims(p)==2) && array_type(p)!=NPY_OBJECT|| PyObjectHasClassName(p,"csr_matrix") || PyObjectHasClassName(p,"DMatrix")) || meta< double >::couldbe_sequence(p);
}

%}
#endif //SWIGPYTHON


/// CasADi::Matrix<double>
#ifdef SWIGOCTAVE
%inline %{
template<> char meta< CasADi::Matrix<double> >::expected_message[] = "Expecting numpy.array2D, numpy.matrix, csr_matrix, DMatrix";

template <>
int meta< CasADi::Matrix<double> >::as(const octave_value& p,CasADi::Matrix<double> &m) {
  NATIVERETURN(CasADi::Matrix<double>,m)
  if(p.is_real_matrix()){
    Matrix mat = p.matrix_value();
    m = CasADi::DMatrix(mat.rows(),mat.cols(),0);
    for(int i=0; i<mat.rows(); ++i){
      for(int j=0; j<mat.cols(); ++j){
        m(i,j) = mat(i,j);
      }
    }
    return true;
  }
  if (p.is_real_scalar()) {
    m = CasADi::DMatrix(1,1,p.double_value());
    return true;
  } 
  return false;
}

// Disallow 1D numpy arrays. Allowing them may introduce conflicts with other typemaps or overloaded methods
template <>
bool meta< CasADi::Matrix<double> >::couldbe(const octave_value& p) { return meta< CasADi::Matrix<double> >::isa(p) || p.is_real_matrix() || p.is_real_scalar();}

%}
#endif //SWIGOCTAVE

/// CasADi::Slice
%inline %{
template<> char meta< CasADi::Slice >::expected_message[] = "Expecting Slice or number";
%}

#ifdef SWIGPYTHON
%inline %{
template <>
int meta< CasADi::Slice >::as(PyObject * p,CasADi::Slice &m) {
  NATIVERETURN(CasADi::Slice,m)

  if (PyInt_Check(p)) {
    m.start_ = PyInt_AsLong(p);
    m.stop_ = m.start_+1;
    return true;
  } else if (PySlice_Check(p)) {
    PySliceObject *r = (PySliceObject*)(p);
    if(r->start!=Py_None) m.start_ = PyInt_AsLong(r->start);
    if(r->stop !=Py_None) m.stop_  = PyInt_AsLong(r->stop);
    if(r->step !=Py_None) m.step_  = PyInt_AsLong(r->step);
    return true;
  } else {
    return false;
  }

}

template <>
bool meta<  CasADi::Slice >::couldbe(PyObject * p) {
  return meta< CasADi::Slice >::isa(p) || PyInt_Check(p) || PySlice_Check(p);
}
%}
#endif //SWIGPYTHON


/// CasADi::Slice
#ifdef SWIGOCTAVE
%inline %{

template <>
int meta< CasADi::Slice >::as(const octave_value& p,CasADi::Slice &m) {
  if (p.is_range()) {
    Range r = p.range_value();
    m.start_ = r.base()-1;
    m.stop_ = r.limit();
    m.step_ = r.inc();
  } else if (p.is_magic_colon()) {
    m.start_ = 0;
  } else if (p.is_numeric_type()) {
    m.start_ = p.int_value()-1;
    m.stop_ = m.start_+1;
  } else {
    return false;
  }
  return true;
}

template <>
bool meta<  CasADi::Slice >::couldbe(const octave_value& p) {
  return p.is_range() || p.is_magic_colon()|| (p.is_real_scalar() && p.is_numeric_type());
}

%}
#endif // SWIGOCTAVE

/// CasADi::IndexList
%inline %{
template<> char meta< CasADi::IndexList >::expected_message[] = "Expecting Slice or number or list of ints";
%}

/// CasADi::IndexList
#ifdef SWIGPYTHON
%inline %{

template <>
int meta< CasADi::IndexList >::as(PyObject * p,CasADi::IndexList &m) {
  
  if (meta< int >::couldbe(p)) {
    m.type = CasADi::IndexList::INT;
    meta< int >::as(p,m.i);
  } else if (meta< std::vector<int> >::couldbe(p)) {
    m.type = CasADi::IndexList::IVECTOR;
    return meta< std::vector<int> >::as(p,m.iv);
  } else if (meta< CasADi::Slice>::couldbe(p)) {
    m.type = CasADi::IndexList::SLICE;
    return meta< CasADi::Slice >::as(p,m.slice);
  } else {
    return false;
  }
  return true;
}


template <>
bool meta<  CasADi::IndexList >::couldbe(PyObject * p) {
  return meta< CasADi::Slice >::couldbe(p) || meta< std::vector<int> >::couldbe(p) || meta< int >::couldbe(p);
}
%}
#endif //SWIGPYTHON

/// CasADi::IndexList
#ifdef SWIGOCTAVE
%inline %{

template <>
int meta< CasADi::IndexList >::as(const octave_value& p,CasADi::IndexList &m) {
  if ((p.is_real_scalar() && p.is_numeric_type())) {
    m.type = CasADi::IndexList::INT;
    m.i = p.int_value()-1;
  } else if (meta< std::vector<int> >::couldbe(p)) {
    m.type = CasADi::IndexList::IVECTOR;
    bool result = meta< std::vector<int> >::as(p,m.iv);
    if (!result) return false;
    for (int k=0; k < m.iv.size();k++) m.iv[k]--;
  } else if (meta< CasADi::Slice>::couldbe(p)) {
    m.type = CasADi::IndexList::SLICE;
    return meta< CasADi::Slice >::as(p,m.slice);
  } else {
    return false;
  }
  return true;
}


template <>
bool meta<  CasADi::IndexList >::couldbe(const octave_value& p) {
  return meta< CasADi::Slice >::couldbe(p) || meta< std::vector<int> >::couldbe(p) || (p.is_real_scalar() && p.is_numeric_type());
}
%}
#endif //SWIGOCTAVE

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

binopsFull(const CasADi::Matrix<CasADi::SX> & b,,CasADi::Matrix<CasADi::SX>,CasADi::Matrix<CasADi::SX>)
binopsFull(const CasADi::SX & b,CasADi::Matrix<CasADi::SX>,CasADi::Matrix<CasADi::SX>,CasADi::Matrix<CasADi::SX>)
binopsFull(const CasADi::MX & b,,CasADi::MX,CasADi::MX)
binopsFull(double b,CasADi::Matrix<double>,,CasADi::Matrix<double>)

concat_operator(CasADi::Matrix<CasADi::SX>,CasADi::Matrix<double>,CasADi::Matrix<CasADi::SX>)
concat_operator(CasADi::MX,CasADi::Matrix<double>,CasADi::MX)

}; // extend Matrix<double>
} // namespace CasADi
#endif // SWIGOCTAVE


namespace CasADi{
%extend Matrix<double> {

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
  __array_priority__ = 999
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

%pythoncode %{
  def __array_wrap__(self,out_arr,context=None):
    name = context[0].__name__
    conversion = {"multiply": "mul", "divide": "div", "subtract":"sub","power":"pow"}
    if name in conversion:
      name = conversion[name]
    if len(context[1])==2 and context[1][1] is self:
      name = 'r' + name
    if not(hasattr(self,name)):
      name = '__' + name + '__'
    fun=getattr(self, name)
    return fun(*context[1][0:-1])
%}

// The following code has some trickery to fool numpy ufunc.
// Normally, because of the presence of __array__, an ufunctor like nump.sqrt
// will unleash its activity on the output of __array__
// However, we wish DMatrix to remain a DMatrix
// So when we receive a call from a functor, we return a dummy empty array
// and return the real result during the postprocessing (__array_wrap__) of the functor.
%pythoncode %{
  def __array__(self,*args,**kwargs):
    import numpy as n
    if len(args) > 1 and isinstance(args[1],tuple) and isinstance(args[1][0],n.ufunc):
      return n.array([])
    else:
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


# define binopsT(T) \
T __pow__ (const T& b) const{ return std::pow(T(*$self),b);} \
T __rpow__(const T& b) const{ return std::pow(b,T(*$self));}

binopsT(CasADi::Matrix<CasADi::SX>)
binopsT(CasADi::MX)

CasADi::Matrix<CasADi::SX> __pow__ (const CasADi::SX& b) const{ return CasADi::Matrix<CasADi::SX>(*$self).__pow__(b);}
CasADi::Matrix<CasADi::SX> __rpow__(const CasADi::SX& b) const{ return CasADi::Matrix<CasADi::SX>(b).__pow__(CasADi::Matrix<CasADi::SX>(*$self));}

#undef binopsT
    

}; // extend Matrix<double>
} // namespace CasADi
#endif // SWIGPYTHON

  
#ifdef SWIGPYTHON
%meta_vector(CasADi::Matrix<double>)
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
