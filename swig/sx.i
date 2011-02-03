%{
#include "casadi/matrix/crs_sparsity.hpp"
#include "casadi/matrix/matrix.hpp"
#include "casadi/sx/sx.hpp"
#include "casadi/sx/sx_tools.hpp"
%}

%include "typemaps.i"
%include "casadi/matrix/crs_sparsity.hpp"
%include "casadi/matrix/matrix.hpp"
%include "casadi/sx/sx.hpp"

%extend std::vector<CasADi::SX>{
  std::string __repr__(){ return CasADi::getRepresentation(*$self); }
  std::string __str__(){ return CasADi::getDescription(*$self); }
};

%extend std::vector<CasADi::Matrix< CasADi::SX> >{
  std::string __repr__(){ return CasADi::getRepresentation(*$self); }
  std::string __str__(){ return CasADi::getDescription(*$self); }
};

%extend std::vector<std::vector< CasADi::SX> >{
  std::string __repr__(){ return CasADi::getRepresentation(*$self); }
  std::string __str__(){ return CasADi::getDescription(*$self); }
};

namespace CasADi {
  %extend Matrix<double> {
  

    
    %pythoncode %{
    def __getitem__(self,s):
      if isinstance(s,tuple) and len(s)==2 and (isinstance(s[1],slice) or isinstance(s[0],slice)):
        s = list(s)
        for k in range(2):
          if isinstance(s[k],slice):
            J = s[k].indices(self.size1())
            s[k] = range(J[0],J[1],J[2])
          elif isinstance(s[k],int):
            s[k] = [s[k]]
      return self.getitem(s)
    %}
    
    %pythoncode %{
    def __setitem__(self,s,val):
      if isinstance(s,tuple) and len(s)==2 and (isinstance(s[1],slice) or isinstance(s[0],slice)):
        s = list(s)
        for k in range(2):
          if isinstance(s[k],slice):
            J = s[k].indices(self.size1())
            s[k] = range(J[0],J[1],J[2])
          elif isinstance(s[k],int):
            s[k] = [s[k]]
      self.setitem(s,val)
    %}

  };

%extend Matrix<SX>{
    // The constructor has to be added since SX::operator Matrix<SX does not work
    // Matrix<SX>(const SX&){ *$self
    
    
    %pythoncode %{
      def __lt__(self,other):
        return _casadi.__lt__(self,other)
      def __le__(self,other):
        return _casadi.__le__(self,other)
      def __eq__(self,other):
        return _casadi.__eq__(self,other)
      def __ne__(self,other):
        return _casadi.__ne__(self,other)
      def __gt__(self,other):
        return _casadi.__gt__(self,other)
      def __ge__(self,other):
        return _casadi.__ge__(self,other)
    %}
    
    // These methods must be added since the implicit type cast does not work
    Matrix<SX> __pow__ (double b) const{ return $self->__pow__(CasADi::Matrix<CasADi::SX>(b));}
    Matrix<SX> __rpow__(double b) const{ return CasADi::Matrix<CasADi::SX>(b).__pow__(*$self);}
    Matrix<SX> __add__ (double b) const{ return *$self + CasADi::Matrix<CasADi::SX>(b);}
    Matrix<SX> __radd__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) + *$self;}
    Matrix<SX> __sub__ (double b) const{ return *$self - CasADi::Matrix<CasADi::SX>(b);}
    Matrix<SX> __rsub__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) - *$self;}
    Matrix<SX> __mul__ (double b) const{ return *$self * CasADi::Matrix<CasADi::SX>(b);}
    Matrix<SX> __rmul__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) * *$self;}
    Matrix<SX> __div__ (double b) const{ return *$self / CasADi::Matrix<CasADi::SX>(b);}
    Matrix<SX> __rdiv__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) / *$self;}

    %pythoncode %{
        @property
        def shape(self):
            return (self.size1(),self.size2())
    %}

    %pythoncode %{
    def toArray(self):
      import numpy as n
      r = n.array((),dtype=object)
      r.resize(self.size1(),self.size2())
      for i in range(self.size1()):  # loop over rows
        for el in range(self.rowind(i),self.rowind(i+1)): # loop over the non-zero elements
          j=self.col(el)  # column
          r[i,j] = self[el] # add the non-zero element

      return r
    %}
    
    
  %pythoncode %{
  def __array_wrap__(self,*args,**kwargs):
    fun=getattr(self, args[1][0].__name__)
    return fun()
  %}

  %pythoncode %{
    def __array__(self,*args,**kwargs):
      import numpy as n
      if len(args) > 1 and isinstance(args[1],tuple) and isinstance(args[1][0],n.ufunc):
        return n.array([])
      else:
        return self.toArray()
  %}
        
    %pythoncode %{
    def toMatrix(self):
      import numpy as n
      return n.matrix(self.toArray())
    %}
    
    %pythoncode %{
    def __getitem__(self,s):
      if isinstance(s,tuple) and len(s)==2 and (isinstance(s[1],slice) or isinstance(s[0],slice)):
        s = list(s)
        for k in range(2):
          if isinstance(s[k],slice):
            J = s[k].indices(self.size1())
            s[k] = range(J[0],J[1],J[2])
          elif isinstance(s[k],int):
            s[k] = [s[k]]
      return self.getitem(s)
    %}

    %pythoncode %{
    def __setitem__(self,s,val):
      if isinstance(s,tuple) and len(s)==2 and (isinstance(s[1],slice) or isinstance(s[0],slice)):
        s = list(s)
        for k in range(2):
          if isinstance(s[k],slice):
            J = s[k].indices(self.size1())
            s[k] = range(J[0],J[1],J[2])
          elif isinstance(s[k],int):
            s[k] = [s[k]]
      self.setitem(s,val)
    %}
};
} // namespace CasADi


#ifdef WITH_NUMPY
#include <numpy/arrayobject.h>
#endif // WITH_NUMPY

// Template instantiations
%template(vector_PyObject)    std::vector<PyObject*>;


%inline %{

bool istype(PyObject *p, swig_type_info *type) {
  int res = SWIG_ConvertPtr(p, 0, type, 0);
  return SWIG_CheckState(res);
}

bool isSXMatrix(PyObject * p) {
  return istype(p,SWIGTYPE_p_CasADi__MatrixT_CasADi__SX_t);
}

bool isSX(PyObject * p) {
  return istype(p,SWIGTYPE_p_CasADi__SX);
}

int getSXMatrix(PyObject * p, CasADi::SXMatrix * & m) {
  void *pd = 0 ;
  int res = SWIG_ConvertPtr(p, &pd,SWIGTYPE_p_CasADi__MatrixT_CasADi__SX_t, 0 );
  if (!SWIG_IsOK(res)) {
    return false;
  }
  m = reinterpret_cast< CasADi::SXMatrix* >(pd);
  return true;
}

int getSX(PyObject * p, CasADi::SX * & m) {
  void *pd = 0 ;
  int res = SWIG_ConvertPtr(p, &pd,SWIGTYPE_p_CasADi__SX, 0 );
  if (!SWIG_IsOK(res)) {
    return false;
  }
  m = reinterpret_cast< CasADi::SX * >(pd);
  return true;
}

bool VSXMatrix(PyObject *p, std::vector<CasADi::SXMatrix> * v) {
  if (PySequence_Check(p)) {                                                  // Look for a sequence
     PyObject *ite = PyObject_GetIter(p);
     PyObject *pe;
     int i=-1;
     while (pe = PyIter_Next(ite)) {                                           // Iterate over the sequence
        i++;
        if (isSXMatrix(pe)) {
          void *mp = 0;
          SWIG_ConvertPtr(pe, &mp, SWIGTYPE_p_CasADi__MatrixT_CasADi__SX_t, 0);
          CasADi::SXMatrix *m=reinterpret_cast< CasADi::SXMatrix * >(mp);
          v->operator[](i) = *m;
        } else if (PySequence_Check(pe)) {                                        // Look for a sequence inside the sequence
            PyObject *itee = PyObject_GetIter(pe);
            PyObject *pee;
            std::vector<CasADi::SX> sxv(PySequence_Size(pe));
            int j=-1;
            while (pee = PyIter_Next(itee)) {                                // Iterate over the sequence inside the sequence
              j++;
              if (!isSX(pee))  {                                                       // Make sure we have SX elements
                Py_DECREF(pee);Py_DECREF(itee);Py_DECREF(ite);return 0;  
              } else {
                void *sp = 0;
                SWIG_ConvertPtr(pee, &sp, SWIGTYPE_p_CasADi__SX, 0);
                sxv[j] = *(reinterpret_cast< CasADi::SX * >(sp));
              }
              Py_DECREF(pee);
            }
            v->operator[](i) = CasADi::SXMatrix(sxv);
            Py_DECREF(itee);
        }else {
          Py_DECREF(pe);Py_DECREF(ite);return 0;
        }
        Py_DECREF(pe);
    }
    Py_DECREF(ite);
    return 1;
  } else {
    return 0;
  }
}

/**
*  \param   PyObject* p             input: python object
*                                     numpy.ndarray(SX) or SXMatrix
*
*  \param   CasADi::SXMatrix * m    output: pointer to an SXMatrix
*  \param   bool& freearg           output: boolean flag indicating wheter m pointer must be deleted in typemap(freearg)
*
*  \return bool                     indicates succes
*/
bool typemapSXMatrixHelper(PyObject* p, CasADi::Matrix< CasADi::SX > * & m,  bool& freearg) {
  if (is_array(p)) { // Numpy arrays will be cast to dense Matrix<double>
    std::cout << "array" << std::endl;
		if (array_type(p)!=NPY_OBJECT) {
			SWIG_Error(SWIG_TypeError, "asSXMatrix: numpy.ndarray must be of dtype object");
			return false;
		}
		if (array_numdims(p)>2 || array_numdims(p)<1) {
			SWIG_Error(SWIG_TypeError, "asSXMatrix: Number of dimensions must be 1 or 2.");
			return false;
		}
		int nrows = array_size(p,0); // 1D array is cast into column vector
		int ncols  = 1;
		if (array_numdims(p)==2)
			ncols=array_size(p,1); 
		int size=nrows*ncols; // number of elements in the dense matrix
		std::vector<CasADi::SX> v(size);
    PyArrayIterObject* it = (PyArrayIterObject*)PyArray_IterNew(p);
    PyObject *pe;
    int i=-1;
		while (it->index < it->size) { 
		  i++;
		  pe = *((PyObject**) PyArray_ITER_DATA(it));
		  void *pd = 0;
		  int res = SWIG_ConvertPtr(pe, &pd,SWIGTYPE_p_CasADi__SX, 0 );
		  if (!SWIG_IsOK(res)) {
		    Py_DECREF(it);
			  SWIG_Error(SWIG_TypeError, "asSXMatrix: SXMatrix cast problem");
			  return false;
		  }
		  v[i] = *(reinterpret_cast< CasADi::SX * >(pd));
		  PyArray_ITER_NEXT(it);
		}
    Py_DECREF(it);

		m = new CasADi::Matrix< CasADi::SX >(v, nrows, ncols);
		freearg = true; // Memory will be freed in typemap(freearg)
	} else if (isSXMatrix(p)) { // SXMatrix object get passed on as-is.
    int result = getSXMatrix(p,m);
		if (!result) {
			SWIG_Error(SWIG_TypeError, "asSXMatrix: SXMatrix cast problem");
			return false;
		}
	} else if (isSX(p)) {
    CasADi::SX * sx;
    int result = getSX(p,sx);
		if (!result) {
			SWIG_Error(SWIG_TypeError, "asSXMatrix: SX cast problem");
			return false;
		}
    m = new CasADi::Matrix< CasADi::SX >(*sx);
    freearg=true; / Memory will be freed in typemap(freearg)
  } else {
    SWIG_Error(SWIG_TypeError, "asSXMatrix: unrecognised type. Should have been caught by typemap(typecheck)");
    return false;
  }
	return true;
}

%}

namespace CasADi{
/*
Attempts to form its argument into a std::vector<SXMatrix> form
Accepts: sequence(sequence(SX)), sequence(SXMatrix), sequence(numpy.ndarray(SX))
*/
%typemap(in) const std::vector<SXMatrix> &  {
  std::vector<CasADi::SXMatrix> *p = new std::vector<CasADi::SXMatrix>(PySequence_Size($input));
  bool result=VSXMatrix($input,p);
  if (!result) {
    delete p;
    SWIG_exception_fail(SWIG_TypeError,"Expecting sequence(sequence(SX)) or sequence(SXMatrix)");
  }
  $1 = p;
}

%typemap(freearg) const std::vector<SXMatrix> & {
    if ($1) delete $1;
}

%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) const std::vector<SXMatrix> & {
    if (PySequence_Check($input)) {
      $1 = 1;
    } else {
      $1=0;
    }
}

%typemap(in) const Matrix<SX> & (bool freearg = false) {
  CasADi::Matrix<CasADi::SX> * m = 0;
  bool result=typemapSXMatrixHelper($input, m, freearg);
  if (!result) {
    if (freearg && m) delete m;
    SWIG_exception_fail(SWIG_TypeError,"Expecting SXMatrix or nump.ndarray(SX)");
  }
  $1 = m;
}


%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) const Matrix<SX> & {
  PyObject* p = $input;
  if ((is_array(p) && array_numdims(p) < 3) && array_type(p)==NPY_OBJECT || isSXMatrix(p) || isSX(p)) {
    $1=1;
  } else {
    $1=0;
  }
}


%typemap(freearg) const Matrix<SX> & {
   if (freearg$argnum) {
    if ($1) delete $1;
   }
}


}
