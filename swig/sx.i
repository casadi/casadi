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
%}

namespace CasADi{
/*
Attempts to form its argument into a std::vector<SXMatrix> form
Accepts: sequence(sequence(SX))
*/
%typemap(in) const std::vector<SXMatrix> &  {
  std::vector<CasADi::SXMatrix> *p = new std::vector<CasADi::SXMatrix>(PySequence_Size($input));
  int result=VSXMatrix($input,p);
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


}
