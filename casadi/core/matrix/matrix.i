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
#ifndef CASADI_MATRIX_I
#define CASADI_MATRIX_I

%include <casadi/core/printable_object.i>
%include <casadi/core/matrix/generic_matrix.i>
%include <casadi/core/matrix/generic_expression.i>

// FIXME: Placing in printable_object.i does not work
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


#ifdef SWIGPYTHON
namespace casadi{
%extend Matrix<double> {
/// Create a 2D contiguous NP_DOUBLE numpy.ndarray

#ifdef WITH_NUMPY
PyObject* arrayView() {
  if ($self->size()!=$self->numel()) 
    throw  casadi::CasadiException("Matrix<double>::arrayview() can only construct arrayviews for dense DMatrices.");
  npy_intp dims[2];
  dims[0] = $self->size2();
  dims[1] = $self->size1();
  std::vector<double> &v = $self->data();
  PyArrayObject* temp = (PyArrayObject*) PyArray_New(&PyArray_Type, 2, dims, NPY_DOUBLE, NULL, &v[0], 0, NPY_ARRAY_CARRAY, NULL);
  PyObject* ret = PyArray_Transpose(temp,NULL);
  Py_DECREF(temp); 
  return ret;
}
#endif // WITH_NUMPY
    
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
  def toCsc_matrix(self):
    import numpy as n
    import warnings
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      from scipy.sparse import csc_matrix
    return csc_matrix( (list(self.data()),self.row(),self.colind()), shape = self.shape, dtype=n.double )

  def tocsc(self):
    return self.toCsc_matrix()

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
      raise Exception("Only a scalar can be cast to an int")
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

binopsrFull(casadi::Matrix<double>)
binopsFull(const casadi::SX & b,,casadi::SX,casadi::SX)
binopsFull(const casadi::MX & b,,casadi::MX,casadi::MX)

}; // extend Matrix<double>

%extend Matrix<int> {

  %python_array_wrappers(998.0)

  %pythoncode %{
    def toArray(self):
      import numpy as n
      r = n.zeros((self.size1(),self.size2()))
      for i in range(r.shape[0]):
        for j in range(r.shape[1]):
          r[i,j] = self.elem(i,j)
      return r
  %}
  
  %pythoncode %{
    def __float__(self):
      if self.numel()!=1:
        raise Exception("Only a scalar can be cast to a float")
      if self.size()==0:
        return 0.0
      return float(self.toScalar())
  %}

  %pythoncode %{
    def __int__(self):
      if self.numel()!=1:
        raise Exception("Only a scalar can be cast to an int")
      if self.size()==0:
        return 0
      return self.toScalar()
  %}

  binopsrFull(casadi::Matrix<int>)
  binopsFull(const casadi::SX & b,,casadi::SX,casadi::SX)
  binopsFull(const casadi::Matrix<double> & b,,casadi::Matrix<double>,casadi::Matrix<double>)
  binopsFull(const casadi::MX & b,,casadi::MX,casadi::MX)
  %pythoncode %{
    def __abs__(self):
      return int(self.__int__())
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
        return {"sparsity" : self.sparsity().__getstate__(), "data": numpy.array(self.data(),dtype=int)}
  %} 
}

%extend Matrix<double> {

  %pythoncode %{
    def __setstate__(self, state):
        sp = Sparsity.__new__(Sparsity)
        sp.__setstate__(state["sparsity"])
        self.__init__(sp,state["data"])

    def __getstate__(self):
        return {"sparsity" : self.sparsity().__getstate__(), "data": numpy.array(self.data(),dtype=float)}
  %}
  
}


} // namespace casadi
#endif // SWIGPYTHON


#endif // CASADI_MATRIX_I
