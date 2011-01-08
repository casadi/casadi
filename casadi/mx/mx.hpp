/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef MX_HPP
#define MX_HPP

#include "../shared_object.hpp"
#include "../sx/sx.hpp"
#include "slicer.hpp"
#include <vector>

namespace CasADi{

/** \brief  Forward declaration */
class MXNode;

/** \brief The main matrix class.
  Operations on MX are lazy on the matrix level.

  Index flatten happens as follows: (i,j) -> k = j+i*size2()

  \author Joel Andersson 
  \date 2010
*/
class MX : public SharedObject{
public:
 
/** \brief  Default constructor */
  MX();

/** \brief  Construct a symbolic matrix (matrix variable)
Internally represented by SymbolicMatrix
*/
  explicit MX(const std::string& name, int n=1, int m=1);

/** \brief  Create constant */
  MX(double x);
  MX(const std::vector<double> &x);
  MX(const std::vector<double> &x, int n, int m=1, char order='R');
#ifdef WITH_CPP0X
  MX(std::initializer_list<double> list);
#endif

/** \brief  Destructor */
  ~MX();

#ifndef SWIG
/** \brief  Immutable matrix style access
   Get an element of the matrix (without changing this object)
 */  
 const MX operator()(int i, int j) const;
/** \brief  Immutable vector style access
   Get an element of the matrix (without changing this object)
 */  
 const MX operator[](int k) const;          

/** \brief  Element of the matrix which is allowed to change the object (vector style index)
  \author Joel Andersson 
  \date 2010
  Allows to get elements from an MX and changing them.
  Elements are scalars, use other techniques for having sliced access.
  \see MatrixElement	
*/
  class Element : public PrintableObject{
    public:
      /** \brief Vector style access constructor */
      Element(MX& mx, int k);
      
      /** \brief  Print */
      virtual void print(std::ostream &stream=std::cout) const;
      
      /** \brief  Automatic type conversion (? =  A[i], sqrt(A[i]) etc.) */
      operator MX() const;
      //@{
      /** \brief  Objects that modify a part of the parent obejct (A[i] = ?, A[i] += ?, etc.) */
      MX& operator=(const MX &y);
      MX& operator+=(const MX &y);
      MX& operator-=(const MX &y);
      MX& operator*=(const MX &y);
      MX& operator/=(const MX &y);
      //@}
    
    protected:
      MX& mx;
      int k;
  };

/** \brief  Mutable matrix style access
   Get an element of the matrix (possibly changing this object)
 */  
  Element operator()(int i, int j);
/** \brief  Mutable vector style access
   Get an element of the matrix (possibly changing this object)
 */  
  Element operator[](int k);
#endif // SWIG

  /** \brief  Get the number of (structural) non-zero elements */
  int size_new() const;

  /** \brief  Get the number of elements */
  int numel() const;

  /** \brief get the first dimension
  For an n-by-m matrix, returns n
  */
  int size1() const;
/** \brief get the second dimension
For an n-by-m matrix, returns m
  */
  int size2() const;

//@{
/** \brief  Operators that changes the object */
  MX& operator+=(const MX &y);
  MX& operator-=(const MX &y);
  MX& operator*=(const MX &y);
  MX& operator/=(const MX &y);
//@}

#ifndef SWIG
//@{
/** \brief  Operators */
  friend MX operator+(const MX &x, const MX &y);
  friend MX operator-(const MX &x, const MX &y);
  friend MX operator*(const MX &x, const MX &y);
  friend MX operator/(const MX &x, const MX &y);
  MX operator-() const;
//@}
#endif // SWIG

  /** \brief  Check if the matrix expression is empty */
  bool isEmpty() const;
  
/** \brief  Initialize the tree */
/** \brief    void init(); */

#ifndef SWIG
/** \brief  Lowlevel. Get a pointer to the node
A regular user should never use this.
*/
  MXNode* get();
/** \brief  Lowlevel. Get a pointer to the node
A regular user should never use this.
*/
  const MXNode* get() const;
#endif // SWIG

//@{
/** \brief  Quick access a member of the node */
  MXNode* operator->();
  const MXNode* operator->() const;
//@}

//@{
  /** \brief  Create nodes by their ID */
  static MX binary(int op, const MX &x, const MX &y);
  static MX unary(int op, const MX &x);
  static MX scalar_matrix(int op, const MX &x, const MX &y);
  static MX matrix_scalar(int op, const MX &x, const MX &y);
  static MX matrix_matrix(int op, const MX &x, const MX &y);
//@}

  /** \brief  Matrix of all zeros */  
  static MX zeros(int nrow, int ncol);
  
  /** \brief  Matrix of all ones */  
  static MX ones(int nrow, int ncol);
  
  //! \brief delayed setting or getting an element
  MX getElement(int k) const;
  MX& setElement(const MX& el, int k);

#ifndef SWIG
  //! \brief Get a slice of an MX
  MX slice(Slicer i, Slicer j) const;
#endif // SWIG

  //! \brief Get a row slice of an MX
  MX getRow(int i) const;

  //! \brief Get a column slice of an MX
  MX getColumn(int j) const;
  
};

} // namespace CasADi

#ifndef SWIG  
namespace std{
//@{
/** \brief  Functions with c equivalents: The implementation and syntax mirrors the standard c functions in math.h */
#define MX CasADi::MX
MX sqrt(const MX &x);
MX sin(const MX &x);
MX cos(const MX &x);
MX tan(const MX &x);
MX atan(const MX &x);
MX asin(const MX &x);
MX acos(const MX &x);
MX exp(const MX &x);
MX log(const MX &x);
MX pow(const MX &x, const MX &n);
MX abs(const MX &x);
MX fabs(const MX &x); // same as abs
MX floor(const MX &x);
MX ceil(const MX &x);
MX erf(const MX &x);
MX fmin(const MX &a, const MX &b);
MX fmax(const MX &a, const MX &b);
#undef MX
//@}
} // namespace std

#else // SWIG

namespace CasADi {

%extend MX {
std::string __repr__() { return $self->getRepresentation(); }


// Get or set an element
MX __getitem__(int k){ return $self->getElement(k);}
MX __getitem__(const std::vector<int> &I ){if(I.size()!=2) throw CasADi::CasadiException("__getitem__: not 2D"); return $self->operator()(I[0],I[1]);}
//MX __getitem__(PyObject* list){
//              if(PyList_Size(list)!=2) throw CasADi::CasadiException("__getitem__: not 2D");
//              return $self->slice(CasADi::Slicer(PyList_GetItem(list, 0)),CasADi::Slicer(PyList_GetItem(list, 1))
//              );
//}

MX __getitem__(PyObject* list){
    if (!PyTuple_Check(list) && ($self->size1()==1 || $self->size2()==1)) {
      CasADi::Slicer *i;
                  bool succes=true;

                  if (PyInt_Check(list)) {
                          i=new CasADi::Slicer(PyInt_AsLong(list));
                  } else if (PyList_Check(list)) {
                          std::vector<int> arg;
                          for (int l=0;l<PyList_Size(list);l++) {
                                  arg.push_back(PyInt_AsLong(PyList_GetItem(list,l)));
                          }
                          i=new CasADi::Slicer(arg);
                  } else if (PyObject_TypeCheck(list,&PySlice_Type)) {
                          i=new CasADi::Slicer(PySliceObjectToSlicerPrimitiveFromTo((PySliceObject*)list));
                  } else {
                          succes=false;
                  }
                  if (succes) {
                    if ($self->size1()==1)
                            return $self->slice(CasADi::Slicer(0),*i);
                    if ($self->size2()==1)
                            return $self->slice(*i,CasADi::Slicer(0));
                  } else {
                          if (succes) delete i;
                          throw CasADi::CasadiException("__getitem__: wrong arguments");
                  }
                        if (succes) delete i;
    }
                if (!PyTuple_Check(list))  throw CasADi::CasadiException("__getitem__: expecting tuple");
                if(PyTuple_Size(list)!=2) throw CasADi::CasadiException("__getitem__: not 2D");
                CasADi::Slicer *i[2];
                bool delme[2];
                bool succes=true;

                for (int k=0;k<2;k++) {
                        delme[k]=true;
                        if (PyInt_Check(PyTuple_GetItem(list, k))) {
                                i[k]=new CasADi::Slicer(PyInt_AsLong(PyTuple_GetItem(list, k)));
                        } else if (PyList_Check(PyTuple_GetItem(list, k))) {
                                std::vector<int> arg;
                                for (int l=0;l<PyList_Size(PyTuple_GetItem(list, k));l++) {
                                        arg.push_back(PyInt_AsLong(PyList_GetItem(PyTuple_GetItem(list, k),l)));
                                }
                                i[k]=new CasADi::Slicer(arg);
                        } else if (PyObject_TypeCheck(PyTuple_GetItem(list, k),&PySlice_Type)) {
                                i[k]=new CasADi::Slicer(PySliceObjectToSlicerPrimitiveFromTo((PySliceObject*)PyTuple_GetItem(list, k)));
                        } else {
                                succes=false;
                                delme[k]=false;
                        }
                }
                if (succes) {
                        return $self->slice(*i[0],*i[1]);
                } else {
                        if (delme[0]) delete i[0];
                        if (delme[1]) delete i[1];
                        throw CasADi::CasadiException("__getitem__: wrong arguments");
                }
                if (delme[0]) delete i[0];
                if (delme[1]) delete i[1];
}

//MX __setitem__(int k, const MX& el){ return $self->setElement(el,k);}
MX __setitem__(int k, const MX& el){ $self->getElement(k) = el; return *$self;}
MX __setitem__(const std::vector<int> &I, const MX&  el){ if(I.size()!=2) throw CasADi::CasadiException("__setitem__: not 2D"); $self->operator()(I[0],I[1]) = el; return *$self;}

// all binary operations with a particular right argument
#define binops(t) \
MX __add__(t b){    return *$self + b;} \
MX __radd__(t b){   return b + *$self;} \
MX __sub__(t b){    return *$self - b;} \
MX __rsub__(t b){   return b - *$self;} \
MX __mul__(t b){    return *$self * b;} \
MX __rmul__(t b){   return b * *$self;} \
MX __div__(t b){    return *$self / b;} \
MX __rdiv__(t b){   return b / *$self;} \
MX __pow__(t b){    return std::pow(*$self,b);} \
MX __rpow__(t b){   return std::pow(b,*$self);} \
MX fmin(t b){       return std::fmin(*$self,b);} \
MX fmax(t b){       return std::fmax(*$self,b);} \
MX prod(t b){       return prod(*$self,b);} \
MX inner_prod(t b){ return inner_prod(*$self,b);} \
MX outer_prod(t b){ return outer_prod(*$self,b);} \


// Binary operations with all right hand sides
binops(const MX&)
binops(double)
#undef binops

// all unary operations
MX __neg__(){ return - *$self;}
MX exp(){ return std::exp(*$self);}
MX log(){ return std::log(*$self);}
MX sqrt(){ return std::sqrt(*$self);}
MX sin(){ return std::sin(*$self);}
MX cos(){ return std::cos(*$self);}
MX tan(){ return std::tan(*$self);}
MX arcsin(){ return std::asin(*$self);}
MX arccos(){ return std::acos(*$self);}
MX arctan(){ return std::atan(*$self);}
MX floor(){ return std::floor(*$self);}
MX ceil(){ return std::ceil(*$self);}
MX erf(){ return std::erf(*$self);}
}

} // namespace CasADi

// Template instantiations
namespace std {
%template(vector_mx) vector<CasADi::MX>;
} // namespace std;

#endif // SWIG

#endif // MX_HPP
