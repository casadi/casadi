%{
#include "casadi/mx/mx.hpp"
#include <ostream>
#include <iostream>
%}

%include "slicer.i"

namespace CasADi {

#ifdef WITH_IMPLICITCONV
%implicitconv SX;
#endif WITH_IMPLICITCONV

class MatrixSize{
  public:
  MatrixSize();
  MatrixSize(int nrow_, int ncol_);
};

%typemap(in) const MatrixSize & {
       $1 = new CasADi::MatrixSize(PyTupleToMatrixSize($input));
}

%typemap(freearg) const MatrixSize & {
    if ($1) {
        delete $1;
    }
}



};


%inline %{
CasADi::MatrixSize PyTupleToMatrixSize(PyObject* tup) {
		  if (!PyTuple_Check(tup))  throw CasADi::CasadiException("__getitem__: expecting tuple");
      if(PyTuple_Size(tup)!=2) throw CasADi::CasadiException("__getitem__: not 2D");
      return CasADi::MatrixSize(PyInt_AsLong(PyTuple_GetItem(tup,0)),PyInt_AsLong(PyTuple_GetItem(tup,1)));
}
%}

namespace CasADi {

class MX : public SharedObject{
  public:
  MX();
  explicit MX(const std::string& name, int n=1, int m=1);
  MX(double x);
  MX(const std::vector<double> &x);
  MX(const std::vector<double> &x, int n, int m=1, char order='R');
  ~MX();

  int size1() const;       // get the first dimension
  int size2() const;       // get the second dimension
  int numel() const;       // get the total number of elements
  
  //Create nodes by their ID 
  static MX binary(int op, const MX &x, const MX &y);
  static MX unary(int op, const MX &x);
  static MX scalar_matrix(int op, const MX &x, const MX &y);
  static MX matrix_scalar(int op, const MX &x, const MX &y);
  static MX matrix_matrix(int op, const MX &x, const MX &y);

  //Matrix of all zeros   
  static MX zeros(int nrow, int ncol);
  
  //Matrix of all ones   
  static MX ones(int nrow, int ncol);
  
  // delayed setting or getting an element
  MX getElement(int k) const;
  MX& setElement(const MX& el, int k);

  //! \brief Get a row slice of an MX
  MX getRow(int i) const;

  //! \brief Get a column slice of an MX
  MX getColumn(int j) const;

};

%extend MX {
std::string __repr__() { return $self->getRepresentation(); }


// Get or set an element
MX __getitem__(int k){ return $self->getElement(k);}
MX __getitem__(const std::vector<int> &I ){if(I.size()!=2) throw CasADi::CasadiException("__getitem__: not 2D"); return $self->operator()(I[0],I[1]);}
//MX __getitem__(PyObject* list){
//		if(PyList_Size(list)!=2) throw CasADi::CasadiException("__getitem__: not 2D");
//		return $self->slice(CasADi::Slicer(PyList_GetItem(list, 0)),CasADi::Slicer(PyList_GetItem(list, 1))
//		);
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
MX norm_2(){  return norm_2(*$self);}
MX norm_1(){  return norm_1(*$self);}
MX norm_inf(){  return norm_inf(*$self);}
MX trans(){  return trans(*$self);}
}

} // namespace CasADi

// Template instantiations
namespace std {
%template(vector_mx) vector<CasADi::MX>;
} // namespace std;

namespace CasADi{
// concatenate
MX vertcat(const std::vector<MX>& comp);
MX horzcat(const std::vector<MX>& comp);
MX vertcat(const MX& a, const MX& b);
MX horzcat(const MX& a, const MX& b);

// Functions
MX norm_2(const MX &x);
MX norm_1(const MX &x);
MX norm_inf(const MX &x);
MX trans(const MX &x); // transpose
MX flatten(const MX &x); // flatten
MX reshape(const MX &x, const MatrixSize &s); // reshape
MX prod(const MX &x, const MX &y); // matrix product
MX inner_prod(const MX &x, const MX &y); // trans(x)*y with x and y vectors
MX outer_prod(const MX &x, const MX &y); // x*trans(y) with x and y vectors
MX if_else(const MX &cond, const MX &if_true, const MX &if_false); // ternary operator, "cond ? if_true : if_false"


} // namespace CasADi


