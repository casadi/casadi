%{
#include "casadi/sx/sx.hpp"
#include "casadi/sx/sx_matrix.hpp"
#include "casadi/expression_tools.hpp"
%}

%include "typemaps.i"

#ifdef WITH_NUMPY
#include <numpy/arrayobject.h>

//%typemap(in) SXMatrix
// {
//     PyErr_SetString(PyExc_TypeError,"Expected an SXMatrix");
//}

#endif // WITH_NUMPY

namespace CasADi {

#ifdef WITH_IMPLICITCONV
%implicitconv SX;
#endif WITH_IMPLICITCONV

class SX {
  public:
  SX();
  SX(const std::string& name);
  SX(double val);
  ~SX();
};

%extend SX {
std::string __str__()  { return $self->toString(); }
std::string __repr__() { return $self->toString(); }
double __float__() { return (*$self)->getValue();}
int __int__() { return (*$self)->getIntValue();}
bool isConstant() const{return (*$self)->isConstant();}
bool isInteger() const{ return (*$self)->isInteger();}
bool isSymbolic() const{ return (*$self)->isSymbolic();}
bool isBinary() const{ return (*$self)->isBinary();}
bool isZero() const{ return (*$self)->isZero();}
bool isOne() const{ return (*$self)->isOne();}
bool isMinusOne() const{ return (*$self)->isMinusOne();}
bool isNan() const{ return (*$self)->isNan();}
bool isInf() const{ return (*$self)->isInf();}
bool isMinusInf() const{ return (*$self)->isMinusInf();}
const std::string& getName() const{ return (*$self)->getName();}
int getOp() const{ return (*$self)->getOp();}
bool isEqual(const SX& scalar) const{ return (*$self)->isEqual(scalar);}

//  all binary operations with a particular right argument
#define binops(T,t) \
T __add__(t b){  return *$self + b;} \
T __radd__(t b){ return b + *$self;} \
T __sub__(t b){  return *$self - b;} \
T __rsub__(t b){ return b - *$self;} \
T __mul__(t b){  return *$self * b;} \
T __rmul__(t b){ return b * *$self;} \
T __div__(t b){  return *$self / b;} \
T __rdiv__(t b){ return b / *$self;} \
T __pow__(t b){  return std::pow(*$self,b);} \
T __rpow__(t b){ return std::pow(b,*$self);} \
T fmin(t b){     return std::fmin(*$self,b);} \
T fmax(t b){     return std::fmax(*$self,b);}

// Binary operations with all right hand sides
binops(SX, const SX&)
binops(SX, double)

// all unary operations
#define unops(T) \
T __neg__(){ return - *$self;}\
T exp(){ return std::exp(*$self);}\
T log(){ return std::log(*$self);}\
T sqrt(){ return std::sqrt(*$self);}\
T sin(){ return std::sin(*$self);}\
T cos(){ return std::cos(*$self);}\
T tan(){ return std::tan(*$self);}\
T arcsin(){ return std::asin(*$self);}\
T arccos(){ return std::acos(*$self);}\
T arctan(){ return std::atan(*$self);}\
T floor(){ return std::floor(*$self);}\
T ceil(){ return std::ceil(*$self);}\
T erf(){ return std::erf(*$self);}

unops(SX)

}
//SX if_else(const SX &cond, const SX &if_true, const SX &if_false); // ternary operator, "cond ? if_true : if_false"

} // namespace CasADi

// Template instantiations
namespace std {
%template(vector_sx) vector<CasADi::SX>;
%template(vector_vector_sx) vector< vector<CasADi::SX> >;
%template(vector_vector_vector_sx) vector< vector< vector<CasADi::SX> > >;
} // namespace std;

namespace CasADi{

#ifdef WITH_IMPLICITCONV
%implicitconv SXMatrix;
#endif WITH_IMPLICITCONV

class SXMatrix : public std::vector<SX>{
public:

/** \brief  constructors */
SXMatrix();                               // empty 0-by-0 matrix
SXMatrix(int n, int m);                   // empty n-by-m matrix
SXMatrix(int n, int m, const SX& val);    // dense n-by-m matrix filled with val

/** \brief  These constructors enable implicit type conversion */
SXMatrix(const SX &scalar);      // create a variable from a scalar
SXMatrix(double val);            // constant
/*SXMatrix(const std::vector<double>& x);*/
SXMatrix(const std::vector<SX>& x);
/*SXMatrix(const std::vector<double>& x,  int n, int m);*/
SXMatrix(const std::vector<SX>& x,  int n, int m);

/** \brief  Create a matrix of symbolic variables  */
explicit SXMatrix(const std::string& name, int n=1, int m=1);   // variable

/** \brief  destructor */
~SXMatrix();
  
void clear();
void resize(int n, int m);
void reserve(int nnz);

int numel() const;       // get the number of elements
int size1() const;       // get the first dimension
int size2() const;       // get the second dimension
bool empty() const; // is the matrix empty
bool scalar() const; // is the matrix scalar
bool vector() const; // is the matrix a vector

};



%extend SXMatrix {
std::string __str__()  { return $self->getRepresentation(); }
std::string __repr__() { return $self->getRepresentation(); }

// Matrix product (quick fix)
SXMatrix dot(const SXMatrix &y){ return prod(*$self,y); }

// Get and set elements
SX __getitem__(int i){ return $self->getElement(i);}
SX __getitem__(const std::vector<int> &I ){ if(I.size()!=2) throw CasADi::CasadiException("__getitem__: not 2D"); return $self->getElement(I[0],I[1]);}

#define SETTERS(T) \
void __setitem__(int i, T el){ $self->getElementRef(i) = el;}\
void __setitem__(const std::vector<int> &I, T el){ if(I.size()!=2) throw CasADi::CasadiException("__setitem__: not 2D"); $self->getElementRef(I[0],I[1]) = el;}
SETTERS(const SX&)
SETTERS(double)
#undef SETTERS

binops(SXMatrix, const SXMatrix&)
binops(SXMatrix, double)
binops(SXMatrix, const std::vector<SX>&)
unops(SXMatrix)

}

} // namespace CasADi




namespace CasADi{
// Create an stl vector of sx variables
std::vector<SX> create_symbolic(const std::string& name, int n);
std::vector< vector<SX> > create_symbolic(const std::string& name, int n, int m);
std::vector< std::vector< std::vector< SX> > > create_symbolic(const std::string& name, int n, int m, int p);

} // namespace CasADi



#undef unops
#undef binops

//%include "numpy.i"


//namespace CasADi{

// This simple typemap will convert a returned SXMAtrix to a numpy array of SX's

//%typemap(out) SXMatrix
//{
 //    #PyObject * 
 //    PyErr_SetString(PyExc_TypeError,"Hoorah");
 //    #$result = PySXMatrix2Numpy($1);
//}
//}



