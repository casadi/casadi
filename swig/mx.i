%{
#include "casadi/mx/mx.hpp"
%}

namespace CasADi {

#ifdef WITH_IMPLICITCONV
%implicitconv SX;
#endif WITH_IMPLICITCONV

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

};

%extend MX {
std::string __repr__() { return $self->getRepresentation(); }


// Get or set an element
MX __getitem__(int k){ return $self->getElement(k);}
//MX __setitem__(int k, const MX& el){ return $self->setElement(el,k);}
MX __setitem__(int k, const MX& el){ $self->getElement(k) = el; return *$self;}

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
MX prod(const MX &x, const MX &y); // matrix product
MX inner_prod(const MX &x, const MX &y); // trans(x)*y with x and y vectors
MX outer_prod(const MX &x, const MX &y); // x*trans(y) with x and y vectors
MX if_else(const MX &cond, const MX &if_true, const MX &if_false); // ternary operator, "cond ? if_true : if_false"
} // namespace CasADi



