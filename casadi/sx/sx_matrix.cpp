#include "sx_matrix.hpp"

using namespace std;

namespace CasADi{


SXMatrix::SXMatrix(const std::vector<SX>& x) : Matrix<SX>(x){
}

SXMatrix::SXMatrix(const std::vector<double>& x) : Matrix<SX>(x){
}

SXMatrix::SXMatrix(const std::vector<SX>& x,  int n, int m) : Matrix<SX>(x,n,m){
}

SXMatrix::SXMatrix(const std::vector<double>& x,  int n, int m) : Matrix<SX>(x,n,m){
}

SXMatrix::SXMatrix() : Matrix<SX>(){
}

SXMatrix::SXMatrix(double val) : Matrix<SX>(val){
}

SXMatrix::SXMatrix(int n, int m) : Matrix<SX>(n,m){
}

SXMatrix::SXMatrix(int n, int m, const SX& val) : Matrix<SX>(n,m,val){
}

SXMatrix::SXMatrix(const string& name, int n, int m){
  nrow = ncol = 0;
  if(n*m == 0) return;
  resize(n,m);
  reserve(n*m);
  if(scalar())
    getElementRef() = SX(name);
  else
      for(int i=0; i<n; ++i)
	for(int j=0; j<m; ++j){
	  stringstream name_;
 	  name_ << name << "_" << i;
	  if(m>1) name_ << "_" << j;
	  getElementRef(i,j) = SX(name_.str());
	}
}

SXMatrix::SXMatrix(const SX& scalar){
  nrow = ncol = 0;
  resize(1,1);
  getElementRef() = scalar;
}

SXMatrix::SXMatrix(const SXMatrix::Element &el){
  nrow = ncol = 0;
  resize(1,1);
  getElementRef() = el.mat.getElement(el.i,el.j);  
}

  // destructor
SXMatrix::~SXMatrix(){
}

SX& SXMatrix::operator[](int k){
  return at(k); //  return Element(*this,k%size1(), k/size1());
}

const SX& SXMatrix::operator[](int k) const{
  return at(k); //  return getElement(k%size1(), k/size1());
}

const std::vector<SX> SXMatrix::operator[](const std::vector<int>& ind) const{
  std::vector<SX> ret;
  for(std::vector<int>::const_iterator it=ind.begin(); it!=ind.end(); ++it)
    ret.push_back(at(*it));
  return ret;
}

const SX SXMatrix::operator()(int i, int j) const{
  return getElement(i,j);
}


SXMatrix::Element SXMatrix::operator()(int i, int j){
  return Element(*this,i,j);
}

SXMatrix& operator+=(SXMatrix &ex, const SXMatrix &expr){
  return ex = ex + expr;
}

SXMatrix operator+(const SXMatrix &x, const SXMatrix &y){
  return binary(ADD_NODE,x,y);
}

SXMatrix& operator-=(SXMatrix &ex, const SXMatrix &expr){
 return ex = ex - expr;
}

SXMatrix operator-(SXMatrix &ex){
  return unary(NEG_NODE,ex);
}

SXMatrix operator-(const SXMatrix &x, const SXMatrix &y){ 
  return binary(SUB_NODE,x,y);
}

SXMatrix& operator*=(SXMatrix &ex, const SXMatrix &expr){
 return ex = ex * expr;
}

SXMatrix operator*(const SXMatrix &x, const SXMatrix &y){
  return binary(MUL_NODE,x,y);
}


SXMatrix& operator/=(SXMatrix &ex, const SXMatrix &expr){
  return ex = ex / expr;
}

SXMatrix operator/(const SXMatrix &x, const SXMatrix &y){
  return binary(DIV_NODE,x,y);
}

void SXMatrix::reserve(int nnz){
  std::vector<SX>::reserve(nnz);
  col.reserve(nnz);
}

void SXMatrix::clear(){
  *this = SXMatrix();
}

void SXMatrix::fill(const SX& val){
  if(val->isZero())    makeEmpty(size1(),size2());
  else                 makeDense(size1(),size2(),val);
}

void SXMatrix::resize(int n_, int m_){
  if(n_==size1() && m_== size2()) // if the dimensions remain the same
    return; // do nothing
  else if(n_ < size1() || m_ < size2()){ // if we need to remove some elements
    SXMatrix newexpr(n_,m_);
    const SXMatrix &old = *this;
    for(int i=0; i<n_ && i<size1(); ++i)
      for(int j=0; j<m_ && j<size2(); ++j){
        SX temp = old(i,j);
        if(!temp->isZero())
          newexpr.getElementRef(i,j) = temp;
      }

    *this = newexpr;
  } else { // make the object larger (CHEAP!)
    nrow = n_; ncol = m_;
    rowind.resize(size1()+1,size());
  }
}

// const MatrixSize& SXMatrix::size() const{
//   return sz;
// }




SXMatrix::Element::Element(SXMatrix& mat_, int i_, int j_) : mat(mat_), i(i_), j(j_){ 
}

SXNode* const SXMatrix::Element::get() const{
  return mat.getElement(i,j).get();
}

const SXNode* SXMatrix::Element::operator->() const{
  return mat.getElement(i,j).get();  
}

SXNode* SXMatrix::Element::operator->(){
  return mat.getElement(i,j).get();  
}

SX& SXMatrix::Element::operator=(const SXMatrix::Element &y){
  return mat.getElementRef(i,j) = y.mat.getElement(y.i,y.j);
}

SX& SXMatrix::Element::operator=(const SXMatrix &y){
  if(!y.scalar()) throw CasadiException("operator=: argument not scalar");
  return mat.getElementRef(i,j) = y(0);
}

SXMatrix binary(int op, const SXMatrix &x, const SXMatrix &y){
  SXMatrix r;
  if(x.scalar())
    if(y.scalar())
      return sfcn[op](x(0),y(0));
    else
      return scalar_matrix(op,x(0),y);
  else if(y.scalar())
    return matrix_scalar(op,x,y(0));
  else
    return matrix_matrix(op,x,y);
}

SXMatrix unary(int op, const SXMatrix &x){
  if(x.scalar()){
    return sfcn[op](x(0),SX::nan);
  } else {
    SXMatrix r(x.size1(),x.size2());
    SX y; // dummy argument
    r.fill(sfcn[op](0,y));
    for(int i=0; i<r.size1(); ++i){ // loop over rows
      for(int el=x.rowind[i]; el<x.rowind[i+1]; ++el){
        int j = x.col[el];
        r(i,j) = sfcn[op](x[el],y);
      }
    }  
    return r;
  }
}

SXMatrix scalar_matrix(int op, const SX &x, const SXMatrix &y){
  SXMatrix r(y.size1(),y.size2());
  r.fill(sfcn[op](x,0));
  for(int i=0; i<r.size1(); ++i){ // loop over rows
    for(int el=y.rowind[i]; el<y.rowind[i+1]; ++el){
      int j = y.col[el];
      r(i,j) = sfcn[op](x,y[el]);
    }
  }
  return r;
}

SXMatrix matrix_scalar(int op, const SXMatrix &x, const SX &y){
  SXMatrix r(x.size1(),x.size2());
  r.fill(sfcn[op](0,y));
  for(int i=0; i<r.size1(); ++i){ // loop over rows
    for(int el=x.rowind[i]; el<x.rowind[i+1]; ++el){
      int j = x.col[el];
      r(i,j) = sfcn[op](x[el],y);
    }
  }
  return r;
}

SXMatrix matrix_matrix(int op, const SXMatrix &x, const SXMatrix &y){
if(x.size1() != y.size1() || x.size2() != y.size2()) throw CasadiException("matrix_matrix: dimension mismatch");
  SXMatrix r(x.size1(),x.size2());
  r.fill(sfcn[op](0,0));
  for(int i=0; i<r.size1(); ++i){ // loop over rows
    int el1 = x.rowind[i];
    int el2 = y.rowind[i];
    int k1 = x.rowind[i+1];
    int k2 = y.rowind[i+1];
    while(el1 < k1 || el2 < k2){
      int j1 = (el1 < k1) ? x.col[el1] : r.numel() ;
      int j2 = (el2 < k2) ? y.col[el2] : r.numel() ;
      
      if(j1==j2)
        r(i,j1) = sfcn[op](x[el1++],y[el2++]); 
      else if(j1>j2)
        r(i,j2) = sfcn[op](0,y[el2++]);
      else
        r(i,j1) = sfcn[op](x[el1++],0);
      }
    }
  return r;
}

ostream& operator<<(ostream &stream, const SXMatrix &mat){
  mat.print(stream);
  return stream;
}


} //namespace CasADi

namespace std{
using namespace CasADi;


SXMatrix sin(const SXMatrix& x){
  return unary(SIN_NODE,x);
}

SXMatrix cos(const SXMatrix& x){
  return unary(COS_NODE,x);
}

SXMatrix tan(const SXMatrix& x){
  return unary(TAN_NODE,x);
}

SXMatrix atan(const SXMatrix& x){
  return unary(ATAN_NODE,x);
}

SXMatrix asin(const SXMatrix& x){
  return unary(ASIN_NODE,x);
}

SXMatrix acos(const SXMatrix& x){
  return unary(ACOS_NODE,x);
}

SXMatrix exp(const SXMatrix& x){
  return unary(EXP_NODE,x);
}

SXMatrix log(const SXMatrix& x){
  return unary(LOG_NODE,x);
}

SXMatrix pow(const SXMatrix& x, const SXMatrix& n){
  return binary(POW_NODE,x,n);
}

SXMatrix sqrt(const SXMatrix& x){
  return unary(SQRT_NODE,x);
}

SXMatrix fmin(const SXMatrix& x, const SXMatrix& y){
  return binary(FMIN_NODE,x,y);
}

SXMatrix fmax(const SXMatrix& x, const SXMatrix& y){
  return binary(FMAX_NODE,x,y);
}

SXMatrix floor(const SXMatrix& x){
  return unary(FLOOR_NODE,x);
}

SXMatrix ceil(const SXMatrix& x){
  return unary(CEIL_NODE,x);
}

SXMatrix erf(const SXMatrix& x){
  return unary(ERF_NODE,x);
}


} // namespace std

