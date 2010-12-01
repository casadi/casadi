#include "norm.hpp"

namespace CasADi{

Norm::Norm(const MX& x) : MXNode(x){
  sz.nrow = 1;
  sz.ncol = 1;
}

void Norm::evaluate(int fsens_order, int asens_order){
  throw CasadiException("Norm::evaluate not implemented");
}

Norm2::Norm2(const MX& x) : Norm(x){
}

Norm2* Norm2::clone() const{
  return new Norm2(*this);
}

void Norm2::print(std::ostream &stream) const{
  stream << "||" << dep(0) << "||_2"; 
}

Norm1::Norm1(const MX& x) : Norm(x){
}

Norm1* Norm1::clone() const{
  return new Norm1(*this);
}

void Norm1::print(std::ostream &stream) const{
  stream << "||" << dep(0) << "||_1"; 
}

NormInf::NormInf(const MX& x) : Norm(x){
}

NormInf* NormInf::clone() const{
  return new NormInf(*this);
}

void NormInf::print(std::ostream &stream) const{
  stream << "||" << dep(0) << "||_inf"; 
}

} // namespace CasADi

