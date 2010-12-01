#include "integrator.hpp"
#include "integrator_internal.hpp"
#include <cassert>

namespace CasADi{

  
  
Integrator::Integrator(){
}

void Integrator::printStats(ostream &stream) const{
  (*this)->printStats(stream);
}
  
IntegratorInternal* Integrator::operator->(){
  return (IntegratorInternal*)(FX::operator->());
}

const IntegratorInternal* Integrator::operator->() const{
   return (const IntegratorInternal*)(FX::operator->()); 
}
  
void Integrator::reset(int fsens_order, int asens_order){
  (*this)->reset(fsens_order, asens_order);
}

void Integrator::integrate(double t_out){
  (*this)->integrate(t_out);
}
  
void Integrator::setStopTime(double tf){
  (*this)->setStopTime(tf);
}
  
void Integrator::assertNode() const{
  if(!dynamic_cast<const IntegratorInternal*>(get()))
    throw CasadiException("Integrator::assertNode");
}
  
 
} // namespace CasADi

