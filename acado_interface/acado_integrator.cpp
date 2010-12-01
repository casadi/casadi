#include "acado_integrator.hpp"
#include "acado_integrator_internal.hpp"

using namespace std;
namespace CasADi{

ACADOIntegrator::ACADOIntegrator(){ 
}

ACADOIntegrator::ACADOIntegrator(const FX& f){
  assignNode(new ACADOIntegratorInternal(f));
}

ACADOIntegratorInternal* ACADOIntegrator::operator->(){
  return (ACADOIntegratorInternal*)(FX::operator->());
}

const ACADOIntegratorInternal* ACADOIntegrator::operator->() const{
  return (const ACADOIntegratorInternal*)(FX::operator->());
}


} // namespace CasADi


