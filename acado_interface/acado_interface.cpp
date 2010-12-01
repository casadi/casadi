#include "acado_internal.hpp"

#include <acado_optimal_control.hpp>

#include <cassert>
#include <limits>

using namespace std;

namespace CasADi{

AcadoInterface::AcadoInterface(){ 
}

AcadoInterface::AcadoInterface(const FX& ffcn, const FX& mfcn, const FX& cfcn, const FX& rfcn){
  assignNode(new AcadoInternal(ffcn, mfcn, cfcn, rfcn));
}

AcadoInternal* AcadoInterface::operator->(){
  return (AcadoInternal*)(FX::operator->());
}

const AcadoInternal* AcadoInterface::operator->() const{
   return (const AcadoInternal*)(FX::operator->()); 
}

void AcadoInterface::assertNode() const{
  if(!dynamic_cast<const AcadoInternal*>(get()))
    throw CasadiException("AcadoInterface::assertNode");
}



} // namespace CasADi

