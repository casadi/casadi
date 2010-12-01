#include "idas_integrator.hpp"
#include "idas_internal.hpp"

using namespace std;
namespace CasADi{
namespace Sundials{

IdasIntegrator::IdasIntegrator(){ 
}

IdasIntegrator::IdasIntegrator(const FX& f, const FX& q, const FX& jacx, const FX& jacp){
  assignNode(new IdasInternal(f,q,jacx,jacp));
}

IdasInternal* IdasIntegrator::operator->(){
  return (IdasInternal*)(FX::operator->());
}

const IdasInternal* IdasIntegrator::operator->() const{
  return (const IdasInternal*)(FX::operator->());
}
} // namespace Sundials
} // namespace CasADi


