#include "cvodes_integrator.hpp"
#include "cvodes_internal.hpp"

using namespace std;
namespace CasADi{
namespace Sundials{

CVodesIntegrator::CVodesIntegrator(){ 
}

CVodesIntegrator::CVodesIntegrator(const FX& f,const FX& q){
  assignNode(new CVodesInternal(f,q));
}

CVodesInternal* CVodesIntegrator::operator->(){
  return (CVodesInternal*)(FX::operator->());
}

const CVodesInternal* CVodesIntegrator::operator->() const{
  return (const CVodesInternal*)(FX::operator->());
}


} // namespace Sundials
} // namespace CasADi


