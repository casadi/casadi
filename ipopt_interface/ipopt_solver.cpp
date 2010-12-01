#include "ipopt_internal.hpp"
#include "ipopt_nlp.hpp"

using namespace std;

#include <coin/IpIpoptApplication.hpp>
namespace CasADi{

IpoptSolver::IpoptSolver(){
}
  
IpoptSolver::IpoptSolver(const FX& F, const FX& G, const FX& H, const FX& J){
  assignNode(new IpoptInternal(F,G,H,J));
}

IpoptInternal* IpoptSolver::operator->(){
  return (IpoptInternal*)(NLPSolver::operator->());
}

const IpoptInternal* IpoptSolver::operator->() const{
  return (const IpoptInternal*)(NLPSolver::operator->());
}
    
void IpoptSolver::assertNode() const{
  if(!dynamic_cast<const IpoptInternal*>(get()))
    throw CasadiException("IpoptSolver::assertNode");
}

} // namespace CasADi
