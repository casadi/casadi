#include "nlp_solver.hpp"
#include "nlp_solver_internal.hpp"

namespace CasADi{

NLPSolver::NLPSolver(){
}

NLPSolverInternal* NLPSolver::operator->(){
  return (NLPSolverInternal*)(FX::operator->());
}

const NLPSolverInternal* NLPSolver::operator->() const{
  return (const NLPSolverInternal*)(FX::operator->());
}
    
void NLPSolver::assertNode() const{
  if(!dynamic_cast<const NLPSolverInternal*>(get()))
    throw CasadiException("NLPSolver::assertNode");
}


} // namespace CasADi

