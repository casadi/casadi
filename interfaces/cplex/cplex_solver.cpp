#include "cplex_solver.hpp"
#include "cplex_internal.hpp"

using namespace std;

namespace CasADi{

CplexSolver::CplexSolver(){
}

CplexSolver::CplexSolver(const FX& F, const FX& G, const FX& H, const FX& J, const FX& GF){
  assignNode(new CplexInternal(F,G,H,J,GF));
}

CplexInternal* CplexSolver::operator->(){
  return (CplexInternal*)(NLPSolver::operator->());
}

const CplexInternal* CplexSolver::operator->() const{
  return (const CplexInternal*)(NLPSolver::operator->());
}

bool CplexSolver::checkNode() const{
  return dynamic_cast<const CplexInternal*>(get());
}

void CplexSolver::setIntParam(const std::string& name, int val){
  (*this)->int_param_[name] = val;
}

void CplexSolver::setDoubleParam(const std::string& name, double val){
  (*this)->double_param_[name] = val;
}


} // namespace CasADi
