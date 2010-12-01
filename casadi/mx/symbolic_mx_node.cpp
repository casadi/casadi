#include "symbolic_mx_node.hpp"
#include <cassert>
#include <vector>

using namespace std;

namespace CasADi{

SymbolicMatrix::SymbolicMatrix(const std::string& name_, int n, int m) : name(name_) {
  sz = MatrixSize(n,m);
}

SymbolicMatrix* SymbolicMatrix::clone() const{
  return new SymbolicMatrix(*this);
}

void SymbolicMatrix::print(std::ostream &stream) const{
  stream << name;
}

void SymbolicMatrix::evaluate(int fsens_order, int asens_order){
}

// void SymbolicMatrix::evaluateAdj(){
// }

bool SymbolicMatrix::isSymbolic() const{
  return true;
}

const std::string& SymbolicMatrix::getName() const{
  return name;
}

} // namespace CasADi

