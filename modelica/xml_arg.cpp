#include "xml_arg.hpp"
#include "casadi/expression_tools.hpp"

namespace CasADi{
namespace Modelica{

StrArg::StrArg(const std::string& str) : str(str){
}

StrArg::operator std::string() const{
  return str;
}

StrArg::operator int() const{
  std::istringstream buffer(str);
  int ret;
  buffer >> ret;
  return ret;
}

StrArg::operator double() const{
  std::istringstream buffer(str);
  double ret;
  buffer >> ret;
  return ret;
}

StrArg::operator SXMatrix() const{
  std::istringstream buffer(str);
  SXMatrix ret;
  buffer >> ret;
  return ret;
}

StrArg::operator SX() const{
  SXMatrix m = *this;
  return SX(m.getElement());
}

StrArg::operator MX() const{
  return MX(str);
}


std::ostream& operator<<(std::ostream &stream, const StrArg& arg){
  return stream << arg.str;
}

} // namespace Modelica
} // namespace CasADi

