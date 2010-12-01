#ifndef XML_ARG_HPP
#define XML_ARG_HPP

#include <string>
#include <sstream>
#include <casadi/sx/sx_matrix.hpp>
#include <casadi/mx/mx.hpp>

namespace CasADi{
namespace Modelica{

class StrArg{
  public:
  explicit StrArg(const std::string& str);

/** \brief  automatic typeconversion to various types */
  operator std::string() const;
  operator int() const;
  operator double() const;
  operator SXMatrix() const;
  operator MX() const;
  operator SX() const;

/** \brief  print */
  friend std::ostream& operator<<(std::ostream &stream, const StrArg& arg);

  std::string str;
};

} // namespace Modelica
} // namespace CasADi

#endif //XML_ARG_HPP
