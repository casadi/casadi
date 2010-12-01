#ifndef XML_PARSER_HPP
#define XML_PARSER_HPP

#include <iostream>
#include <fstream>
#include <string>

#include "xml_node.hpp"
#include "casadi/printable_object.hpp"

namespace CasADi{
namespace Modelica{

/** \brief  Forward declarations */
class TiXmlElement;
class TiXmlNode;

class XMLParser : public PrintableObject{

public:
  XMLParser(const std::string& filename);  // constructor
  virtual ~XMLParser()=0; // destructor

  virtual void print(std::ostream &stream=std::cout) const;
  
protected:

  XMLNode document;

};

} // namespace Modelica
} // namespace CasADi

#endif //XML_PARSER_HPP
