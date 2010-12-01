#include "xml_parser.hpp"
#include "external_packages/tinyxml/tinyxml.h"
#include <sstream>

namespace CasADi{
namespace Modelica{

XMLParser::XMLParser(const std::string& filename){
  TiXmlDocument doc;
  bool flag = doc.LoadFile(filename.data());

  if(!flag){
    throw CasadiException("XMLParser::loadFile: Cound not open " + filename);
  }

  // parse
  document.name = filename;
  document.addNode(&doc);
}

XMLParser::~XMLParser(){

}

void XMLParser::print(std::ostream &stream) const{
  stream << document;
}


} // namespace Modelica
} // namespace CasADi
