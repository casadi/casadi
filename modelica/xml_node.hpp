#ifndef XML_NODE_HPP
#define XML_NODE_HPP

#include <string>
#include <vector>
#include <map>
#include "xml_arg.hpp"

/** \brief  Forward declarations */
class TiXmlElement;
class TiXmlNode;

namespace CasADi{
namespace Modelica{

class XMLNode{
friend class XMLParser;
public:
  XMLNode();
  XMLNode(const std::string& name);
  ~XMLNode();

/** \brief  Add an attribute */
  void setAttribute(const std::string& attribute_name, const std::string& attribute);

/** \brief  Add a child */
  void addChild(XMLNode *child);

/** \brief  Get an attribute by its name */
  StrArg attribute(const std::string& attribute_name) const;

/** \brief  Get a reference to a child by its index */
  XMLNode& operator[](int i) const;

/** \brief  Get a reference to a child by its name */
  XMLNode& operator[](const std::string& childname) const;

/** \brief  Check if a child is present */
  bool hasChild(const std::string& childname) const;
  
/** \brief  Check if an attribute is present */
  bool hasAttribute(const std::string& attribute_name) const;

/** \brief  Get the number of children */
  int size() const;

/** \brief  Get the name of the node */
  const std::string& getName() const;

/** \brief  check if the name is equal to something */
  bool checkName(const std::string& str) const;

/** \brief  Get the value of the text field */
  StrArg getText() const;

  void addAttributes(TiXmlElement* el);
  void addNode(TiXmlNode* node);

  friend std::ostream& operator<<(std::ostream &stream, const XMLNode& node);

  void dump(std::ostream &stream, int indent=0) const;

  protected:
/** \brief  Attributes and children (binary search tree) */
  std::map<std::string, std::string>  attributes;
  std::vector<XMLNode*>               children;
  std::map<std::string,int>           child_indices; // the index of the children sorted by their name

  std::string name;
  std::string comment;
  std::string text;

};

} // namespace Modelica
} // namespace CasADi


#endif //XML_NODE_HPP
