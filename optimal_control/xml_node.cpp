/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "xml_node.hpp"
#include "external_packages/tinyxml/tinyxml.h"
#include "../casadi/casadi_exception.hpp"

using namespace std;
namespace CasADi{
namespace OptimalControl{

XMLNode::XMLNode(){
}

XMLNode::XMLNode(const string& name) : name(name) {  
}

XMLNode::~XMLNode(){
  // delete all children
  for(int i=0; i<children.size(); ++i)
    delete children[i];
}

bool XMLNode::hasAttribute(const string& attribute_name) const{
  map<string, string>::const_iterator it = attributes.find(attribute_name);
  return it!=attributes.end();
}

StrArg XMLNode::attribute(const string& attribute_name) const{
  // find the attribute
  map<string, string>::const_iterator it = attributes.find(attribute_name);

  // check that the attribute was indeed found
  if(it == attributes.end()){
    throw CasadiException("Error in XMLNode::attribute: could not find " + attribute_name);
  }

  // Return an instance of the XMLArg class (that automatically typeconverts into other types)
  return StrArg(it->second);
}

XMLNode& XMLNode::operator[](int i) const{
  casadi_assert_message(i>=0 && i < size(), "XMLNode::operator[]: index out of bounds for element " << i << " of node " << getName());
  
  return *children.at(i);
}

bool XMLNode::hasChild(const string& childname) const{
  map<string,int>::const_iterator it = child_indices.find(childname);
  return it!=child_indices.end();
}
  

XMLNode& XMLNode::operator[](const string& childname) const{
  // Find the child
  map<string,int>::const_iterator it = child_indices.find(childname);

  // check that the child was indeed found
  if(it == child_indices.end()){
    throw CasadiException("Error in XMLNode::operator[]: could not find " + childname);
  }

  // Return an index to the child
  return *children[it->second];
}

void XMLNode::setAttribute(const string& attribute_name, const string& attribute){
  attributes[attribute_name] = attribute;
}

void XMLNode::addChild(XMLNode* child){
  children.push_back(child);
  child_indices[child->name] = children.size()-1;
}

void XMLNode::addAttributes(TiXmlElement* el){

  if ( !el) return;
  for(TiXmlAttribute* pAttrib=el->FirstAttribute(); pAttrib; pAttrib=pAttrib->Next()){
      setAttribute(pAttrib->Name(),pAttrib->Value());
   }
}

void XMLNode::addNode(TiXmlNode* n){
  if (!n) throw CasadiException("Error in XMLNode::addNode: Node is 0");

  // Save attributes
  int type = n->Type();  
  if(type == TiXmlNode::ELEMENT){
    addAttributes(n->ToElement());
  } else if(type == TiXmlNode::DOCUMENT) {
    // do nothing
  } else {
    throw CasadiException("XMLNode::addNode");
  }

  // add children
  for ( TiXmlNode* child = n->FirstChild(); child != 0; child= child->NextSibling()){ 
      int childtype = child->Type();  

      if(childtype == TiXmlNode::ELEMENT){
	XMLNode *newnode = new XMLNode(child->Value());
	newnode->addNode(child);
	addChild(newnode);
      } else if(childtype == TiXmlNode::COMMENT){
	comment = child->Value();
      } else if(childtype == TiXmlNode::TEXT){
	text = child->ToText()->Value();
      } else if (childtype == TiXmlNode::DECLARATION){
	cout << "Warning: Skipped TiXmlNode::DECLARATION" << endl;
      } else {
	throw CasadiException("Error in XMLNode::addNode: Unknown node type");
      }
  }
}

ostream& operator<<(ostream &stream, const XMLNode& node){
  node.dump(stream);
  return stream;
}

int XMLNode::size() const{
  return children.size();
}

const string& XMLNode::getName() const{
  return name;
}

void XMLNode::setName(const string& name){
  this->name = name;
}

void XMLNode::dump(ostream &stream, int indent) const{
  // Print name
  stream << string( indent,' ') << "Node: " << name << endl;

  // Print comment
  if(!comment.empty()){
    stream << string( indent,' ') << "----- comment starts ----- "  << endl;
    stream << comment << endl;
    stream << string( indent,' ') << "----- comment ends ----- "  << endl;
  }

  // Print text
  if(!text.empty())
    stream << string( indent+2,' ') << "Text: " << text << endl;

  // Print attributes
  for(map<string, string>::const_iterator it=attributes.begin(); it != attributes.end(); ++it)
    stream << string( indent+2,' ') << "attribute " << it->first << " = " << it->second << endl;

  // Print Children
  for(int i=0; i<size(); ++i){
    stream << string( indent,' ') << "Child " << i << ":" << endl;
    (*this)[i].dump(stream,indent+2);
  }
}

StrArg XMLNode::getText() const{
  return StrArg(text);
}

bool XMLNode::checkName(const string& str) const{
  return name.compare(str) == 0;
}

} // namespace OptimalControl
} // namespace CasADi
