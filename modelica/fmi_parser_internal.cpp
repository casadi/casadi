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

#include "fmi_parser_internal.hpp"
#include "variable_internal.hpp"
#include <map>
#include <string>
#include <sstream>
#include <cassert>
#include "casadi/stl_vector_tools.hpp"
#include "external_packages/tinyxml/tinyxml.h"

using namespace std;
namespace CasADi{
namespace Modelica{

FMIParserInternal::FMIParserInternal(const std::string& filename) : filename_(filename){
  TiXmlDocument doc;
  bool flag = doc.LoadFile(filename.data());

  if(!flag){
    throw CasadiException("XMLParser::loadFile: Cound not open " + filename);
  }

  // parse
  document_.setName(filename);
  document_.addNode(&doc);
}

FMIParserInternal::~FMIParserInternal(){

}

const OCP& FMIParserInternal::ocp() const{
  return ocp_;
}

OCP& FMIParserInternal::ocp(){
  return ocp_;
}

OCP& FMIParserInternal::parse(){	

  // Initialize function lookup tables
  unary_["Exp"]  = exp;
  unary_["Sin"]  = sin;
  unary_["Cos"]  = cos;
  unary_["Tan"]  = tan;
  unary_["Log"]  = log;
  unary_["Sqrt"]  = sqrt;
  unary_["Asin"]  = asin;
  unary_["Acos"]  = acos;
  unary_["Atan"]  = atan;
  
  binary_["Pow"]  = pow;

  // Clear the ocp
  ocp_ = OCP();
  
 // Add model variables
  addModelVariables();
  
  // Add binding equations
  addBindingEquations();

  // Add dynamic equations
  addDynamicEquations();

  // Add initial equations
  addInitialEquations();

  // Add optimization
  if(document_[0].hasChild("opt:Optimization"))
    addOptimization();

  // Return a reference to the created ocp
  return ocp_;

}

void FMIParserInternal::addModelVariables(){
cout << "vars" << endl;

  map<string,int> catCounter;
  
  // Time
  Variable time("time");
  time.setExpression(SX("time"));
  time.setIndependent(true);
  ocp_.variables->add(time);
  
  // Get a reference to the ModelVariables node
  const XMLNode& modvars = document_[0]["ModelVariables"];

  // Add variables
  for(int i=0; i<modvars.size(); ++i){
    
    // Get a reference to the variable
    const XMLNode& vnode = modvars[i];

    // Get the attributes
    string name        = vnode.attribute("name");
    int valueReference = vnode.attribute("valueReference");
    string variability = vnode.attribute("variability");
    string causality   = vnode.attribute("causality");
    string alias       = vnode.attribute("alias");
    
    // Skip to the next variable if its an alias
    if(alias.compare("alias") == 0 || alias.compare("negatedAlias") == 0)
      continue;
    
    if(modvars[i].hasChild("VariableCategory")){
      string variableCategory = modvars[i]["VariableCategory"].getText();
      map<string,int>::iterator it=catCounter.find(variableCategory);
      if(it!=catCounter.end())
        it->second++;
      else
        catCounter[variableCategory] = 1;
    }
    
    // Add to ocp
    Variable var = ocp_.variables;
    
    //cout << "adding " << var << endl;
    const XMLNode& nn = vnode["QualifiedName"];
    for(int i=0; i<nn.size(); ++i){
      string namepart = nn[i].attribute("name");
      if(!var->has(namepart)){
        stringstream fullname;
        fullname << var << "." << namepart;
        var->add(Variable(fullname.str()),namepart);
      }
      
      var = var(namepart);

      // Add index
      if(nn[i].size()>0){
        // Find the index
        int ind = nn[i]["exp:ArraySubscripts"]["exp:IndexExpression"]["exp:IntegerLiteral"].getText();
        
        // Allocate size
        if(var->col.size()<ind)
          var->col.resize(ind);
       
        if(var[ind].isNull()){
          stringstream fullname;
          fullname << var << "[" << ind << "]";
          var->col.at(ind-1) = Variable(fullname.str());
        }

        // Go to the sub-collection
        var = var[ind];
      }
    }

    // Set fields if not already added
    if(var.getExpression()->isNan()){
      
      // Expression
      var.setExpression(SX(name));

      // Symbolic variable
      var.setName(name);

      // Value reference
      var.setValueReference(valueReference);
      
      // Variability
      if(variability.compare("constant")==0)
        var.setVariability(CONSTANT);
      else if(variability.compare("parameter")==0)
        var.setVariability(PARAMETER);
      else if(variability.compare("discrete")==0)
        var.setVariability(DISCRETE);
      else if(variability.compare("continuous")==0)
        var.setVariability(CONTINUOUS);
      else throw CasadiException("Unknown variability");
  
      // Causality
      if(causality.compare("input")==0)
        var.setCausality(INPUT);
      else if(causality.compare("output")==0)
        var.setCausality(OUTPUT);
      else if(causality.compare("internal")==0)
        var.setCausality(INTERNAL);
      else throw CasadiException("Unknown causality");
      
      // Alias
      if(alias.compare("noAlias")==0)
        var.setAlias(NO_ALIAS);
      else if(alias.compare("alias")==0)
        var.setAlias(ALIAS);
      else if(alias.compare("negatedAlias")==0)
        var.setAlias(NEGATED_ALIAS);
      else throw CasadiException("Unknown alias");
      

      // Other properties
      const XMLNode& props = vnode[0];
      if(props.hasAttribute("unit"))         var.setUnit(props.attribute("unit"));
      if(props.hasAttribute("displayUnit"))  var.setDisplayUnit(props.attribute("displayUnit"));
      if(props.hasAttribute("min"))          var.setMin(props.attribute("min"));
      if(props.hasAttribute("max"))          var.setMax(props.attribute("max"));
      if(props.hasAttribute("start"))        var.setStart(props.attribute("start"));
      if(props.hasAttribute("nominal"))      var.setNominal(props.attribute("nominal"));

    } else {
      // Mark the variable differentiated
      var.setDerivative(SX(name)); // NOTE: assumes that the derivative is added after the variable
      
    }
  }
  
  ocp_.variables->print(cout,0);
  
  cout << "VariableCategories:" << endl;
  for(map<string,int>::const_iterator it=catCounter.begin(); it!=catCounter.end(); ++it)
    cout << "  " << it->first << ": " << it->second << endl;
}


void FMIParserInternal::addBindingEquations(){
  // Get a reference to the BindingEquations node
  const XMLNode& bindeqs = document_[0]["equ:BindingEquations"];
  
  for(int i=0; i<bindeqs.size(); ++i){
    const XMLNode& beq = bindeqs[i];

    // Get the variable
    Variable var = readVariable(beq[0]);

    // Get the binding equation
    SX bexpr = readExpr_new(beq[1][0]);

    // Pass to variable
    var.setBindingEquation(bexpr);
  }
}

void FMIParserInternal::addDynamicEquations(){
  // Get a reference to the DynamicEquations node
  const XMLNode& dyneqs = document_[0]["equ:DynamicEquations"];

  // Add equations
  for(int i=0; i<dyneqs.size(); ++i){

    // Get a reference to the variable
    const XMLNode& dnode = dyneqs[i];

    // Add the differential equation
    SX de_new = readExpr_new(dnode[0]);
    ocp_.dae.push_back(de_new);
  }
}

void FMIParserInternal::addInitialEquations(){
  // Get a reference to the DynamicEquations node
  const XMLNode& initeqs = document_[0]["equ:InitialEquations"];

  // Add equations
  for(int i=0; i<initeqs.size(); ++i){

    // Get a reference to the node
    const XMLNode& inode = initeqs[i];

    // Add the differential equations
    for(int i=0; i<inode.size(); ++i){
      ocp_.initeq.push_back(readExpr_new(inode[i]));
    }
  }
}

void FMIParserInternal::addOptimization(){
  // Get a reference to the DynamicEquations node
  const XMLNode& opts = document_[0]["opt:Optimization"];
  
  // Start time
  ocp_.t0  = opts["opt:IntervalStartTime"]["opt:Value"].getText();

  // Terminal time
  ocp_.tf = opts["opt:IntervalFinalTime"]["opt:Value"].getText();

  for(int i=0; i<opts.size(); ++i){
    
    // Get a reference to the node
    const XMLNode& onode = opts[i];

    // Get the type
    if(onode.checkName("opt:ObjectiveFunction")){ // mayer term
      addObjectiveFunction(onode);
    } else if(onode.checkName("opt:IntegrandObjectiveFunction")){
      addIntegrandObjectiveFunction(onode);
    } else if(onode.checkName("opt:IntervalStartTime")) {
       addIntervalStartTime(onode);
    } else if(onode.checkName("opt:IntervalFinalTime")) {
       addIntervalFinalTime(onode);
    } else if(onode.checkName("opt:TimePoints")) {
//       addTimePoints(onode);
    } else if(onode.checkName("opt:Constraints")) {
      addConstraints(onode);
    } else throw "FMIParserInternal::addOptimization: Unknown node";
  }
}

void FMIParserInternal::addObjectiveFunction(const XMLNode& onode){
  // Add components
  for(int i=0; i<onode.size(); ++i){
    const XMLNode& var = onode[i];

    if(var.checkName("exp:TimedVariable")){
      SX v = readExpr_new(var[0]);
      SX tp = readExpr_new(var["exp:Instant"]);
      ocp_.mterm.push_back(v);
      ocp_.mtp.push_back(tp->getValue());  // always a constant???
    }
  }
}

void FMIParserInternal::addIntegrandObjectiveFunction(const XMLNode& onode){
  for(int i=0; i<onode.size(); ++i){
    const XMLNode& var = onode[i];
    SX v = readExpr_new(var[0]);
    ocp_.lterm.push_back(v);
  }
}

void FMIParserInternal::addIntervalStartTime(const XMLNode& onode){
  
}

void FMIParserInternal::addIntervalFinalTime(const XMLNode& onode){
  
}


void FMIParserInternal::addConstraints(const XMLNode& onode){
  for(int i=0; i<onode.size(); ++i){

    const XMLNode& constr_i = onode[i];
    if(constr_i.checkName("opt:ConstraintLeq")){
      SX ex = readExpr_new(constr_i[0]);
      SX ub = readExpr_new(constr_i[1]);
      ocp_.cfcn.push_back(ex);
      ocp_.cfcn_lb.push_back(-numeric_limits<double>::infinity());
      ocp_.cfcn_ub.push_back(ub);
    } else if(constr_i.checkName("opt:ConstraintGeq")){
      SX ex = readExpr_new(constr_i[0]);
      SX lb = readExpr_new(constr_i[1]);
      ocp_.cfcn.push_back(ex);
      ocp_.cfcn_lb.push_back(lb);
      ocp_.cfcn_ub.push_back(numeric_limits<double>::infinity());
    } else {
      cerr << "unknown constraint type" << constr_i.getName() << endl;
      throw "FMIParserInternal::addConstraints";
    }
  }
}

Variable FMIParserInternal::readVariable(const XMLNode& node) const{
  // Get a pointer to the variable collection
  Variable col = ocp_.variables;
  
  // Find the variable
  for(int i=0; i<node.size(); ++i){
    // Go to the right variable
    string namepart = node[i].attribute("name");
    col = col(namepart);
    
    // Go to the right index, if there is any
    if(node[i].size()>0){
      // Go to the sub-collection
      int ind = node[i]["exp:ArraySubscripts"]["exp:IndexExpression"]["exp:IntegerLiteral"].getText();
      col = col[ind];
    }
  }
  
  // Return the variables
  return col;
}


SX FMIParserInternal::readExpr_new(const XMLNode& node){
  const string& fullname = node.getName();
  if (fullname.find("exp:")== string::npos) {
    stringstream ss;
    ss << "FMIParserInternal::readExpr: unknown - expression is supposed to start with 'exp:' , got " << fullname;
    throw CasadiException(ss.str());
  }
  
  // Chop the 'exp:'
  string name = fullname.substr(4);

  // The switch below is alphabetical, and can be thus made more efficient
  if(name.compare("Add")==0){
    return readExpr_new(node[0]) + readExpr_new(node[1]);
  } else if(name.compare("Der")==0){
    return readVariable(node[0]).der();
  } else if(name.compare("Div")==0){
    return readExpr_new(node[0]) / readExpr_new(node[1]);
  } else if(name.compare("Identifier")==0){
    return readVariable(node);
  } else if(name.compare("IntegerLiteral")==0){
    return int(node.getText());
  } else if(name.compare("Instant")==0){
    return node.getText();
  } else if(name.compare("LogLt")==0){ // Logical less than
    return readExpr_new(node[0]) < readExpr_new(node[1]);
  } else if(name.compare("Mul")==0){ // Multiplication
    return readExpr_new(node[0]) * readExpr_new(node[1]);
  } else if(name.compare("Neg")==0){
    return -readExpr_new(node[0]);
  } else if(name.compare("NoEvent")==0) {
    // NOTE: This is a workaround, we assume that whenever NoEvent occurs, what is meant is a switch
    int n = node.size();
    
    // Default-expression
    SX ex = readExpr_new(node[n-1]);
    
    // Evaluate ifs
    for(int i=n-3; i>=0; i -= 2) ex = if_else(readExpr_new(node[i]),readExpr_new(node[i+1]),ex);
    
    return ex;
  } else if(name.compare("RealLiteral")==0){
    return double(node.getText());
  } else if(name.compare("Sub")==0){
    return readExpr_new(node[0]) - readExpr_new(node[1]);
  } else if(name.compare("Time")==0){
    return ocp_.variables("time");
  }

  // Check if it is a unary function
  map<string,SX (*)(const SX&)>::iterator uit;
  uit = unary_.find(name);
  if (uit!=unary_.end()) {
    // Make sure that there is exactly one child
    if ((node.size()) != 1)
      throw CasadiException(string("FMIParserInternal::readExpr: Node \"") + fullname + "\" does not have one child");

    // Call the function
    return (*uit->second)(readExpr_new(node[0]));
  }

  // Check if it is a binary function
  map<string,SX (*)(const SX&, const SX&)>::iterator bit;
  bit = binary_.find(name);
  if (bit!=binary_.end()) {
    // Make sure that there are exactly two children
    if ((node.size()) != 2)
      throw CasadiException(string("FMIParserInternal::readExpr: Node \"") + fullname + "\" does not have two children");

    // Call the function
    return (*bit->second)(readExpr_new(node[0]),readExpr_new(node[1]));
  }

  // throw error if reached this point
  throw CasadiException(string("FMIParserInternal::readExpr: unknown node: ") + name);
  
}

void FMIParserInternal::repr(std::ostream &stream) const{
  stream << "FMI parser (XML file: \"" << filename_ << "\")";
}

void FMIParserInternal::print(std::ostream &stream) const{
  stream << document_;
}

} // namespace Modelica
} // namespace CasADi
