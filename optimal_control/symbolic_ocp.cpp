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

#include "symbolic_ocp.hpp"
#include "variable_internal.hpp"
#include "xml_node.hpp"

#include <map>
#include <string>
#include <sstream>
#include <ctime>

#include "casadi/stl_vector_tools.hpp"
#include "external_packages/tinyxml/tinyxml.h"

#include "../casadi/casadi_exception.hpp"
#include "../casadi/stl_vector_tools.hpp"
#include "variable_tools.hpp"
#include "../casadi/matrix/matrix_tools.hpp"
#include "../casadi/sx/sx_tools.hpp"
#include "../casadi/fx/integrator.hpp"
#include "../casadi/casadi_calculus.hpp"

using namespace std;
namespace CasADi{

SymbolicOCP::SymbolicOCP(){
  t = ssym("t");
  t0 = numeric_limits<double>::quiet_NaN();
  tf = numeric_limits<double>::quiet_NaN();
  t0_free = false;
  tf_free = false;
}

void SymbolicOCP::parseFMI(const std::string& filename, const Dictionary& options){
  // Default options
  bool verbose = false;
  bool scale_variables = false;
  bool scale_equations = false;
  bool eliminate_dependent = true;
  bool sort_equations = true;
  bool make_explicit = false;
  
  // Read user options
  for(Dictionary::const_iterator it=options.begin(); it!=options.end(); ++it){
    if(it->first.compare("verbose")==0){
      verbose = it->second;
    } else if(it->first.compare("scale_variables")==0){
      scale_variables = it->second;
    } else if(it->first.compare("scale_equations")==0){
      scale_equations = it->second;
    } else if(it->first.compare("eliminate_dependent")==0){
      eliminate_dependent = it->second;
    } else if(it->first.compare("sort_equations")==0){
      sort_equations = it->second;
    } else if(it->first.compare("make_explicit")==0){
      make_explicit = it->second;
    } else {
      stringstream ss;
      ss << "Unknown option \"" << it->first << "\"" << endl;
      throw CasadiException(ss.str());
    }
  }
  
  // Load 
  TiXmlDocument doc;
  bool flag = doc.LoadFile(filename.c_str());
  casadi_assert_message(flag, "Cound not open " << filename);

  // parse
  XMLNode document;
  document.setName(filename);
  document.addNode(&doc);

  double time1 = clock();

  // **** Add model variables ****
  {
    if(verbose) cout << "Adding model variables." << endl;
  
    // Get a reference to the ModelVariables node
    const XMLNode& modvars = document[0]["ModelVariables"];

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
	  
      // Get the name
      const XMLNode& nn = vnode["QualifiedName"];
      string qn = qualifiedName(nn);
      
      // Find variable
      map<string,Variable>::iterator it = varmap_.find(qn);
      
      // Add variable, if not already added
      if(it == varmap_.end()){
	
	// Create variable
	Variable var(name);
	
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
	if(vnode.hasChild("Real")){
	  const XMLNode& props = vnode["Real"];
	  if(props.hasAttribute("unit"))         var.setUnit(props.attribute("unit"));
	  if(props.hasAttribute("displayUnit"))  var.setDisplayUnit(props.attribute("displayUnit"));
	  if(props.hasAttribute("min"))          var.setMin(props.attribute("min"));
	  if(props.hasAttribute("max"))          var.setMax(props.attribute("max"));
	  if(props.hasAttribute("start"))        var.setStart(props.attribute("start"));
	  if(props.hasAttribute("nominal"))      var.setNominal(props.attribute("nominal"));
	  if(props.hasAttribute("free"))         var.setFree(string(props.attribute("free")).compare("true") == 0);
	  if(props.hasAttribute("initialGuess")) var.setInitialGuess(props.attribute("initialGuess"));
	}
	
	// Variable category
	if(vnode.hasChild("VariableCategory")){
	  string cat = vnode["VariableCategory"].getText();
	  if(cat.compare("derivative")==0)
	    var.setCategory(CAT_DERIVATIVE);
	  else if(cat.compare("state")==0)
	    var.setCategory(CAT_STATE);
	  else if(cat.compare("dependentConstant")==0)
	    var.setCategory(CAT_DEPENDENT_CONSTANT);
	  else if(cat.compare("independentConstant")==0)
	    var.setCategory(CAT_INDEPENDENT_CONSTANT);
	  else if(cat.compare("dependentParameter")==0)
	    var.setCategory(CAT_DEPENDENT_PARAMETER);
	  else if(cat.compare("independentParameter")==0)
	    var.setCategory(CAT_INDEPENDENT_PARAMETER);
	  else if(cat.compare("algebraic")==0)
	    var.setCategory(CAT_ALGEBRAIC);
	  else throw CasadiException("Unknown variable category: " + cat);
	}
	
	// Add to list of variables
	addVariable(qn,var);
      }
    }
  }
  
  // **** Add binding equations ****
  {
    if(verbose) cout << "Adding binding equations." << endl;
  
    // Get a reference to the BindingEquations node
    const XMLNode& bindeqs = document[0]["equ:BindingEquations"];
  
    for(int i=0; i<bindeqs.size(); ++i){
      const XMLNode& beq = bindeqs[i];

      // Get the variable
      Variable var = readVariable(beq[0]);

      // Get the binding equation
      SX bexpr = readExpr(beq[1][0]);
      
      // Add binding equation
      var.setBinding(bexpr);
      y.push_back(var); // legacy
      dep.append(bexpr); // legacy
    }
  }

  // **** Add dynamic equations ****
  {
    // Get a reference to the DynamicEquations node
    const XMLNode& dyneqs = document[0]["equ:DynamicEquations"];

    // Add equations
    for(int i=0; i<dyneqs.size(); ++i){

      // Get a reference to the variable
      const XMLNode& dnode = dyneqs[i];

      // Add the differential equation
      SX de_new = readExpr(dnode[0]);
      dae.append(de_new);
    }
  }
  
  // **** Add initial equations ****
  {
    // Get a reference to the DynamicEquations node
    const XMLNode& initeqs = document[0]["equ:InitialEquations"];

    // Add equations
    for(int i=0; i<initeqs.size(); ++i){

      // Get a reference to the node
      const XMLNode& inode = initeqs[i];

      // Add the differential equations
      for(int i=0; i<inode.size(); ++i){
	initial.append(readExpr(inode[i]));
      }
    }
  }
  
  // **** Add optimization ****
  if(document[0].hasChild("opt:Optimization")){
    
    // Get a reference to the DynamicEquations node
    const XMLNode& opts = document[0]["opt:Optimization"];
    
    // Start time
    t0  = opts["opt:IntervalStartTime"]["opt:Value"].getText();

    // Terminal time
    tf = opts["opt:IntervalFinalTime"]["opt:Value"].getText();

    for(int i=0; i<opts.size(); ++i){
      
      // Get a reference to the node
      const XMLNode& onode = opts[i];

      // Get the type
      if(onode.checkName("opt:ObjectiveFunction")){ // mayer term
	try{
	  // Add components
	  for(int i=0; i<onode.size(); ++i){
	    const XMLNode& var = onode[i];
	    SX v = readExpr(var);
	    mterm.append(v);
	  }
	} catch(exception& ex){
	  if(verbose){
	    cout << "WARNING: addObjectiveFunction" << ex.what() << endl;
	  }
	}
      } else if(onode.checkName("opt:IntegrandObjectiveFunction")){
	try{
	  for(int i=0; i<onode.size(); ++i){
	    const XMLNode& var = onode[i];
	    SX v = readExpr(var);
	    lterm.append(v);
	  }
	} catch(exception& ex){
	  cout << "WARNING: addIntegrandObjectiveFunction" << ex.what() << endl;
	}
      } else if(onode.checkName("opt:IntervalStartTime")) {
	// TODO
      } else if(onode.checkName("opt:IntervalFinalTime")) {
	// TODO
      } else if(onode.checkName("opt:TimePoints")) {
	// TODO
      } else if(onode.checkName("opt:Constraints")) {
	
	for(int i=0; i<onode.size(); ++i){

	  const XMLNode& constr_i = onode[i];
	  if(constr_i.checkName("opt:ConstraintLeq")){
	    SX ex = readExpr(constr_i[0]);
	    SX ub = readExpr(constr_i[1]);
	    path.append(ex-ub);
	    path_min.append(-numeric_limits<double>::infinity());
	    path_max.append(0.);
	  } else if(constr_i.checkName("opt:ConstraintGeq")){
	    SX ex = readExpr(constr_i[0]);
	    SX lb = readExpr(constr_i[1]);
	    path.append(ex-lb);
	    path_min.append(0.);
	    path_max.append(numeric_limits<double>::infinity());
	  } else if(constr_i.checkName("opt:ConstraintEq")){
	    SX ex = readExpr(constr_i[0]);
	    SX eq = readExpr(constr_i[1]);
	    path.append(ex-eq);
	    path_min.append(0.);
	    path_max.append(0.);
	  } else {
	    cerr << "unknown constraint type" << constr_i.getName() << endl;
	    throw "SymbolicOCP::addConstraints";
	  }
	}
	
      } else throw "SymbolicOCP::addOptimization: Unknown node";
    }
  }
  
  // Sort the variables according to type
  sortType();
  
  // Make sure that the dimensions are consistent at this point
  casadi_assert(x.size()==dae.size());
  casadi_assert(xd.size()==ode.size());
  casadi_assert(xa.size()==alg.size());
  casadi_assert(xq.size()==quad.size());
  casadi_assert(y.size()==dep.size());
  
  // Return a reference to the created ocp
  double time2 = clock();
  double tparse = double(time2-time1)/CLOCKS_PER_SEC;
  if(verbose) cout << "... parsing complete after " << tparse << " seconds" << endl;

  // Scale the variables
  if(scale_variables)
    scaleVariables();

  if(eliminate_dependent){
    // Eliminate interdependencies
    eliminateInterdependencies();
  
    // Eliminate the dependent variables
    eliminateDependent();
  }

  // Sort the equations and variables
  if(sort_equations){
    casadi_assert(eliminate_dependent);
    sortDAE();
  }

  // Scale the equations
  if(scale_equations)
    scaleEquations();
    
  // Transform to semi-explicit form
  if(make_explicit){
    casadi_assert(eliminate_dependent);
    
    // Save the old number of dependent states
    int ny = y.size();
    
    // Solve for the highest order derivatives
    makeExplicit();
    
    // Eliminate dependents again if necessary
    if(y.size()!=ny){
      eliminateDependent();
    }
  }
}

Variable& SymbolicOCP::readVariable(const XMLNode& node){
  // Qualified name
  string qn = qualifiedName(node);
  
  // Find and return the variable
  return variable(qn);
}

SX SymbolicOCP::readExpr(const XMLNode& node){
  const string& fullname = node.getName();
  if (fullname.find("exp:")== string::npos) {
    casadi_error("SymbolicOCP::readExpr: unknown - expression is supposed to start with 'exp:' , got " << fullname);
  }
  
  // Chop the 'exp:'
  string name = fullname.substr(4);

  // The switch below is alphabetical, and can be thus made more efficient, for example by using a switch statement of the first three letters, if it would ever become a bottleneck
  if(name.compare("Add")==0){
    return readExpr(node[0]) + readExpr(node[1]);
  } else if(name.compare("Acos")==0){
    return acos(readExpr(node[0]));
  } else if(name.compare("Asin")==0){
    return asin(readExpr(node[0]));
  } else if(name.compare("Atan")==0){
    return atan(readExpr(node[0]));
  } else if(name.compare("Cos")==0){
    return cos(readExpr(node[0]));
  } else if(name.compare("Der")==0){
    Variable v = readVariable(node[0]);
    v.setDifferential(true);
    return v.der();
  } else if(name.compare("Div")==0){
    return readExpr(node[0]) / readExpr(node[1]);
  } else if(name.compare("Exp")==0){
    return exp(readExpr(node[0]));
  } else if(name.compare("Identifier")==0){
    return readVariable(node).var();
  } else if(name.compare("IntegerLiteral")==0){
    return int(node.getText());
  } else if(name.compare("Instant")==0){
    return double(node.getText());
  } else if(name.compare("Log")==0){
    return log(readExpr(node[0]));
  } else if(name.compare("LogLt")==0){ // Logical less than
    return readExpr(node[0]) < readExpr(node[1]);
  } else if(name.compare("LogGt")==0){ // Logical less than
    return readExpr(node[0]) > readExpr(node[1]);
  } else if(name.compare("Mul")==0){ // Multiplication
    return readExpr(node[0]) * readExpr(node[1]);
  } else if(name.compare("Neg")==0){
    return -readExpr(node[0]);
  } else if(name.compare("NoEvent")==0) {
    // NOTE: This is a workaround, we assume that whenever NoEvent occurs, what is meant is a switch
    int n = node.size();
    
    // Default-expression
    SX ex = readExpr(node[n-1]);
    
    // Evaluate ifs
    for(int i=n-3; i>=0; i -= 2) ex = if_else(readExpr(node[i]),readExpr(node[i+1]),ex);
    
    return ex;
  } else if(name.compare("Pow")==0){
    return pow(readExpr(node[0]),readExpr(node[1]));
  } else if(name.compare("RealLiteral")==0){
    return double(node.getText());
  } else if(name.compare("Sin")==0){
    return sin(readExpr(node[0]));
  } else if(name.compare("Sqrt")==0){
    return sqrt(readExpr(node[0]));
  } else if(name.compare("StringLiteral")==0){
    throw CasadiException(string(node.getText()));
  } else if(name.compare("Sub")==0){
    return readExpr(node[0]) - readExpr(node[1]);
  } else if(name.compare("Tan")==0){
    return tan(readExpr(node[0]));
  } else if(name.compare("Time")==0){
    return t.toScalar();
  } else if(name.compare("TimedVariable")==0){
    return readVariable(node[0]).atTime(double(node[1].getText()),true);
  }

  // throw error if reached this point
  throw CasadiException(string("SymbolicOCP::readExpr: unknown node: ") + name);
  
}

void SymbolicOCP::repr(std::ostream &stream) const{
  stream << "Flat OCP";
}

void SymbolicOCP::print(ostream &stream) const{
  stream << "Dimensions: "; 
  stream << "s = " << x.size() << ", ";
  stream << "#xd = " << xd.size() << ", ";
  stream << "#z = " << xa.size() << ", ";
  stream << "#q = " << xq.size() << ", ";
  stream << "#y = " << y.size() << ", ";
  stream << "#p = " << p.size() << ", ";
  stream << "#u = " << u.size() << ", ";
  stream << endl << endl;

  // Variables in the class hierarchy
  stream << "Variables" << endl;

  // Print the variables
  stream << "{" << endl;
  stream << "  t = " << t << endl;
  stream << "  s =  " << x << endl;
  stream << "  xd = " << xd << endl;
  stream << "  z =  " << xa << endl;
  stream << "  q =  " << xq << endl;
  stream << "  y =  " << y << endl;
  stream << "  p =  " << p << endl;
  stream << "  pi =  " << pi << endl;
  stream << "  pd =  " << pd << endl;
  stream << "  p_free =  " << p_free << endl;
  stream << "  ci =  " << ci << endl;
  stream << "  cd =  " << cd << endl;
  stream << "  u =  " << u << endl;
  stream << "}" << endl;
  
  // Print the differential-algebraic equation
  stream << "Implicit dynamic equations" << endl;
  for(vector<SX>::const_iterator it=dae.begin(); it!=dae.end(); it++){
    stream << "0 == "<< *it << endl;
  }
  stream << endl;
  
  stream << "Explicit differential equations" << endl;
  for(int k=0; k<xd.size(); ++k){
    stream << xd.at(k).der() << " == " << ode.at(k) << endl;
  }
  stream << endl;

  stream << "Algebraic equations" << endl;
  for(int k=0; k<xa.size(); ++k){
    stream << xa.at(k) << " == " << alg.at(k) << endl;
  }
  stream << endl;
  
  stream << "Quadrature equations" << endl;
  for(int k=0; k<xq.size(); ++k){
    stream << xq.at(k).der() << " == " << quad.at(k) << endl;
  }
  stream << endl;

  stream << "Initial equations" << endl;
  for(vector<SX>::const_iterator it=initial.begin(); it!=initial.end(); it++){
    stream << "0 == " << *it << endl;
  }
  stream << endl;

  // Dependent equations
  stream << "Dependent equations" << endl;
  for(int i=0; i<y.size(); ++i)
    stream << y.at(i) << " == " << dep.at(i) << endl;
  stream << endl;

  // Mayer terms
  stream << "Mayer objective terms" << endl;
  for(int i=0; i<mterm.size(); ++i)
    stream << mterm.at(i) << endl;
  stream << endl;
  
  // Lagrange terms
  stream << "Lagrange objective terms" << endl;
  for(int i=0; i<lterm.size(); ++i)
    stream << lterm.at(i) << endl;
  stream << endl;
  
  // Constraint functions
  stream << "Constraint functions" << endl;
  for(int i=0; i<path.size(); ++i)
    stream << path_min.at(i) << " <= " << path.at(i) << " <= " << path_max.at(i) << endl;
  stream << endl;
  
  // Constraint functions
  stream << "Time horizon" << endl;
  stream << "t0 = " << t0 << endl;
  stream << "tf = " << tf << endl;
  
}

void SymbolicOCP::eliminateInterdependencies(){
  substituteInPlace(var(y),dep,false);
  
  // Make sure that the dependent variables have been properly eliminated from the dependent expressions
  casadi_assert(!dependsOn(dep,var(y)));
}

vector<SXMatrix> SymbolicOCP::substituteDependents(const vector<SXMatrix>& x) const{
  return substitute(x,var(y),dep);
}

void SymbolicOCP::eliminateDependent(bool eliminate_dependents_with_bounds){
  // All the functions to be replaced
  vector<SXMatrix> fcn(8);
  fcn[0] = dae;
  fcn[1] = ode;
  fcn[2] = alg;
  fcn[3] = quad;
  fcn[4] = initial;
  fcn[5] = path;
  fcn[6] = mterm;
  fcn[7] = lterm;
  
  // Replace all at once
  vector<SXMatrix> fcn_new = substituteDependents(fcn);
  
  // Save the new expressions
  dae = fcn_new[0];
  ode = fcn_new[1];
  alg = fcn_new[2];
  quad = fcn_new[3];
  initial = fcn_new[4];
  path    = fcn_new[5];
  mterm   = fcn_new[6];
  lterm   = fcn_new[7];
}

void SymbolicOCP::eliminateLagrangeTerms(){
  // Index for the names
  int ind = 0;
  // For every integral term in the objective function
  for(vector<SX>::iterator it=lterm.begin(); it!=lterm.end(); ++it){
    
    // Give a name to the quadrature state
    stringstream q_name;
    q_name << "q_" << ind++;
    
    // Create a new quadrature state
    Variable qv(q_name.str());
    qv.setVariability(CONTINUOUS);
    qv.setCausality(INTERNAL);
    qv.setStart(0.0);
    if(tf==tf) qv.setNominal(tf); // if not not-a-number
  
    // Add to the list of variables
    addVariable(q_name.str(),qv);
    
    // Add to the quadrature states
    xq.push_back(qv);

    // Add the Lagrange term to the list of quadratures
    quad.append(*it);
    
    // Add to the list of Mayer terms
    mterm.append(qv.var());
  }
  
  // Remove the Lagrange terms
  lterm.clear();
}

void SymbolicOCP::eliminateQuadratureStates(){
  
  // Move all the quadratures to the list of differential states
  xd.insert(xd.end(),xq.begin(),xq.end());
  xq.clear();
  
  // Move the equations to the list of ODEs
  ode.append(quad);
  quad.clear();
}

void SymbolicOCP::sortType(bool sort_by_variable_category){
  // Clear variables
  x.clear();
  xd.clear();
  xa.clear();
  u.clear();
  cd.clear();
  ci.clear();
  pd.clear();
  pi.clear();
  p_free.clear();
  
  // To be removed
  p.clear();
  
  if(sort_by_variable_category){

    // Loop over variables
    for(map<string,Variable>::iterator it=varmap_.begin(); it!=varmap_.end(); ++it){
      // Get the variable
      Variable& v = it->second;
      
      // Sort by category
      switch(v.getCategory()){
	case CAT_DERIVATIVE:
	  // Skip derivatives
	  break;
	case CAT_STATE:
	  x.push_back(v);
	  break;
	case CAT_DEPENDENT_CONSTANT:
	  cd.push_back(v);
	  break;
	case CAT_INDEPENDENT_CONSTANT:
	  ci.push_back(v);
	  break;
	case CAT_DEPENDENT_PARAMETER:
	  pd.push_back(v);
	  break;
	case CAT_INDEPENDENT_PARAMETER:
	  if(v.getFree()){
	    p_free.push_back(v);
	  } else {
	    pi.push_back(v);
	  }
	  break;
	case CAT_ALGEBRAIC:
	  if(v.getCausality() == INTERNAL){
	    x.push_back(v);
	  } else if(v.getCausality() == INPUT){
	    u.push_back(v);
	  }
	  break;
	default:
	  casadi_assert_message(0,"Unknown category");
      }
    }
   
   // Old implementation
  } else {
  
    // Mark all dependent variables
    for(vector<Variable>::iterator it=y.begin(); it!=y.end(); ++it){
      it->var().setTemp(1);
    }
    
    // Loop over variables
    for(map<string,Variable>::iterator it=varmap_.begin(); it!=varmap_.end(); ++it){
      // Get the variable
      Variable& v = it->second;
      
      // If not dependent
      if(v.var().getTemp()!=1){
	// Try to determine the type
	if(v.getVariability() == PARAMETER){
	  if(v.getFree()){
	    p.push_back(v); 
	  } else {
	    casadi_assert(0);
	  }
	} else if(v.getVariability() == CONTINUOUS) {
	  if(v.getCausality() == INTERNAL){
	    x.push_back(v);
	  } else if(v.getCausality() == INPUT){
	    u.push_back(v);
	  }
	} else if(v.getVariability() == CONSTANT){
	  y.push_back(v);
	  dep.append(v.getNominal());
	}
      }
    }

    // Unmark all dependent variables
    for(vector<Variable>::iterator it=y.begin(); it!=y.end(); ++it){
      it->var().setTemp(0);
    }
  }
}

void SymbolicOCP::scaleVariables(){
  cout << "Scaling variables ..." << endl;
  double time1 = clock();
  
  // Variables
  Matrix<SX> _x = var(x);
  Matrix<SX> _xdot = der(x);
  Matrix<SX> _xd = var(xd);
  Matrix<SX> _xa = var(xa);
  Matrix<SX> _p = var(p);
  Matrix<SX> _u = var(u);
  
  // Collect all the variables
  Matrix<SX> v;
  v.append(t);
  v.append(_x);
  v.append(_xdot);
  v.append(_xd);
  v.append(_xa);
  v.append(_p);
  v.append(_u);
  
  // Nominal values
  Matrix<SX> t_n = 1.;
  Matrix<SX> x_n = getNominal(x);
  Matrix<SX> xd_n = getNominal(xd);
  Matrix<SX> xa_n = getNominal(xa);
  Matrix<SX> p_n = getNominal(p);
  Matrix<SX> u_n = getNominal(u);
  
  // Get all the old variables in expressed in the nominal ones
  Matrix<SX> v_old;
  v_old.append(t*t_n);
  v_old.append(_x*x_n);
  v_old.append(_xdot*x_n);
  v_old.append(_xd*xd_n);
  v_old.append(_xa*xa_n);
  v_old.append(_p*p_n);
  v_old.append(_u*u_n);
  
  // Temporary variable
  Matrix<SX> temp;

  // Substitute equations
  dae = substitute(dae,v,v_old);
  ode = substitute(ode,v,v_old);
  alg = substitute(alg,v,v_old);
  quad = substitute(quad,v,v_old);
  dep = substitute(dep,v,v_old);
  initial = substitute(initial,v,v_old);
  path    = substitute(path,v,v_old);
  mterm   = substitute(mterm,v,v_old);
  lterm   = substitute(lterm,v,v_old);
  
  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... variable scaling complete after " << dt << " seconds." << endl;
}
    
void SymbolicOCP::scaleEquations(){
  
  // Quick return if no implicit equations
  if(dae.empty())
    return;

  cout << "Scaling equations ..." << endl;
  double time1 = clock();

  // Variables
  enum Variables{T,X,XDOT,Z,P,U,NUM_VAR};
  vector<Matrix<SX> > v(NUM_VAR); // all variables
  v[T] = t;
  v[X] = var(xd);
  v[XDOT] = der(xd);
  v[Z] = var(xa);
  v[P] = var(p);
  v[U] = var(u);

  // Create the jacobian of the implicit equations with respect to [x,z,p,u] 
  Matrix<SX> xz;
  xz.append(v[X]);
  xz.append(v[Z]);
  xz.append(v[P]);
  xz.append(v[U]);
  SXFunction fcn = SXFunction(xz,dae);
  SXFunction J(v,fcn.jac());

  // Evaluate the Jacobian in the starting point
  J.init();
  J.setInput(0.0,T);
  J.setInput(getStart(xd,true),X);
  J.input(XDOT).setAll(0.0);
  J.setInput(getStart(xa,true),Z);
  J.setInput(getStart(p,true),P);
  J.setInput(getStart(u,true),U);
  J.evaluate();
  
  // Get the maximum of every row
  Matrix<double> &J0 = J.output();
  vector<double> scale(J0.size1(),0.0); // scaling factors
  for(int i=0; i<J0.size1(); ++i){
    // Loop over non-zero entries of the row
    for(int el=J0.rowind(i); el<J0.rowind(i+1); ++el){
      // Column
      //int j=J0.col(el);
      
      // The scaling factor is the maximum norm, ignoring not-a-number entries
      if(!isnan(J0.at(el))){
        scale[i] = max(scale[i],fabs(J0.at(el)));
      }
    }
    
    // Make sure 
    if(scale[i]==0){
      cout << "Warning: Could not generate a scaling factor for equation " << i << "(0 == " << dae.at(i) << "), selecting 1." << endl;
      scale[i]=1.;
    }
  }
  
  // Scale the equations
  for(int i=0; i<dae.size(); ++i){
    dae[i] /= scale[i];
  }
  
  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... equation scaling complete after " << dt << " seconds." << endl;
}

void SymbolicOCP::sortDAE(){
  // Get the sparsity of the Jacobian df/fx + tau*df/xdot
  SXMatrix dae_only_x = substitute(dae,der(x),ssym("tau")*var(x));
  SXFunction f(var(x),dae_only_x);
  f.init();
  CRSSparsity sp = f.jacSparsity();
  
  // BLT transformation
  vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
  int nb = sp.dulmageMendelsohn(rowperm,colperm,rowblock,colblock,coarse_rowblock,coarse_colblock);

  // Permute equations
  vector<SX> dae_new(dae.size());
  for(int i=0; i<dae.size(); ++i){
    dae_new[i] = dae.at(rowperm[i]);
  }
  dae_new.swap(dae.data());
  
  // Permute variables
  vector<Variable> x_new(x.size());
  for(int i=0; i<x.size(); ++i){
    x_new[i]= x[colperm[i]];
  }
  x_new.swap(x);
}

void SymbolicOCP::makeExplicit(){
  
  // Quick return if there are no implicitly defined states
  if(x.empty()) return;
  
  // Write the DAE as a function of the highest unknown derivatives (algebraic states and state derivatives)
  SXFunction f(highest(x),dae);
  f.init();

  // Get the sparsity of the Jacobian which can be used to determine which variable can be calculated from which other
  CRSSparsity sp = f.jacSparsity();

  // BLT transformation
  vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
  int nb = sp.dulmageMendelsohn(rowperm,colperm,rowblock,colblock,coarse_rowblock,coarse_colblock);

  // Permute equations
  vector<SX> dae_new(dae.size());
  for(int i=0; i<dae.size(); ++i){
    dae_new[i] = dae.at(rowperm[i]);
  }
  dae_new.swap(dae.data());
  dae_new.clear();
  
  // Permute variables
  vector<Variable> x_new(x.size());
  for(int i=0; i<x.size(); ++i){
    x_new[i]= x[colperm[i]];
  }
  x_new.swap(x);
  x_new.clear();

  // Rewrite the sorted DAE as a function of the highest unknown derivatives
  f = SXFunction(highest(x),dae);
  f.init();

  // Get the Jacobian
  SXMatrix J = f.jac();
  
  // Block variables and equations
  vector<Variable> xb, xdb, xab;
  vector<SX> fb;

  // Variables where we have found an explicit expression
  vector<Variable> x_exp;
  
  // Explicit equations
  vector<SX> f_exp;
  
  // Loop over blocks
  for(int b=0; b<nb; ++b){
    
    // Block size
    int bs = rowblock[b+1] - rowblock[b];
    
    // Get local variables
    xb.clear();
    xdb.clear();
    xab.clear();
    for(int i=colblock[b]; i<colblock[b+1]; ++i){
      xb.push_back(x[i]);
      if(x[i].isDifferential()){
        xdb.push_back(x[i]);
      } else {
        xab.push_back(x[i]);
      }
    }

    // Get local equations
    fb.clear();
    for(int i=rowblock[b]; i<rowblock[b+1]; ++i)
      fb.push_back(dae.at(i));

    // Get local Jacobian
    SXMatrix Jb = J(range(rowblock[b],rowblock[b+1]),range(colblock[b],colblock[b+1]));

    // If Jb depends on xb, then we cannot solve for it explicitly
    if(dependsOn(Jb,highest(xb))){
      
      // If the block only contains algebraic states ...
      if(xdb.empty()){
        // ... we can simply add the equations to the list of algebraic equations ...
        alg.append(fb);
        
        // ... and the variables accordingly
        xa.insert(xa.end(),xab.begin(),xab.end());
      } else { // The block contains differential states
        casadi_error("Cannot find an explicit expression for variable(s) " << xdb);
      }
    } else { // The variables that we wish to determine enter linearly
      
      // Divide fb into a part which depends on vb and a part which doesn't according to "fb == mul(Jb,vb) + fb_res"
      SXMatrix fb_res = substitute(fb,highest(xb),SXMatrix(xb.size(),1,0)).data();
      SXMatrix fb_exp;
      
      // Solve for vb
      if (bs <= 3){
        // Calculate inverse and multiply for very small matrices
        fb_exp = mul(inv(Jb),-fb_res);
      } else {
        // QR factorization
        fb_exp = solve(Jb,-fb_res);
      }

      // Add to explicitly determined equations and variables
      x_exp.insert(x_exp.end(),xb.begin(),xb.end());
      f_exp.insert(f_exp.end(),fb_exp.data().begin(),fb_exp.data().end());
    }
  }
  
  // Eliminate inter-dependencies in fb_exp
  SXMatrix f_expmat = f_exp;
  substituteInPlace(highest(x_exp),f_expmat,false);
  f_exp = f_expmat.data();

  // New dependent variables and binding equations
  vector<Variable> y_new;
  vector<SX> dep_new;
  
  // Split the dependent variables from the state derivative expressions
  for(int k=0; k<x_exp.size(); ++k){
    
    // Check if differential state
    if(x_exp[k].isDifferential()){
      // Add to the ODE
      xd.push_back(x_exp[k]);
      ode.append(f_exp[k]);
    } else {
      // Add to the list of the dependent variables
      y_new.push_back(x_exp[k]);
      dep_new.push_back(f_exp[k]);
    }
  }
  
  // Add to the beginning of the dependent variables (since the other dependent variable might depend on them)
  y.insert(y.begin(),y_new.begin(),y_new.end());
  dep = vertcat(SXMatrix(dep_new),dep);
  
  // Remove the eliminated variables and equations
  x.clear();
  dae.clear();
}

vector<Variable> SymbolicOCP::x_all() const{
  vector<Variable> ret;
  ret.insert(ret.end(),x.begin(),x.end());
  ret.insert(ret.end(),xd.begin(),xd.end());
  ret.insert(ret.end(),xa.begin(),xa.end());
  return ret;
}

vector<SXMatrix> SymbolicOCP::daeArg() const{
  // All states
  vector<Variable> _x = x_all();
  
  // Return value
  vector<SXMatrix> ret(DAE_NUM_IN);
  ret[DAE_T] = t;
  ret[DAE_Y] = var(_x);
  ret[DAE_YDOT] = der(_x);
  ret[DAE_P] = vertcat<SX>(var(p),var(u));
  return ret;
}

void SymbolicOCP::makeAlgebraic(const std::string& name){
  makeAlgebraic(variable(name));
}

void SymbolicOCP::makeAlgebraic(const Variable& v){
  // Find variable among the explicit variables
  for(int k=0; k<xd.size(); ++k){
    if(xd[k].get()==v.get()){
      
      // Add to list of algebraic variables and to the list of algebraic equations
      xa.push_back(v);
      alg.append(ode.at(k));
      
      // Remove from list of differential variables and the list of differential equations
      xd.erase(xd.begin()+k);
      vector<SX> ode_ = ode.data();
      ode_.erase(ode_.begin()+k);
      ode = ode_;

      // Successfull return
      return;
    }
  }
  
  // Find the variable among the implicit variables
  for(int k=0; k<x.size(); ++k){
    if(x[k].get()==v.get()){
      
      // Substitute the state derivative with zero
      dae = substitute(dae,x[k].der(),0.0);

      // Remove the highest state derivative expression from the variable
      x[k].setDifferential(false);

      // Successfull return
      return;
    }
  }
  
  // Error if this point reached
  throw CasadiException("v not a differential state");
}

Variable& SymbolicOCP::variable(const std::string& name){
  // Find the variable
  map<string,Variable>::iterator it = varmap_.find(name);
  if(it==varmap_.end()){
    casadi_error("No such variable: \"" << name << "\".");
  }
  
  // Return the variable
  return it->second;
}

void SymbolicOCP::addVariable(const std::string& name, const Variable& var){
  // Try to find the name
  map<string,Variable>::iterator it = varmap_.find(name);
  if(it!=varmap_.end()){
    stringstream ss;
    ss << "Variable \"" << name << "\" has already been added.";
    throw CasadiException(ss.str());
  }
  
  varmap_[name] = var;
}

std::string SymbolicOCP::qualifiedName(const XMLNode& nn){
  // Stringstream to assemble name
  stringstream qn;
  
  for(int i=0; i<nn.size(); ++i){
    // Add a dot
    if(i!=0) qn << ".";
    
    // Get the name part
    string namepart = nn[i].attribute("name");
    qn << namepart;

    // Get the index, if any
    if(nn[i].size()>0){
      int ind = nn[i]["exp:ArraySubscripts"]["exp:IndexExpression"]["exp:IntegerLiteral"].getText();
      qn << "[" << ind << "]";
    }
  }
  
  // Return the name
  return qn.str();
}

void SymbolicOCP::generateMuscodDatFile(const std::string& filename, const Dictionary& mc2_ops) const{
  // Make sure that the OCP is in semi-explicit form
  casadi_assert_message(x.empty(), "The DAE must be in semi-explicit form");
  
  // Print
  cout << "Generating: " << filename << endl;
  
  // Create the datfile
  ofstream datfile;
  datfile.open(filename.c_str());
  datfile.precision(numeric_limits<double>::digits10+2);
  datfile << scientific; // This is really only to force a decimal dot, would be better if it can be avoided
  
  // Print header
  datfile << "* This function was automatically generated by CasADi" << endl;
  datfile << endl;

  // User set options
  for(Dictionary::const_iterator it=mc2_ops.begin(); it!=mc2_ops.end(); ++it){
    // Print the name
    datfile << it->first << endl;
    
    // Get the value
    const GenericType& val = it->second;
    
    // Print value
    if(val.isInt()){
      datfile << int(val) << endl;
    } else if(val.isDouble()){
      datfile << double(val) << endl;
    } else if(val.isString()){
      datfile << string(val) << endl;
    } else if(val.isIntVector()){
      vector<int> valv = val;
      for(int k=0; k<valv.size(); ++k){
        datfile << k << ": " << valv[k] << endl;
      }
    } else if(val.isDoubleVector()){
      vector<double> valv = val;
      for(int k=0; k<valv.size(); ++k){
        datfile << k << ": " << valv[k] << endl;
      }
    } else if(val.isStringVector()){
      vector<string> valv = val;
      for(int k=0; k<valv.size(); ++k){
        datfile << k << ": " << valv[k] << endl;
      }
    }
    datfile << endl;
  }
  
  // Get the stage duration
  double h = tf-t0;

  // Is the stage duration fixed?
  bool h_fix = !t0_free && !tf_free;
  
  // Get bounds on the stage duration
  double h_min=h, h_max=h;
  
  // Set to some dummy variables if stage duration not fixed
  if(!h_fix){
    casadi_warning("h_min and h_max being set to dummy variables!");
    h_min = 0;
    h_max = numeric_limits<double>::infinity();
  }
  
  datfile << "* model stage duration start values, scale factors, and bounds" << endl;
  datfile << "h" << endl;
  datfile << "0: " << h << endl;
  datfile << endl;

  datfile << "h_sca" << endl;
  datfile << "0: " << h << endl;
  datfile << endl;
  
  datfile << "h_min" << endl;
  datfile << "0: " << h_min << endl;
  datfile << endl;
  
  datfile << "h_max" << endl;
  datfile << "0: " << h_max << endl;
  datfile << endl;
  
  datfile << "h_fix" << endl;
  datfile << "0: " << h_fix << endl;
  datfile << endl;
    
  // Parameter properties
  if(!p.empty()){
    datfile << "*  global model parameter start values, scale factors, and bounds" << endl;
    datfile << "p" << endl;
    for(int k=0; k<p.size(); ++k){
      datfile << k << ": " << p[k].getStart() << endl;
    }
    datfile << endl;
    
    datfile << "p_sca" << endl;
    for(int k=0; k<p.size(); ++k){
      datfile << k << ": " << p[k].getNominal() << endl;
    }
    datfile << endl;
    
    datfile << "p_min" << endl;
    for(int k=0; k<p.size(); ++k){
      datfile << k << ": " << p[k].getMin() << endl;
    }
    datfile << endl;
    
    datfile << "p_max" << endl;
    for(int k=0; k<p.size(); ++k){
      datfile << k << ": " << p[k].getMax() << endl;
    }
    datfile << endl;
    
    datfile << "p_fix" << endl;
    for(int k=0; k<p.size(); ++k){
      datfile << k << ": " << (p[k].getMin()==p[k].getMax()) << endl;
    }
    datfile << endl;
    
    datfile << "p_name" << endl;
    for(int k=0; k<p.size(); ++k){
      datfile << k << ": " << p[k].getName() << endl;
    }
    datfile << endl;
    
    datfile << "p_unit" << endl;
    for(int k=0; k<p.size(); ++k){
      datfile << k << ": " << p[k].getUnit() << endl;
    }
    datfile << endl;
  }

  // Differential state properties
  if(!xd.empty()){
    datfile << "*  differential state start values, scale factors, and bounds" << endl;
    datfile << "sd(*,*)" << endl;
    for(int k=0; k<xd.size(); ++k){
      datfile << k << ": " << xd[k].getStart() << endl;
    }
    datfile << endl;
    
    datfile << "sd_sca(*,*)" << endl;
    for(int k=0; k<xd.size(); ++k){
      datfile << k << ": " << xd[k].getNominal() << endl;
    }
    datfile << endl;
    
    datfile << "sd_min(*,*)" << endl;
    for(int k=0; k<xd.size(); ++k){
      datfile << k << ": " << xd[k].getMin() << endl;
    }
    datfile << endl;
    
    datfile << "sd_max(*,*)" << endl;
    for(int k=0; k<xd.size(); ++k){
      datfile << k << ": " << xd[k].getMax() << endl;
    }
    datfile << endl;
    
    datfile << "sd_fix(*,*)" << endl;
    for(int k=0; k<xd.size(); ++k){
      datfile << k << ": " << (xd[k].getMin()==xd[k].getMax()) << endl;
    }
    datfile << endl;

    datfile << "xd_name" << endl;
    for(int k=0; k<xd.size(); ++k){
      datfile << k << ": " << xd[k].getName() << endl;
    }
    datfile << endl;

    datfile << "xd_unit" << endl;
    for(int k=0; k<xd.size(); ++k){
      datfile << k << ": " << xd[k].getUnit() << endl;
    }
    datfile << endl;
  }
  
  // Algebraic state properties
  if(!xa.empty()){
    datfile << "*  algebraic state start values, scale factors, and bounds" << endl;
    datfile << "sa(*,*)" << endl;
    for(int k=0; k<xa.size(); ++k){
      datfile << k << ": " << xa[k].getStart() << endl;
    }
    datfile << endl;
    
    datfile << "sa_sca(*,*)" << endl;
    for(int k=0; k<xa.size(); ++k){
      datfile << k << ": " << xa[k].getNominal() << endl;
    }
    datfile << endl;
    
    datfile << "sa_min(*,*)" << endl;
    for(int k=0; k<xa.size(); ++k){
      datfile << k << ": " << xa[k].getMin() << endl;
    }
    datfile << endl;
    
    datfile << "sa_max(*,*)" << endl;
    for(int k=0; k<xa.size(); ++k){
      datfile << k << ": " << xa[k].getMax() << endl;
    }
    datfile << endl;
    
    datfile << "sa_fix(*,*)" << endl;
    for(int k=0; k<xa.size(); ++k){
      datfile << k << ": " << (xa[k].getMin()==xa[k].getMax()) << endl;
    }
    datfile << endl;

    datfile << "xa_name" << endl;
    for(int k=0; k<xa.size(); ++k){
      datfile << k << ": " << xa[k].getName() << endl;
    }
    datfile << endl;
    
    datfile << "xa_unit" << endl;
    for(int k=0; k<xa.size(); ++k){
      datfile << k << ": " << xa[k].getUnit() << endl;
    }
    datfile << endl;
  }
  
  // Control properties
  if(!u.empty()){
    datfile << "* control start values, scale factors, and bounds" << endl;
    datfile << "u(*,*)" << endl;
    for(int k=0; k<u.size(); ++k){
      datfile << k << ": " << u[k].getStart() << endl;
    }
    datfile << endl;
    
    datfile << "u_sca(*,*)" << endl;
    for(int k=0; k<u.size(); ++k){
      datfile << k << ": " << u[k].getNominal() << endl;
    }
    datfile << endl;
    
    datfile << "u_min(*,*)" << endl;
    for(int k=0; k<u.size(); ++k){
      datfile << k << ": " << u[k].getMin() << endl;
    }
    datfile << endl;
    
    datfile << "u_max(*,*)" << endl;
    for(int k=0; k<u.size(); ++k){
      datfile << k << ": " << u[k].getMax() << endl;
    }
    datfile << endl;
    
    datfile << "u_fix(*,*)" << endl;
    for(int k=0; k<u.size(); ++k){
      datfile << k << ": " << (u[k].getMin()==u[k].getMax()) << endl;
    }
    datfile << endl;
    
    datfile << "u_name" << endl;
    for(int k=0; k<u.size(); ++k){
      datfile << k << ": " << u[k].getName() << endl;
    }
    datfile << endl;
    
    datfile << "u_unit" << endl;
    for(int k=0; k<u.size(); ++k){
      datfile << k << ": " << u[k].getUnit() << endl;
    }
    datfile << endl;
  }

  // Close the datfile
  datfile.close();
}

} // namespace CasADi

