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

#include "flat_ocp_internal.hpp"
#include "variable_internal.hpp"
#include <map>
#include <string>
#include <sstream>
#include <ctime>

#include "casadi/stl_vector_tools.hpp"
#include "external_packages/tinyxml/tinyxml.h"

#include "../casadi/casadi_exception.hpp"
#include "../casadi/pre_c99_support.hpp"
#include "../casadi/stl_vector_tools.hpp"
#include "variable_tools.hpp"
#include "../casadi/matrix/matrix_tools.hpp"
#include "../casadi/sx/sx_tools.hpp"
#include "../casadi/fx/integrator.hpp"

using namespace std;
namespace CasADi{
namespace OptimalControl{

FlatOCPInternal::FlatOCPInternal(const std::string& filename) : filename_(filename){
  addOption("scale_variables",          OT_BOOLEAN,      false, "Scale the variables so that they get unity order of magnitude");
  addOption("eliminate_dependent",      OT_BOOLEAN,      true,  "Eliminate variables that can be expressed as an expression of other variables");
  addOption("scale_equations",          OT_BOOLEAN,      false,  "Scale the implicit equations so that they get unity order of magnitude");
  addOption("sort_equations",           OT_BOOLEAN,      true,  "Sort the dynamic equations");
  addOption("make_explicit",            OT_BOOLEAN,      false, "Make the DAE semi-explicit");
  addOption("eliminate_algebraic",      OT_BOOLEAN,      false, "Completely eliminate algebraic states");
  addOption("verbose",                  OT_BOOLEAN,      true,  "Verbose parsing");
  addOption("eliminate_dependents_with_bounds",  OT_BOOLEAN,      true,  "Verbose parsing");
  
  TiXmlDocument doc;
  bool flag = doc.LoadFile(filename.data());

  if(!flag){
    throw CasadiException("XMLParser::loadFile: Cound not open " + filename);
  }

  // parse
  document_.setName(filename);
  document_.addNode(&doc);

  t_ = SX("t");
  t0_ = numeric_limits<double>::quiet_NaN();
  tf_ = numeric_limits<double>::quiet_NaN();
}

void FlatOCPInternal::init(){

  // Read options
  verbose_ = getOption("verbose");
  bool scale_variables = getOption("scale_variables");
  bool scale_equations = getOption("scale_equations");
  bool eliminate_dependent = getOption("eliminate_dependent");
  bool sort_equations = getOption("sort_equations");
  bool make_explicit = getOption("make_explicit");
  
  // Obtain the symbolic representation of the OCP
  parse();

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
    int ny = y_.size();
    
    // Solve for the highest order derivatives
    makeExplicit();
    
    // Eliminate dependents again if necessary
    if(y_.size()!=ny){
      eliminateDependent();
    }
  }
}

FlatOCPInternal::~FlatOCPInternal(){

}

void FlatOCPInternal::parse(){
  cout << "Parsing XML ..." << endl;
  double time1 = clock();

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

  // Sort the variables according to type
  sortType();
  
  // Make sure that the dimensions are consistent at this point
  casadi_assert(x_.size()==dae_.size());
  casadi_assert(xd_.size()==ode_.size());
  casadi_assert(xa_.size()==alg_.size());
  casadi_assert(q_.size()==quad_.size());
  casadi_assert(y_.size()==dep_.size());
  
  // Return a reference to the created ocp
  double time2 = clock();
  double tparse = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... parsing complete after " << tparse << " seconds" << endl;
}

void FlatOCPInternal::addModelVariables(){

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
      const XMLNode& props = vnode[0];
      if(props.hasAttribute("unit"))         var.setUnit(props.attribute("unit"));
      if(props.hasAttribute("displayUnit"))  var.setDisplayUnit(props.attribute("displayUnit"));
      if(props.hasAttribute("min"))          var.setMin(props.attribute("min"));
      if(props.hasAttribute("max"))          var.setMax(props.attribute("max"));
      if(props.hasAttribute("start"))        var.setStart(props.attribute("start"));
      if(props.hasAttribute("nominal"))      var.setNominal(props.attribute("nominal"));
      if(props.hasAttribute("free"))         var.setFree(string(props.attribute("free")).compare("true") == 0);
      
      // Add to list of variables
      addVariable(qn,var);
    }
  }
}

void FlatOCPInternal::addBindingEquations(){
  if(verbose_){
    cout << "Adding binding equations." << endl;
  }
  
  // Get a reference to the BindingEquations node
  const XMLNode& bindeqs = document_[0]["equ:BindingEquations"];
  
  for(int i=0; i<bindeqs.size(); ++i){
    const XMLNode& beq = bindeqs[i];

    // Get the variable
    Variable var = readVariable(beq[0]);

    // Get the binding equation
    SX bexpr = readExpr(beq[1][0]);
    
    // Add binding equation
    y_.push_back(var);
    dep_.push_back(bexpr);
  }
}

void FlatOCPInternal::addDynamicEquations(){
  // Get a reference to the DynamicEquations node
  const XMLNode& dyneqs = document_[0]["equ:DynamicEquations"];

  // Add equations
  for(int i=0; i<dyneqs.size(); ++i){

    // Get a reference to the variable
    const XMLNode& dnode = dyneqs[i];

    // Add the differential equation
    SX de_new = readExpr(dnode[0]);
    dae_.push_back(de_new);
  }
}

void FlatOCPInternal::addInitialEquations(){
  // Get a reference to the DynamicEquations node
  const XMLNode& initeqs = document_[0]["equ:InitialEquations"];

  // Add equations
  for(int i=0; i<initeqs.size(); ++i){

    // Get a reference to the node
    const XMLNode& inode = initeqs[i];

    // Add the differential equations
    for(int i=0; i<inode.size(); ++i){
      initial_.push_back(readExpr(inode[i]));
    }
  }
}

void FlatOCPInternal::addOptimization(){
  // Get a reference to the DynamicEquations node
  const XMLNode& opts = document_[0]["opt:Optimization"];
  
  // Start time
  cout << "starttime (string) = " << string(opts["opt:IntervalStartTime"]["opt:Value"].getText()) << endl;
  t0_  = opts["opt:IntervalStartTime"]["opt:Value"].getText();
  cout << "starttime (double) = " << t0_ << endl;

  // Terminal time
  cout << "endtime (string) = " << string(opts["opt:IntervalFinalTime"]["opt:Value"].getText()) << endl;
  tf_ = opts["opt:IntervalFinalTime"]["opt:Value"].getText();
  cout << "endtime (double) = " << tf_ << endl;

  for(int i=0; i<opts.size(); ++i){
    
    // Get a reference to the node
    const XMLNode& onode = opts[i];

    // Get the type
    if(onode.checkName("opt:ObjectiveFunction")){ // mayer term
      try{
        addObjectiveFunction(onode);
      } catch(exception& ex){
        cout << "WARNING: addObjectiveFunction" << ex.what() << endl;
      }
    } else if(onode.checkName("opt:IntegrandObjectiveFunction")){
      try{
        addIntegrandObjectiveFunction(onode);
      } catch(exception& ex){
        cout << "WARNING: addIntegrandObjectiveFunction" << ex.what() << endl;
      }
    } else if(onode.checkName("opt:IntervalStartTime")) {
       addIntervalStartTime(onode);
    } else if(onode.checkName("opt:IntervalFinalTime")) {
       addIntervalFinalTime(onode);
    } else if(onode.checkName("opt:TimePoints")) {
//       addTimePoints(onode);
    } else if(onode.checkName("opt:Constraints")) {
      addConstraints(onode);
    } else throw "FlatOCPInternal::addOptimization: Unknown node";
  }
}

void FlatOCPInternal::addObjectiveFunction(const XMLNode& onode){
  // Add components
  for(int i=0; i<onode.size(); ++i){
    const XMLNode& var = onode[i];
    SX v = readExpr(var);
    mterm_.push_back(v);
  }
}

void FlatOCPInternal::addIntegrandObjectiveFunction(const XMLNode& onode){
  for(int i=0; i<onode.size(); ++i){
    const XMLNode& var = onode[i];
    SX v = readExpr(var);
    lterm_.push_back(v);
  }
}

void FlatOCPInternal::addIntervalStartTime(const XMLNode& onode){
  
}

void FlatOCPInternal::addIntervalFinalTime(const XMLNode& onode){
  
}


void FlatOCPInternal::addConstraints(const XMLNode& onode){
  for(int i=0; i<onode.size(); ++i){

    const XMLNode& constr_i = onode[i];
    if(constr_i.checkName("opt:ConstraintLeq")){
      SX ex = readExpr(constr_i[0]);
      SX ub = readExpr(constr_i[1]);
      path_.push_back(ex-ub);
      path_min_.push_back(-numeric_limits<double>::infinity());
      path_max_.push_back(0.);
    } else if(constr_i.checkName("opt:ConstraintGeq")){
      SX ex = readExpr(constr_i[0]);
      SX lb = readExpr(constr_i[1]);
      path_.push_back(ex-lb);
      path_min_.push_back(0.);
      path_max_.push_back(numeric_limits<double>::infinity());
    } else if(constr_i.checkName("opt:ConstraintEq")){
      SX ex = readExpr(constr_i[0]);
      SX eq = readExpr(constr_i[1]);
      path_.push_back(ex-eq);
      path_min_.push_back(0.);
      path_max_.push_back(0.);
    } else {
      cerr << "unknown constraint type" << constr_i.getName() << endl;
      throw "FlatOCPInternal::addConstraints";
    }
  }
}

Variable& FlatOCPInternal::readVariable(const XMLNode& node){
  // Qualified name
  string qn = qualifiedName(node);
  
  // Find and return the variable
  return variable(qn);
}


SX FlatOCPInternal::readExpr(const XMLNode& node){
  const string& fullname = node.getName();
  if (fullname.find("exp:")== string::npos) {
    stringstream ss;
    ss << "FlatOCPInternal::readExpr: unknown - expression is supposed to start with 'exp:' , got " << fullname;
    throw CasadiException(ss.str());
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
    return t_;
  } else if(name.compare("TimedVariable")==0){
    return readVariable(node[0]).atTime(double(node[1].getText()),true);
  }

  // throw error if reached this point
  throw CasadiException(string("FlatOCPInternal::readExpr: unknown node: ") + name);
  
}

void FlatOCPInternal::repr(std::ostream &stream) const{
  stream << "Flat OCP (XML file: \"" << filename_ << "\")";
}

void FlatOCPInternal::print(ostream &stream) const{
  stream << "Dimensions: "; 
  stream << "#s = " << x_.size() << ", ";
  stream << "#xd = " << xd_.size() << ", ";
  stream << "#z = " << xa_.size() << ", ";
  stream << "#q = " << q_.size() << ", ";
  stream << "#y = " << y_.size() << ", ";
  stream << "#p = " << p_.size() << ", ";
  stream << "#u = " << u_.size() << ", ";
  stream << endl << endl;

  // Variables in the class hierarchy
  stream << "Variables" << endl;

  // Print the variables
  stream << "{" << endl;
  stream << "  t = " << t_ << endl;
  stream << "  s =  " << x_ << endl;
  stream << "  xd = " << xd_ << endl;
  stream << "  z =  " << xa_ << endl;
  stream << "  q =  " << q_ << endl;
  stream << "  y =  " << y_ << endl;
  stream << "  p =  " << p_ << endl;
  stream << "  u =  " << u_ << endl;
  stream << "}" << endl;
  
  // Print the differential-algebraic equation
  stream << "Implicit dynamic equations" << endl;
  for(vector<SX>::const_iterator it=dae_.begin(); it!=dae_.end(); it++){
    stream << "0 == "<< *it << endl;
  }
  stream << endl;
  
  stream << "Explicit differential equations" << endl;
  for(int k=0; k<xd_.size(); ++k){
    stream << xd_[k].der() << " == " << ode_[k] << endl;
  }
  stream << endl;

  stream << "Algebraic equations" << endl;
  for(int k=0; k<xa_.size(); ++k){
    stream << xa_[k] << " == " << alg_[k] << endl;
  }
  stream << endl;
  
  stream << "Quadrature equations" << endl;
  for(int k=0; k<q_.size(); ++k){
    stream << q_[k].der() << " == " << quad_[k] << endl;
  }
  stream << endl;

  stream << "Initial equations" << endl;
  for(vector<SX>::const_iterator it=initial_.begin(); it!=initial_.end(); it++){
    stream << "0 == " << *it << endl;
  }
  stream << endl;

  // Dependent equations
  stream << "Dependent equations" << endl;
  for(int i=0; i<y_.size(); ++i)
    stream << y_[i] << " == " << dep_[i] << endl;
  stream << endl;

  // Mayer terms
  stream << "Mayer objective terms" << endl;
  for(int i=0; i<mterm_.size(); ++i)
    stream << mterm_[i] << endl;
  stream << endl;
  
  // Lagrange terms
  stream << "Lagrange objective terms" << endl;
  for(int i=0; i<lterm_.size(); ++i)
    stream << lterm_[i] << endl;
  stream << endl;
  
  // Constraint functions
  stream << "Constraint functions" << endl;
  for(int i=0; i<path_.size(); ++i)
    stream << path_min_[i] << " <= " << path_[i] << " <= " << path_max_[i] << endl;
  stream << endl;
  
  // Constraint functions
  stream << "Time horizon" << endl;
  stream << "t0 = " << t0_ << endl;
  stream << "tf = " << tf_ << endl;
  
}

void FlatOCPInternal::eliminateInterdependencies(){
  bool eliminate_constants = true; // also simplify constant expressions
  dep_ = substituteInPlace(var(y_),dep_,eliminate_constants).data();
}

vector<SXMatrix> FlatOCPInternal::substituteDependents(const vector<SXMatrix>& x) const{
  return substitute(x,var(y_),dep_);
}

void FlatOCPInternal::eliminateDependent(){
  bool eliminate_dependents_with_bounds=getOption("eliminate_dependents_with_bounds");
  if(verbose_)
    cout << "eliminateDependent ..." << endl;
  double time1 = clock();
  
  // All the functions to be replaced
  vector<SXMatrix> fcn(8);
  fcn[0] = dae_;
  fcn[1] = ode_;
  fcn[2] = alg_;
  fcn[3] = quad_;
  fcn[4] = initial_;
  fcn[5] = path_;
  fcn[6] = mterm_;
  fcn[7] = lterm_;
  
  // Replace all at once
  vector<SXMatrix> fcn_new = substituteDependents(fcn);
  
  // Save the new expressions
  dae_= fcn_new[0].data();
  ode_= fcn_new[1].data();
  alg_= fcn_new[2].data();
  quad_= fcn_new[3].data();
  initial_= fcn_new[4].data();
  path_    = fcn_new[5].data();
  mterm_   = fcn_new[6].data();
  lterm_   = fcn_new[7].data();

  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  if(verbose_)
    cout << "... eliminateDependent complete after " << dt << " seconds." << endl;
}

void FlatOCPInternal::sortType(){
  // Clear variables
  x_.clear();
  xd_.clear();
  xa_.clear();
  u_.clear();
  p_.clear();
  
  // Mark all dependent variables
  for(vector<Variable>::iterator it=y_.begin(); it!=y_.end(); ++it){
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
          p_.push_back(v); 
        } else {
          casadi_assert(0);
        }
      } else if(v.getVariability() == CONTINUOUS) {
        if(v.getCausality() == INTERNAL){
          x_.push_back(v);
        } else if(v.getCausality() == INPUT){
          u_.push_back(v);
        }
      } else if(v.getVariability() == CONSTANT){
        y_.push_back(v);
        dep_.push_back(v.getNominal());
      }
    }
  }

  // Unmark all dependent variables
  for(vector<Variable>::iterator it=y_.begin(); it!=y_.end(); ++it){
    it->var().setTemp(0);
  }
}

void FlatOCPInternal::scaleVariables(){
  cout << "Scaling variables ..." << endl;
  double time1 = clock();
  
  // Variables
  Matrix<SX> t = t_;
  Matrix<SX> x = var(x_);
  Matrix<SX> xdot = der(x_);
  Matrix<SX> xd = var(xd_);
  Matrix<SX> xa = var(xa_);
  Matrix<SX> p = var(p_);
  Matrix<SX> u = var(u_);
  
  // Collect all the variables
  Matrix<SX> v;
  append(v,t);
  append(v,x);
  append(v,xdot);
  append(v,xd);
  append(v,xa);
  append(v,p);
  append(v,u);
  
  // Nominal values
  Matrix<SX> t_n = 1.;
  Matrix<SX> x_n = getNominal(x_);
  Matrix<SX> xd_n = getNominal(xd_);
  Matrix<SX> xa_n = getNominal(xa_);
  Matrix<SX> p_n = getNominal(p_);
  Matrix<SX> u_n = getNominal(u_);
  
  // Get all the old variables in expressed in the nominal ones
  Matrix<SX> v_old;
  append(v_old,t*t_n);
  append(v_old,x*x_n);
  append(v_old,xdot*x_n);
  append(v_old,xd*xd_n);
  append(v_old,xa*xa_n);
  append(v_old,p*p_n);
  append(v_old,u*u_n);
  
  // Temporary variable
  Matrix<SX> temp;

  // Substitute equations
  dae_= substitute(dae_,v,v_old).data();
  ode_= substitute(ode_,v,v_old).data();
  alg_= substitute(alg_,v,v_old).data();
  quad_= substitute(quad_,v,v_old).data();
  dep_= substitute(dep_,v,v_old).data();
  initial_= substitute(initial_,v,v_old).data();
  path_    = substitute(path_,v,v_old).data();
  mterm_   = substitute(mterm_,v,v_old).data();
  lterm_   = substitute(lterm_,v,v_old).data();
  
  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... variable scaling complete after " << dt << " seconds." << endl;
}
    
void FlatOCPInternal::scaleEquations(){
  
  // Quick return if no implicit equations
  if(dae_.empty())
    return;

  cout << "Scaling equations ..." << endl;
  double time1 = clock();

  // Variables
  enum Variables{T,X,XDOT,Z,P,U,NUM_VAR};
  vector<Matrix<SX> > v(NUM_VAR); // all variables
  v[T] = t_;
  v[X] = var(xd_);
  v[XDOT] = der(xd_);
  v[Z] = var(xa_);
  v[P] = var(p_);
  v[U] = var(u_);

  // Create the jacobian of the implicit equations with respect to [x,z,p,u] 
  Matrix<SX> xz;
  append(xz,v[X]);
  append(xz,v[Z]);
  append(xz,v[P]);
  append(xz,v[U]);
  SXFunction fcn = SXFunction(xz,dae_);
  SXFunction J(v,fcn.jac());

  // Evaluate the Jacobian in the starting point
  J.init();
  J.setInput(0.0,T);
  J.setInput(getStart(xd_,true),X);
  J.input(XDOT).setAll(0.0);
  J.setInput(getStart(xa_,true),Z);
  J.setInput(getStart(p_,true),P);
  J.setInput(getStart(u_,true),U);
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
      cout << "Warning: Could not generate a scaling factor for equation " << i << "(0 == " << dae_[i] << "), selecting 1." << endl;
      scale[i]=1.;
    }
  }
  
  // Scale the equations
  for(int i=0; i<dae_.size(); ++i){
    dae_[i] /= scale[i];
  }
  
  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... equation scaling complete after " << dt << " seconds." << endl;
}

void FlatOCPInternal::sortDAE(){
  if(verbose_)
    cout << "Sorting DAE" << endl;

  // Get the sparsity of the Jacobian df/fx + tau*df/xdot
  SXMatrix dae_only_x = substitute(dae_,der(x_),symbolic("tau")*var(x_));
  SXFunction f(var(x_),dae_only_x);
  f.init();
  CRSSparsity sp = f.jacSparsity();
  
  // BLT transformation
  vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
  int nb = sp.dulmageMendelsohn(rowperm,colperm,rowblock,colblock,coarse_rowblock,coarse_colblock);

  // Permute equations
  vector<SX> dae_new(dae_.size());
  for(int i=0; i<dae_.size(); ++i){
    dae_new[i] = dae_[rowperm[i]];
  }
  dae_new.swap(dae_);
  
  // Permute variables
  vector<Variable> x_new(x_.size());
  for(int i=0; i<x_.size(); ++i){
    x_new[i]= x_[colperm[i]];
  }
  x_new.swap(x_);
}

void FlatOCPInternal::makeExplicit(){
  
  // Quick return if there are no implicitly defined states
  if(x_.empty()) return;
  
  // Write the DAE as a function of the highest unknown derivatives (algebraic states and state derivatives)
  SXFunction f(highest(x_),dae_);
  f.init();

  // Get the sparsity of the Jacobian which can be used to determine which variable can be calculated from which other
  CRSSparsity sp = f.jacSparsity();

  // BLT transformation
  vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
  int nb = sp.dulmageMendelsohn(rowperm,colperm,rowblock,colblock,coarse_rowblock,coarse_colblock);

  // Permute equations
  vector<SX> dae_new(dae_.size());
  for(int i=0; i<dae_.size(); ++i){
    dae_new[i] = dae_[rowperm[i]];
  }
  dae_new.swap(dae_);
  dae_new.clear();
  
  // Permute variables
  vector<Variable> x_new(x_.size());
  for(int i=0; i<x_.size(); ++i){
    x_new[i]= x_[colperm[i]];
  }
  x_new.swap(x_);
  x_new.clear();

  // Rewrite the sorted DAE as a function of the highest unknown derivatives
  f = SXFunction(highest(x_),dae_);
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
      xb.push_back(x_[i]);
      if(x_[i].isDifferential()){
        xdb.push_back(x_[i]);
      } else {
        xab.push_back(x_[i]);
      }
    }

    // Get local equations
    fb.clear();
    for(int i=rowblock[b]; i<rowblock[b+1]; ++i)
      fb.push_back(dae_[i]);

    // Get local Jacobian
    SXMatrix Jb = J(range(rowblock[b],rowblock[b+1]),range(colblock[b],colblock[b+1]));

    // If Jb depends on xb, then we cannot solve for it explicitly
    if(dependsOn(Jb,highest(xb))){
      
      // If the block only contains algebraic states ...
      if(xdb.empty()){
        // ... we can simply add the equations to the list of algebraic equations ...
        alg_.insert(alg_.end(),fb.begin(),fb.end());
        
        // ... and the variables accordingly
        xa_.insert(xa_.end(),xab.begin(),xab.end());
      } else { // The block contains differential states
        stringstream ss;
        ss << "Cannot find an explicit expression for variable(s) " << xdb;
        throw CasadiException(ss.str());
      }
    } else { // The variables that we wish to determine enter linearly
      
      // Divide fb into a part which depends on vb and a part which doesn't according to "fb == prod(Jb,vb) + fb_res"
      SXMatrix fb_res = substitute(fb,highest(xb),SXMatrix(xb.size(),1,0)).data();
      SXMatrix fb_exp;
      
      // Solve for vb
      if (bs <= 3){
        // Calculate inverse and multiply for very small matrices
        fb_exp = prod(inv(Jb),-fb_res);
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
  bool eliminate_constants = true; // also simplify constant expressions
  f_exp = substituteInPlace(highest(x_exp),f_exp,eliminate_constants).data();

  // New dependent variables and binding equations
  vector<Variable> y_new;
  vector<SX> dep_new;
  
  // Split the dependent variables from the state derivative expressions
  for(int k=0; k<x_exp.size(); ++k){
    
    // Check if differential state
    if(x_exp[k].isDifferential()){
      // Add to the ODE
      xd_.push_back(x_exp[k]);
      ode_.push_back(f_exp[k]);
    } else {
      // Add to the list of the dependent variables
      y_new.push_back(x_exp[k]);
      dep_new.push_back(f_exp[k]);
    }
  }
  
  // Add to the beginning of the dependent variables (since the other dependent variable might depend on them)
  y_.insert(y_.begin(),y_new.begin(),y_new.end());
  dep_.insert(dep_.begin(),dep_new.begin(),dep_new.end());
  
  // Remove the eliminated variables and equations
  x_.clear();
  dae_.clear();
}

vector<Variable> FlatOCPInternal::x_all() const{
  vector<Variable> ret;
  ret.insert(ret.end(),x_.begin(),x_.end());
  ret.insert(ret.end(),xd_.begin(),xd_.end());
  ret.insert(ret.end(),xa_.begin(),xa_.end());
  return ret;
}

vector<SXMatrix> FlatOCPInternal::daeArg() const{
  // All states
  vector<Variable> x = x_all();
  
  // Return value
  vector<SXMatrix> ret(DAE_NUM_IN);
  ret[DAE_T] = t_;
  ret[DAE_Y] = var(x);
  ret[DAE_YDOT] = der(x);
  ret[DAE_P] = vertcat<SX>(var(p_),var(u_));
  return ret;
}
    
FX FlatOCPInternal::daeFcn() const{
  
  // Get the DAE arguments
  vector<SXMatrix> dae_in = daeArg();
  
  // Get the DAE right hand side
  SXMatrix dae_out = dae_;
  if(!xd_.empty()){
    
    // Explicit DAE part
    SXMatrix ode_explicit = ode_;
    ode_explicit -= der(xd_);
    
    // Append to DAE
    dae_out = vertcat(dae_out,ode_explicit);
  }
  if(!xa_.empty()){
    // Append algebraic equation
    dae_out = vertcat<SX>(dae_out,alg_);
  }
    
  // Create the DAE right hand side function
  SXFunction daefcn(dae_in,dae_out);
  return daefcn;
}

void FlatOCPInternal::makeAlgebraic(const Variable& v){
  // Find variable among the explicit variables
  for(int k=0; k<xd_.size(); ++k){
    if(xd_[k].get()==v.get()){
      
      // Add to list of algebraic variables and to the list of algebraic equations
      xa_.push_back(v);
      alg_.push_back(ode_[k]);
      
      // Remove from list of differential variables and the list of diffential equations
      xd_.erase(xd_.begin()+k);
      ode_.erase(ode_.begin()+k);

      // Successfull return
      return;
    }
  }
  
  // Find the variable among the implicit variables
  for(int k=0; k<x_.size(); ++k){
    if(x_[k].get()==v.get()){
      
      // Substitute the state derivative with zero
      dae_ = substitute(dae_,x_[k].der(),0.0).data();

      // Remove the highest state derivative expression from the variable
      x_[k].setDifferential(false);

      // Successfull return
      return;
    }
  }
  
  // Error if this point reached
  throw CasadiException("v not a differential state");
}

Variable& FlatOCPInternal::variable(const std::string& name){
  // Find the variable
  map<string,Variable>::iterator it = varmap_.find(name);
  if(it==varmap_.end()){
    stringstream ss;
    ss << "No such variable: \"" << name << "\".";
    throw CasadiException(ss.str());
  }
  
  // Return the variable
  return it->second;
}

void FlatOCPInternal::addVariable(const std::string& name, const Variable& var){
  // Try to find the name
  map<string,Variable>::iterator it = varmap_.find(name);
  if(it!=varmap_.end()){
    stringstream ss;
    ss << "Variable \"" << name << "\" has already been added.";
    throw CasadiException(ss.str());
  }
  
  varmap_[name] = var;
}

std::string FlatOCPInternal::qualifiedName(const XMLNode& nn){
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

void FlatOCPInternal::generateMuscodDatFile(const std::string& filename, const Dictionary& mc2_ops) const{
  // Make sure that the OCP is in semi-explicit form
  casadi_assert_message(x_.empty(), "The DAE must be in semi-explicit form");
  
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
  double h = tf_-t0_;

  // Is the stage duration fixed?
  bool h_fix = !t0_free_ && !tf_free_;
  
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
  if(!p_.empty()){
    datfile << "*  global model parameter start values, scale factors, and bounds" << endl;
    datfile << "p" << endl;
    for(int k=0; k<p_.size(); ++k){
      datfile << k << ": " << p_[k].getStart() << endl;
    }
    datfile << endl;
    
    datfile << "p_sca" << endl;
    for(int k=0; k<p_.size(); ++k){
      datfile << k << ": " << p_[k].getNominal() << endl;
    }
    datfile << endl;
    
    datfile << "p_min" << endl;
    for(int k=0; k<p_.size(); ++k){
      datfile << k << ": " << p_[k].getMin() << endl;
    }
    datfile << endl;
    
    datfile << "p_max" << endl;
    for(int k=0; k<p_.size(); ++k){
      datfile << k << ": " << p_[k].getMax() << endl;
    }
    datfile << endl;
    
    datfile << "p_fix" << endl;
    for(int k=0; k<p_.size(); ++k){
      datfile << k << ": " << (p_[k].getMin()==p_[k].getMax()) << endl;
    }
    datfile << endl;
    
    datfile << "p_name" << endl;
    for(int k=0; k<p_.size(); ++k){
      datfile << k << ": " << p_[k].getName() << endl;
    }
    datfile << endl;
    
    datfile << "p_unit" << endl;
    for(int k=0; k<p_.size(); ++k){
      datfile << k << ": " << p_[k].getUnit() << endl;
    }
    datfile << endl;
  }

  // Differential state properties
  if(!xd_.empty()){
    datfile << "*  differential state start values, scale factors, and bounds" << endl;
    datfile << "sd(*,*)" << endl;
    for(int k=0; k<xd_.size(); ++k){
      datfile << k << ": " << xd_[k].getStart() << endl;
    }
    datfile << endl;
    
    datfile << "sd_sca(*,*)" << endl;
    for(int k=0; k<xd_.size(); ++k){
      datfile << k << ": " << xd_[k].getNominal() << endl;
    }
    datfile << endl;
    
    datfile << "sd_min(*,*)" << endl;
    for(int k=0; k<xd_.size(); ++k){
      datfile << k << ": " << xd_[k].getMin() << endl;
    }
    datfile << endl;
    
    datfile << "sd_max(*,*)" << endl;
    for(int k=0; k<xd_.size(); ++k){
      datfile << k << ": " << xd_[k].getMax() << endl;
    }
    datfile << endl;
    
    datfile << "sd_fix(*,*)" << endl;
    for(int k=0; k<xd_.size(); ++k){
      datfile << k << ": " << (xd_[k].getMin()==xd_[k].getMax()) << endl;
    }
    datfile << endl;

    datfile << "xd_name" << endl;
    for(int k=0; k<xd_.size(); ++k){
      datfile << k << ": " << xd_[k].getName() << endl;
    }
    datfile << endl;

    datfile << "xd_unit" << endl;
    for(int k=0; k<xd_.size(); ++k){
      datfile << k << ": " << xd_[k].getUnit() << endl;
    }
    datfile << endl;
  }
  
  // Algebraic state properties
  if(!xa_.empty()){
    datfile << "*  algebraic state start values, scale factors, and bounds" << endl;
    datfile << "sa(*,*)" << endl;
    for(int k=0; k<xa_.size(); ++k){
      datfile << k << ": " << xa_[k].getStart() << endl;
    }
    datfile << endl;
    
    datfile << "sa_sca(*,*)" << endl;
    for(int k=0; k<xa_.size(); ++k){
      datfile << k << ": " << xa_[k].getNominal() << endl;
    }
    datfile << endl;
    
    datfile << "sa_min(*,*)" << endl;
    for(int k=0; k<xa_.size(); ++k){
      datfile << k << ": " << xa_[k].getMin() << endl;
    }
    datfile << endl;
    
    datfile << "sa_max(*,*)" << endl;
    for(int k=0; k<xa_.size(); ++k){
      datfile << k << ": " << xa_[k].getMax() << endl;
    }
    datfile << endl;
    
    datfile << "sa_fix(*,*)" << endl;
    for(int k=0; k<xa_.size(); ++k){
      datfile << k << ": " << (xa_[k].getMin()==xa_[k].getMax()) << endl;
    }
    datfile << endl;

    datfile << "xa_name" << endl;
    for(int k=0; k<xa_.size(); ++k){
      datfile << k << ": " << xa_[k].getName() << endl;
    }
    datfile << endl;
    
    datfile << "xa_unit" << endl;
    for(int k=0; k<xa_.size(); ++k){
      datfile << k << ": " << xa_[k].getUnit() << endl;
    }
    datfile << endl;
  }
  
  // Control properties
  if(!u_.empty()){
    datfile << "* control start values, scale factors, and bounds" << endl;
    datfile << "u(*,*)" << endl;
    for(int k=0; k<u_.size(); ++k){
      datfile << k << ": " << u_[k].getStart() << endl;
    }
    datfile << endl;
    
    datfile << "u_sca(*,*)" << endl;
    for(int k=0; k<u_.size(); ++k){
      datfile << k << ": " << u_[k].getNominal() << endl;
    }
    datfile << endl;
    
    datfile << "u_min(*,*)" << endl;
    for(int k=0; k<u_.size(); ++k){
      datfile << k << ": " << u_[k].getMin() << endl;
    }
    datfile << endl;
    
    datfile << "u_max(*,*)" << endl;
    for(int k=0; k<u_.size(); ++k){
      datfile << k << ": " << u_[k].getMax() << endl;
    }
    datfile << endl;
    
    datfile << "u_fix(*,*)" << endl;
    for(int k=0; k<u_.size(); ++k){
      datfile << k << ": " << (u_[k].getMin()==u_[k].getMax()) << endl;
    }
    datfile << endl;
    
    datfile << "u_name" << endl;
    for(int k=0; k<u_.size(); ++k){
      datfile << k << ": " << u_[k].getName() << endl;
    }
    datfile << endl;
    
    datfile << "u_unit" << endl;
    for(int k=0; k<u_.size(); ++k){
      datfile << k << ": " << u_[k].getUnit() << endl;
    }
    datfile << endl;
  }

  // Close the datfile
  datfile.close();
}

} // namespace OptimalControl
} // namespace CasADi
