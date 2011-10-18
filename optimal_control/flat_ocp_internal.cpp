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
  addOption("scale_equations",          OT_BOOLEAN,      true,  "Scale the implicit equations so that they get unity order of magnitude");
  addOption("semi_explicit",            OT_BOOLEAN,      false, "Make the DAE semi-explicit");
  addOption("fully_explicit",           OT_BOOLEAN,      false, "Make the DAE fully explicit (not always possible)");

  TiXmlDocument doc;
  bool flag = doc.LoadFile(filename.data());

  if(!flag){
    throw CasadiException("XMLParser::loadFile: Cound not open " + filename);
  }

  // parse
  document_.setName(filename);
  document_.addNode(&doc);

  scaled_variables_ = false;
  scaled_equations_ = false;
  eliminated_dependents_ = false;
  blt_sorted_ = false;
  t_ = SX("t");
  t0_ = numeric_limits<double>::quiet_NaN();
  tf_ = numeric_limits<double>::quiet_NaN();
}

void FlatOCPInternal::init(){

  // Obtain the symbolic representation of the OCP
  parse();

  // Scale the variables
  scaleVariables();

  // Eliminate the dependent variables
  eliminateDependent();

  // Scale the equations
  scaleEquations();
  
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
    map<string,int>::iterator it = varname_.find(qn);
    
    // Add variable, if not already added
    if(it == varname_.end()){
      
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
  // Get a reference to the BindingEquations node
  const XMLNode& bindeqs = document_[0]["equ:BindingEquations"];
  
  for(int i=0; i<bindeqs.size(); ++i){
    const XMLNode& beq = bindeqs[i];

    // Get the variable
    Variable var = readVariable(beq[0]);

    // Get the binding equation
    SX bexpr = readExpr(beq[1][0]);
    
    // Add binding equation
    addExplicitEquation(var.var(),bexpr);
  }
  
  // Mark all dependent variables
  for(vector<SX>::iterator it=explicit_var_.begin(); it!=explicit_var_.end(); ++it){
    it->setTemp(1);
  }
  
  // Add binding equations to constant variables lacking this
  vector<SX> beq_var;
  vector<SX> beq_exp;
  for(vector<Variable>::iterator it=vars_.begin(); it!=vars_.end(); ++it){
    if(it->getVariability()==CONSTANT && it->var().getTemp()!=1){
      // Save binding equation
      beq_var.push_back(it->var());
      beq_exp.push_back(it->getStart());
    }
  }
  
  // Unmark all dependent variables
  for(vector<SX>::iterator it=explicit_var_.begin(); it!=explicit_var_.end(); ++it){
    it->setTemp(0);
  }
  
  // Add to the beginning of the list of binding equations (in case some binding equations depend on them)
  explicit_var_.insert(explicit_var_.begin(),beq_var.begin(),beq_var.end());
  explicit_fcn_.insert(explicit_fcn_.begin(),beq_exp.begin(),beq_exp.end());
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
    implicit_fcn_.push_back(de_new);
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
      initial_eq_.push_back(readExpr(inode[i]));
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
    mterm.push_back(v);
  }
}

void FlatOCPInternal::addIntegrandObjectiveFunction(const XMLNode& onode){
  for(int i=0; i<onode.size(); ++i){
    const XMLNode& var = onode[i];
    SX v = readExpr(var);
    lterm.push_back(v);
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
      path_fcn_.push_back(ex-ub);
      path_min_.push_back(-numeric_limits<double>::infinity());
      path_max_.push_back(0.);
    } else if(constr_i.checkName("opt:ConstraintGeq")){
      SX ex = readExpr(constr_i[0]);
      SX lb = readExpr(constr_i[1]);
      path_fcn_.push_back(ex-lb);
      path_min_.push_back(0.);
      path_max_.push_back(numeric_limits<double>::infinity());
    } else if(constr_i.checkName("opt:ConstraintEq")){
      SX ex = readExpr(constr_i[0]);
      SX eq = readExpr(constr_i[1]);
      path_fcn_.push_back(ex-eq);
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
    return readVariable(node[0]).der(true);
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
  stream << "FMI parser (XML file: \"" << filename_ << "\")";
}

void FlatOCPInternal::print(ostream &stream) const{
  // Variables in the class hierarchy
  stream << "Variables" << endl;

  // Print the variables
  stream << "{" << endl;
  stream << "  t = " << t_ << endl;
  stream << "  x =  " << x_ << endl;
  stream << "  z =  " << z_ << endl;
  stream << "  p =  " << p_ << endl;
  stream << "  u =  " << u_ << endl;
  stream << "  y =  " << y_ << endl;
  stream << "}" << endl;
  stream << "Dimensions: "; 
  stream << "#x = " << x_.size() << ", ";
  stream << "#z = " << z_.size() << ", ";
  stream << "#p = " << p_.size() << ", ";
  stream << "#u = " << u_.size() << ", ";
  stream << endl << endl;
  
  // Print the differential-algebraic equation
  stream << "Dynamic equations" << endl;
  for(vector<SX>::const_iterator it=implicit_fcn_.begin(); it!=implicit_fcn_.end(); it++){
    stream << "0 == "<< *it << endl;
  }
  stream << endl;

  stream << "Initial equations" << endl;
  for(vector<SX>::const_iterator it=initial_eq_.begin(); it!=initial_eq_.end(); it++){
    stream << "0 == " << *it << endl;
  }
  stream << endl;

  // Print the explicit differential equations
/*  stream << "Differential equations (explicit)" << endl;
  for(vector<Variable>::const_iterator it=x_.begin(); it!=x_.end(); it++){
    SX de = it->rhs();
    if(!de->isNan())
      stream << "der(" << *it << ") == " << de << endl;
  }
  stream << endl;*/
  
  // Dependent equations
  stream << "Dependent equations" << endl;
  for(int i=0; i<explicit_var_.size(); ++i)
    stream << explicit_var_[i] << " == " << explicit_fcn_[i] << endl;
  stream << endl;

  // Mayer terms
  stream << "Mayer objective terms" << endl;
  for(int i=0; i<mterm.size(); ++i)
    stream << mterm[i] << endl;
  stream << endl;
  
  // Lagrange terms
  stream << "Lagrange objective terms" << endl;
  for(int i=0; i<lterm.size(); ++i)
    stream << lterm[i] << endl;
  stream << endl;
  
  // Constraint functions
  stream << "Constraint functions" << endl;
  for(int i=0; i<path_fcn_.size(); ++i)
    stream << path_min_[i] << " <= " << path_fcn_[i] << " <= " << path_max_[i] << endl;
  stream << endl;
  
  // Constraint functions
  stream << "Time horizon" << endl;
  stream << "t0 = " << t0_ << endl;
  stream << "tf = " << tf_ << endl;
  
}

void FlatOCPInternal::eliminateDependent(){
  cout << "eliminateDependent ..." << endl;
  double time1 = clock();

  Matrix<SX> v = explicit_var_;
  Matrix<SX> v_old = explicit_fcn_;
  
  implicit_fcn_= substitute(implicit_fcn_,v,v_old).data();
  initial_eq_= substitute(initial_eq_,v,v_old).data();
  path_fcn_    = substitute(path_fcn_,v,v_old).data();
  mterm   = substitute(mterm,v,v_old).data();
  lterm   = substitute(lterm,v,v_old).data();
  eliminated_dependents_ = true;

  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... eliminateDependent complete after " << dt << " seconds." << endl;
}

void FlatOCPInternal::addExplicitEquation(const Matrix<SX>& var, const Matrix<SX>& bind_eq, bool to_front){
  if(to_front){
    // Eliminate expression from the current binding equations
    Matrix<SX> explicit_fcn_eliminated = substitute(explicit_fcn_,var, bind_eq);

    // Save in the old location
    explicit_fcn_.swap(explicit_fcn_eliminated.data());
    
    // Add to list of equations
    explicit_var_.insert(explicit_var_.end(),var.begin(),var.end());
    explicit_fcn_.insert(explicit_fcn_.end(),bind_eq.begin(),bind_eq.end());
    
  } else {
    // Eliminate previous binding equations from the expression
    Matrix<SX> bind_eq_eliminated = substitute(bind_eq, explicit_var_, explicit_fcn_);
    
    // Add to list of equations
    explicit_var_.insert(explicit_var_.end(),var.begin(),var.data().end());
    explicit_fcn_.insert(explicit_fcn_.end(),bind_eq_eliminated.begin(),bind_eq_eliminated.data().end());
  }
}

void FlatOCPInternal::sortType(){
  // Clear variables
  x_.clear();
  z_.clear();
  u_.clear();
  p_.clear();
  y_.clear();
  
  // Mark all dependent variables
  for(vector<SX>::iterator it=explicit_var_.begin(); it!=explicit_var_.end(); ++it){
    it->setTemp(1);
  }
  
  // Get implicit variables
  bool find_implicit = implicit_var_.size() != implicit_fcn_.size();
  if(find_implicit){
    implicit_var_.clear();
    implicit_var_.reserve(implicit_fcn_.size());
  }
  
  // Loop over variables
  for(vector<Variable>::iterator it=vars_.begin(); it!=vars_.end(); ++it){
    // If not dependent
    if(it->highest().getTemp()!=1){
      // Try to determine the type
      if(it->getVariability() == PARAMETER){
        if(it->getFree()){
          p_.push_back(*it); 
        } else {
          casadi_assert(0);
        }
      } else if(it->getVariability() == CONTINUOUS) {
        if(it->getCausality() == INTERNAL){
          if(it->isDifferential()){
            x_.push_back(*it);
          } else {
            z_.push_back(*it);
          }
          
          // Add to list of implicit variables
          if(find_implicit) implicit_var_.push_back(it->highest());

        } else if(it->getCausality() == INPUT){
          u_.push_back(*it);
        }
      } else if(it->getVariability() == CONSTANT){
        casadi_assert(0);
      }
    } else {
      y_.push_back(*it);
    }
  }

  // Unmark all dependent variables
  for(vector<SX>::iterator it=explicit_var_.begin(); it!=explicit_var_.end(); ++it){
    it->setTemp(0);
  }
  
  // Assert consistent number of equations and variables
/*  cout << implicit_var_ << endl;
  cout << implicit_fcn_ << endl;
  cout << implicit_var_.size() << endl;
  cout << implicit_fcn_.size() << endl;*/
  casadi_assert(implicit_var_.size() == implicit_fcn_.size());
}


void FlatOCPInternal::scaleVariables(){
  cout << "Scaling variables ..." << endl;
  double time1 = clock();
  
  //Make sure that the variables has not already been scaled
  casadi_assert(!scaled_variables_);

  // Variables
  Matrix<SX> t = t_;
  Matrix<SX> x = var(x_);
  Matrix<SX> xdot = der(x_);
  Matrix<SX> z = var(z_);
  Matrix<SX> p = var(p_);
  Matrix<SX> u = var(u_);
  
  // Collect all the variables
  Matrix<SX> v;
  append(v,t);
  append(v,x);
  append(v,xdot);
  append(v,z);
  append(v,p);
  append(v,u);
  
  // Nominal values
  Matrix<SX> t_n = 1.;
  Matrix<SX> x_n = getNominal(x_);
  Matrix<SX> xdot_n = getNominal(x_);
  Matrix<SX> z_n = getNominal(z_);
  Matrix<SX> p_n = getNominal(p_);
  Matrix<SX> u_n = getNominal(u_);
  
  // Get all the old variables in expressed in the nominal ones
  Matrix<SX> v_old;
  append(v_old,t*t_n);
  append(v_old,x*x_n);
  append(v_old,xdot*xdot_n);
  append(v_old,z*z_n);
  append(v_old,p*p_n);
  append(v_old,u*u_n);
  
  // Temporary variable
  Matrix<SX> temp;

  // Substitute equations
  explicit_fcn_= substitute(explicit_fcn_,v,v_old).data();
  implicit_fcn_= substitute(implicit_fcn_,v,v_old).data();
  initial_eq_= substitute(initial_eq_,v,v_old).data();
  path_fcn_    = substitute(path_fcn_,v,v_old).data();
  mterm   = substitute(mterm,v,v_old).data();
  lterm   = substitute(lterm,v,v_old).data();
  
  scaled_variables_ = true;
  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... variable scaling complete after " << dt << " seconds." << endl;
}
    
void FlatOCPInternal::scaleEquations(){
  
  // Make sure that the equations has not already been scaled
  casadi_assert(!scaled_equations_);
  
  // Make sure that the dependents have been eliminated
  casadi_assert(eliminated_dependents_);
  
  // Make sure that the variables have been scaled
  casadi_assert(scaled_variables_);
  
  // Quick return if no implicit equations
  if(implicit_fcn_.empty())
    return;

  cout << "Scaling equations ..." << endl;
  double time1 = clock();

  // Variables
  enum Variables{T,X,XDOT,Z,P,U,NUM_VAR};
  vector<Matrix<SX> > v(NUM_VAR); // all variables
  v[T] = t_;
  v[X] = var(x_);
  v[XDOT] = der(x_);
  v[Z] = var(z_);
  v[P] = var(p_);
  v[U] = var(u_);

  // Create the jacobian of the implicit equations with respect to [x,z,p,u] 
  Matrix<SX> xz;
  append(xz,v[X]);
  append(xz,v[Z]);
  append(xz,v[P]);
  append(xz,v[U]);
  SXFunction fcn = SXFunction(xz,implicit_fcn_);
  SXFunction J(v,fcn.jac());

  // Evaluate the Jacobian in the starting point
  J.init();
  J.setInput(0.0,T);
  J.setInput(getStart(x_,true),X);
  J.input(XDOT).setAll(0.0);
  J.setInput(getStart(z_,true),Z);
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
      cout << "Warning: Could not generate a scaling factor for equation " << i << "(0 == " << implicit_fcn_[i] << "), selecting 1." << endl;
      scale[i]=1.;
    }
  }
  
  // Scale the equations
  for(int i=0; i<implicit_fcn_.size(); ++i){
    implicit_fcn_[i] /= scale[i];
  }
  
  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... equation scaling complete after " << dt << " seconds." << endl;
  scaled_equations_ = true;
}

void FlatOCPInternal::sortBLT(bool with_x){
  cout << "BLT sorting ..." << endl;
  double time1 = clock();
  
  // Sparsity pattern
  CRSSparsity sp;
  
  if(with_x){
    // inverse time constant
    SX invtau("invtau");

    // Replace x with invtau*xdot in order to get a Jacobian which also includes x
    SXMatrix implicit_fcn_with_x = substitute(implicit_var_,var(x_),invtau*SXMatrix(var(x_)));
    
    // Create Jacobian in order to find the sparsity
    SXFunction fcn(implicit_var_,implicit_fcn_with_x);
    Matrix<SX> J = fcn.jac();
    sp = J.sparsity();
  } else {
    // Create Jacobian in order to find the sparsity
    SXFunction fcn(implicit_var_,implicit_fcn_);
    Matrix<SX> J = fcn.jac();
    sp = J.sparsity();
  }
  
  // BLT transformation
  std::vector<int> blt_rowperm; // row permutations
  std::vector<int> blt_colperm; // column permutations
  std::vector<int> blt_rowblock;  // block k is rows r[k] to r[k+1]-1
  std::vector<int> blt_colblock;  // block k is cols s[k] to s[k+1]-1
  int blt_nb;
  std::vector<int> blt_coarse_rowblock;  // coarse row decomposition
  std::vector<int> blt_coarse_colblock;  //coarse column decomposition

  blt_nb = sp.dulmageMendelsohn(blt_rowperm,blt_colperm,blt_rowblock,blt_colblock,blt_coarse_rowblock,blt_coarse_colblock);

  // Permute equations
  vector<SX> implicit_fcn_new(implicit_fcn_.size());
  for(int i=0; i<implicit_fcn_.size(); ++i){
    implicit_fcn_new[i] = implicit_fcn_[blt_rowperm[i]];
  }
  implicit_fcn_new.swap(implicit_fcn_);
  
  // Permute variables
  vector<SX> implicit_var_new(implicit_var_.size());
  for(int i=0; i<implicit_var_.size(); ++i){
    implicit_var_new[i]= implicit_var_[blt_colperm[i]];
  }
  implicit_var_new.swap(implicit_var_);
  
  // Save blocks
  rowblock_ = blt_rowblock;
  colblock_ = blt_colblock;
  nb_ = blt_nb;
  
  blt_sorted_ = true;

  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... BLT sorting complete after " << dt << " seconds." << endl;
}

void FlatOCPInternal::makeExplicit(){
  casadi_assert_message(blt_sorted_,"OCP has not been BLT sorted, call sortBLT()");

  cout << "Making explicit..." << endl;
  double time1 = clock();
  
  // Create Jacobian
  SXFunction fcn(implicit_var_,implicit_fcn_);
  SXMatrix J = fcn.jac();

  // Mark the algebraic variables
  for(int i=0; i<z_.size(); ++i){
    z_[i].var().setTemp(i+1);
  }
  
  // Get initial values for all implicit variables
  vector<double> implicit_var_guess(implicit_var_.size(),0);
  for(int i=0; i<implicit_var_.size(); ++i){
    int ind = implicit_var_[i].getTemp()-1;
    if(ind>=0){
      implicit_var_guess[i] = z_[ind].getStart()/z_[ind].getNominal();
    }
  }
  
  // Unmark the algebraic variables
  for(int i=0; i<z_.size(); ++i){
    z_[i].var().setTemp(0);
  }
  
  // Block variables and equations
  vector<SX> vb, fb;

  // Save old number of explicit variables
  int old_nexp = explicit_var_.size();

  // New implicit equation and variables
  vector<SX> implicit_var_new, implicit_fcn_new;
    
  // Loop over blocks
  for(int b=0; b<nb_; ++b){
    
    // Block size
    int bs = rowblock_[b+1] - rowblock_[b];
    
/*    cout << "block " << b << ": size = " << bs << endl;*/
    
    // Get local variables
    vb.clear();
    for(int i=colblock_[b]; i<colblock_[b+1]; ++i)
      vb.push_back(implicit_var_[i]);

    // Get local equations
    fb.clear();
    for(int i=rowblock_[b]; i<rowblock_[b+1]; ++i)
      fb.push_back(implicit_fcn_[i]);

    // Get local Jacobian
    SXMatrix Jb = J(range(rowblock_[b],rowblock_[b+1]),range(colblock_[b],colblock_[b+1]));
    if(dependsOn(Jb,vb)){
      
      // Guess for the solution
      SXMatrix x_k(vb.size(),1,0);
      int offset = rowblock_[b];
      for(int i=0; i<x_k.size(); ++i){
        x_k.at(i) = implicit_var_guess[offset+i];
      }
      
      // Make Newton iterations
      bool exact_newton = true;
      SXFunction newtonIter; // Newton iteration function
      
      if(exact_newton){
        // Use exact Newton (i.e. update Jacobian in every iteration)
        newtonIter = SXFunction(vb,vb-solve(Jb,SXMatrix(fb)));
      } else {
        // Evaluate the Jacobian
        SXMatrix Jb0 = substitute(Jb,vb,x_k);
        
        // Use quasi-Newton
        newtonIter = SXFunction(vb,vb-solve(Jb0,SXMatrix(fb)));
      }
      newtonIter.init();
      
      // Make Newton iterations
      const int n_newton = 3; // NOTE: the number of newton iterations is fixed!
      for(int k=0; k<n_newton; ++k){
        x_k = newtonIter.eval(x_k);
      }
      
      cout << "Using " << (exact_newton ? "an exact " : "a quasi-") << "Newton iteration with "<< n_newton << " iterations to solve block " << b << " for the " << vb.size() << " variables " << vb << endl;
      cout << "the implicit equation has " << countNodes(fb) << " nodes" << endl;
      cout << "the Newton algorithm has " << newtonIter.algorithm().size() << " nodes" << endl;
      cout << "the explicit expression has " << countNodes(x_k) << " nodes" << endl;

      // Add binding equations
      addExplicitEquation(vb,x_k);
      
      // TODO: implement tearing
      casadi_assert_message(0,"not implemented");
      if(false){ // not ready (but make sure it compiles)
      
        // Find out which dependencies are nonlinear and which are linear
        Matrix<int> Jb_lin(Jb.sparsity());
        for(int i=0; i<Jb.size(); ++i){
          Jb_lin.at(i) = dependsOn(Jb.at(i),vb) ? 2 : 1;
        }
        
        // Loop over rows (equations)
        for(int i=0; i<Jb_lin.size1(); ++i){
          // Number of variables appearing linearly
          int n_lin = 0;
          int j_lin = -1; // index of a variable appearing linearily
          
          // Loop over non-zero elements
          for(int el=Jb_lin.rowind(i); el<Jb_lin.rowind(i+1); ++el){
            // Column (variable)
            int j = Jb_lin.col(el);
            
            // Count number of variables appearing linearly
            if(Jb_lin.at(el)==1){
              j_lin = j;
              n_lin ++;
            }
          }
          
          // Make causal if only one variable appears linearily
          if(n_lin==1){
            // 
            
          }
        }
      }
                  
      // Cannot solve for vb, add to list of implicit equations
//      implicit_var_new.insert(implicit_var_new.end(),vb.begin(),vb.end());
//      implicit_fcn_new.insert(implicit_fcn_new.end(),fb.begin(),fb.end());
      
//      cout << "added " << fb << " and " << vb << " to list of implicit equations (" << vb.size() << " equations)" << endl;
      
    } else {
      
      // Divide fb into a part which depends on vb and a part which doesn't according to "fb == prod(Jb,vb) + fb_res"
      SXMatrix fb_res = substitute(fb,vb,SXMatrix(vb.size(),1,0));
      SXMatrix fb_exp;
      
      // Solve for vb
      if (bs <= 3){
        // Calculate inverse and multiply for very small matrices
        fb_exp = prod(inv(Jb),-fb_res);
      } else {
        // QR factorization
        fb_exp = solve(Jb,-fb_res);
      }

      // Add explicit equation
      addExplicitEquation(vb,fb_exp);
    }
  }

  // Update implicit equations
  implicit_var_new.swap(implicit_var_);
  implicit_fcn_new.swap(implicit_fcn_);

  // Mark the variables made explicit
  for(vector<SX>::iterator it=explicit_var_.begin()+old_nexp; it!=explicit_var_.end(); ++it){
    it->setTemp(1);
  }
  
  // New algebraic variables
  vector<Variable> z_new;

  // Loop over algebraic variables
  for(vector<Variable>::iterator it=z_.begin(); it!=z_.end(); ++it){
    // Check if marked
    if(it->var().getTemp()){
      // Make dependent
      y_.push_back(*it);
      
      // If upper or lower bounds are finite, add path constraint
      if(!isinf(it->getMin()) || !isinf(it->getMax())){
        path_fcn_.push_back(it->var());
        path_min_.push_back(it->getMin()/it->getNominal());
        path_max_.push_back(it->getMax()/it->getNominal());
      }
    } else {
      z_new.push_back(*it);
    }
  }
  
  // Update new z_
  z_.swap(z_new);

  // Unark the variables made explicit
  for(vector<SX>::iterator it=explicit_var_.begin()+old_nexp; it!=explicit_var_.end(); ++it){
    it->setTemp(0);
  }
  
  // Eliminate the dependents
  cout << "eliminating dependents" << endl;
  eliminateDependent();

  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... makeExplicit complete after " << dt << " seconds." << endl;
}

void FlatOCPInternal::createFunctions(bool create_dae, bool create_ode, bool create_quad){
  cout << "createFunctions ..." << endl;
  double time1 = clock();

  // no quad if no quadrature states
  if(lterm.empty()) create_quad=false;

  if(create_ode){
    // Create an explicit ODE
    sortBLT(false);

    // Make the OCP explicit by eliminated all the algebraic states
    makeExplicit();
    
    // The equation could be made explicit symbolically
    if(implicit_fcn_.empty()){
      casadi_assert(z_.empty());
      
      // Mark the explicit variables
      for(int i=0; i<explicit_var_.size(); ++i){
        explicit_var_[i].setTemp(i+1);
      }

      // Find the ODE rhs
      vector<SX> ode(x_.size());
      for(int i=0; i<x_.size(); ++i){
        int ind = x_[i].der().getTemp()-1;
        casadi_assert(ind>=0);
        ode[i] = explicit_fcn_[ind];
      }
  
      // Find the dependency equations
      vector<SX> dep(y_.size());
      for(int i=0; i<y_.size(); ++i){
        int ind = y_[i].var().getTemp()-1;
        casadi_assert(ind>=0);
        dep[i] = explicit_fcn_[ind];
      }
      
      // Unmark the explicit variables
      for(int i=0; i<explicit_var_.size(); ++i){
        explicit_var_.at(i).setTemp(0);
      }

      // Evaluate constant expressions
/*      Matrix<SX> ode_elim = evaluateConstants(ode);*/
      Matrix<SX> ode_elim = ode;
      
      // ODE right hand side function
      vector<SXMatrix> ode_in(DAE_NUM_IN);
      ode_in[DAE_T] = t_;
      ode_in[DAE_Y] = var(x_);
      ode_in[DAE_P] = var(u_);
      oderhs_ = SXFunction(ode_in,ode_elim);
      
      // Dependency function
      output_fcn_ = SXFunction(ode_in,dep);
      
      // ODE quadrature function
      if(create_quad){
        quadrhs_ = SXFunction(ode_in,lterm[0]);
        
        // ODE right hand side with lterm
        append(ode_in[DAE_Y],symbolic("lterm"));
        append(ode_elim,SXMatrix(lterm[0]));
        oderhs_with_lterm_ = SXFunction(ode_in,ode_elim);
      }
      
    } else {
      cout << implicit_fcn_.size() << " implicit equations" << endl;
      
      
    // A Newton algorithm is necessary
    casadi_assert(0);
    
//     # Create an implicit function residual
//     impres_in = (DAE_NUM_IN+1) * [[]]
//     impres_in[0] = xdot+z
//     impres_in[1+DAE_T] = t
//     impres_in[1+DAE_Y] = x
//     impres_in[1+DAE_P] = u
//     impres = SXFunction(impres_in,[ocp.implicit_fcn_])
//     impres.setOption("number_of_fwd_dir",len(x)+1)
// 
//     impres.init()
//     impres.setInput(xdot0+z0, 0)
//     impres.setInput(0., 1+DAE_T)
//     impres.setInput(x0, 1+DAE_Y)
//     impres.setInput(u0, 1+DAE_P)
//     impres.evaluate()
// 
//     # Create an implicit function (KINSOL)
//     impsolver = KinsolSolver(impres)
//     linsol = CSparse(CRSSparsity())
//     impsolver.setLinearSolver(linsol)
// 
//     impsolver.setOption("linear_solver","user_defined")
//     impsolver.setOption("max_krylov",100)
//     impsolver.setOption("number_of_fwd_dir",len(x)+1)
//     impsolver.setOption("number_of_adj_dir",0)
//     impsolver.init()
// 
//     impsolver.setInput(0., DAE_T)
//     impsolver.setInput(x0, DAE_Y)
//     impsolver.setInput(u0, DAE_P)
//     impsolver.output().set(xdot0+z0)
//     impsolver.evaluate(1,0)
//     #print "implicit function residual", impsolver.output()
// 
//     # Create the ODE residual
//     ode_in_mx = DAE_NUM_IN * [[]]
//     ode_in_mx[DAE_T] = MX("T")
//     ode_in_mx[DAE_Y] = MX("X",len(x))
//     ode_in_mx[DAE_P] = MX("U",len(u))
//     [xdot_z] = impsolver.call(ode_in_mx)
//     ode = MXFunction(ode_in_mx,[xdot_z[0:len(x)]])
//     ode.setOption("number_of_fwd_dir",len(x)+1)
// 
//     # ODE quadrature function
//     ode_in = DAE_NUM_IN * [[]]
//     ode_in[DAE_T] = t
//     ode_in[DAE_Y] = x
//     ode_in[DAE_P] = u
//     ode_lterm = SXFunction(ode_in,[[ocp.lterm[0]/1e3]])
//     ode_lterm.setOption("number_of_fwd_dir",len(x)+1)
      
    }
  } else if(create_dae) {
    // BLT sorting
    sortBLT(true);

    // Scale the equations
    scaleEquations();
    
    // Time deriative of the algebraic variables
    SXMatrix zdot = symbolic("zdot",z_.size());
    
    // State
    SXMatrix xz = vertcat(SXMatrix(var(x_)),SXMatrix(var(z_)));
    
    // State derivative
    SXMatrix xzdot = vertcat(SXMatrix(der(x_)),zdot);
    
    // DAE residual arguments
    vector<SXMatrix> dae_in(DAE_NUM_IN);
    dae_in[DAE_T] = t_;
    dae_in[DAE_Y] = xz;
    dae_in[DAE_YDOT] = xzdot;
    dae_in[DAE_P] = var(u_);
    
    // DAE residual
    daeres_ = SXFunction(dae_in,implicit_fcn_);

    // DAE quadrature function
    if(create_quad){
      quadrhs_ = SXFunction(dae_in,lterm[0]/1e3);
    }
  }

  // Time deriative of the algebraic variables
  SXMatrix zdot = symbolic("zdot",z_.size());
    
  // State
  SXMatrix xz = vertcat(SXMatrix(var(x_)),SXMatrix(var(z_)));
    
  // State derivative
  SXMatrix xzdot = vertcat(SXMatrix(der(x_)),zdot);
  
  // Mayer objective function
  SXMatrix xf = symbolic("xf",x_.size()+lterm.size(),1);
  costfcn_ = SXFunction(xf, xf.data().back());

  // Path constraint function
  vector<SXMatrix> cfcn_in(DAE_NUM_IN);
  cfcn_in[DAE_T] = t_;
  cfcn_in[DAE_Y] = xz;
  cfcn_in[DAE_YDOT] = xzdot;
  append(cfcn_in[DAE_Y],symbolic("lterm")); // FIXME
  cfcn_in[DAE_P] = var(u_);
  pathfcn_ = SXFunction(cfcn_in,path_fcn_);
    
  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... createFunctions complete after " << dt << " seconds." << endl;
}


void FlatOCPInternal::makeSemiExplicit(){
  
  
//     fcn = SXFunction([v_new],[dae_new])
//   J = fcn.jac()
//   #J.printDense()
// 
//   # Cumulative variables and definitions
//   vb_cum = SXMatrix()
//   def_cum = SXMatrix()
// 
//   for b in range(blt_nb):
//     Jb = J[blt_rowblock[b]:blt_rowblock[b+1],blt_colblock[b]:blt_colblock[b+1]]
//     vb = v_new[blt_colblock[b]:blt_colblock[b+1]]
//     fb = dae_new[blt_rowblock[b]:blt_rowblock[b+1]]
//     
//     # Block size
//     bs = blt_rowblock[b+1] - blt_rowblock[b]
// 
//     #print "block ", b,
// 
//     if dependsOn(Jb,vb):
//       # Cannot solve for vb, add to list of implicit equations
//       raise Exception("Not implemented")
//       vb_cum.append(vb)
//       def_cum.append(vb)
//       
//     else:
//       # Divide fb into a part which depends on vb and a part which doesn't according to "fb == prod(Jb,vb) + fb_res"
//       fb_res = substitute(fb,vb,SXMatrix(len(vb),1,SX(0)))
//       
//       # Solve for vb
//       if bs <= 3:
//         # Calculate inverse and multiply for very small matrices
//         fb_exp = dot(inv(Jb),-fb_res)
//       else:
//         # QR factorization
//         fb_exp = solve(Jb,-fb_res)
//         
//       # Substitute variables that have already been defined
//       if not def_cum.empty():
//         fb_exp = substitute(fb_exp,vb_cum,def_cum)
//       
//       append(vb_cum,vb)
//       append(def_cum,fb_exp)
  
  
  throw CasadiException("FlatOCPInternal::makeSemiExplicit: Commented out");
#if 0  
  // Move the fully implicit dynamic equations to the list of algebraic equations
  algeq.insert(algeq.end(), dyneq.begin(), dyneq.end());
  dyneq.clear();
    
  // Introduce new explicit differential equations describing the relation between states and state derivatives
  xd.insert(xd.end(), x.begin(), x.end());
  diffeq.insert(diffeq.end(), xdot.begin(), xdot.end());
  
  // Put the state derivatives in the algebraic state category (inefficient!!!)
  xa.insert(xa.end(), xdot.begin(), xdot.end());

  // Remove from old location
  xdot.clear();
  x.clear();
#endif
}

void FlatOCPInternal::makeAlgebraic(const Variable& v){
  // Find variable among the variables
  for(vector<Variable>::iterator it=x_.begin(); it!=x_.end(); ++it){
    if(it->get()==v.get()){
      // Add dependent equation, xdot==0
      addExplicitEquation(v.der(),0.0,true);
            
      // Add to list of algebraic variables
      z_.push_back(v);
      
      // Remove from list of differential variables
      x_.erase(it);

      // Replace xdot with x (z) in the list of implicitly defined variables
      for(vector<SX>::iterator jt=implicit_var_.begin(); jt!=implicit_var_.end(); ++jt){
        if(jt->get()==v.der().get()){
          *jt = v.var();
          return;
        }
      }

      // Error if this point reached
      throw CasadiException("variable could not be found in list of implicitly defined variables");
    }
  }
  
  // Error if this point reached
  throw CasadiException("v not in list of differential states");
}

void FlatOCPInternal::findConsistentIC(){
  // Evaluate the ODE functions
  oderhs_.init();
  oderhs_.setInput(0.0,DAE_T);
  oderhs_.setInput(getStart(x_,true),DAE_Y);
  oderhs_.setInput(getStart(u_,true),DAE_P);
  oderhs_.evaluate();
  
  // Save to the variables
  for(int i=0; i<x_.size(); ++i){
    double xdot0 = oderhs_.output().at(i) * x_[i].getNominal();
    x_[i].setDerivativeStart(xdot0);
  }
  
  // Evaluate the output functions
  output_fcn_.init();
  output_fcn_.setInput(0.0,DAE_T);
  output_fcn_.setInput(getStart(x_,true),DAE_Y);
  output_fcn_.setInput(getStart(u_,true),DAE_P);
  output_fcn_.evaluate();
  
  // Save to the variables
  for(int i=0; i<y_.size(); ++i){
    double z0 = output_fcn_.output().at(i) * y_[i].getNominal();
    y_[i].setStart(z0);
  }
}

SX FlatOCPInternal::getExplicit(const SX& v) const{
  SXMatrix x = v;
  x = substitute(x,explicit_var_,explicit_fcn_);
  return x.toScalar();
}

Variable& FlatOCPInternal::variable(const std::string& name){
  // Find the variable
  map<string,int>::iterator it = varname_.find(name);
  if(it==varname_.end()){
    stringstream ss;
    ss << "No such variable: \"" << name << "\".";
    throw CasadiException(ss.str());
  }
  
  // Return the variable
  return vars_.at(it->second);
}

void FlatOCPInternal::addVariable(const std::string& name, const Variable& var){
  // Try to find the name
  map<string,int>::iterator it = varname_.find(name);
  if(it!=varname_.end()){
    stringstream ss;
    ss << "Variable \"" << name << "\" has already been added.";
    throw CasadiException(ss.str());
  }
  
  vars_.push_back(var);
  varname_[name] = vars_.size()-1;
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


} // namespace OptimalControl
} // namespace CasADi
