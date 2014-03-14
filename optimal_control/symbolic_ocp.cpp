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

#include "symbolic/std_vector_tools.hpp"
#include "external_packages/tinyxml/tinyxml.h"

#include "../symbolic/casadi_exception.hpp"
#include "../symbolic/std_vector_tools.hpp"
#include "variable_tools.hpp"
#include "../symbolic/matrix/matrix_tools.hpp"
#include "../symbolic/sx/sx_tools.hpp"
#include "../symbolic/fx/integrator.hpp"
#include "../symbolic/casadi_calculus.hpp"

using namespace std;
namespace CasADi{

SymbolicOCP::SymbolicOCP(){
  t = SX::sym("t");
  t0 = t0_guess = numeric_limits<double>::quiet_NaN();
  tf = tf_guess = numeric_limits<double>::quiet_NaN();
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
      string name        = vnode.getAttribute("name");
      int valueReference;
      vnode.readAttribute("valueReference",valueReference);
      string variability = vnode.getAttribute("variability");
      string causality   = vnode.getAttribute("causality");
      string alias       = vnode.getAttribute("alias");
      
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
          props.readAttribute("unit",var.unit(),false);
          props.readAttribute("displayUnit",var.displayUnit(),false);
          props.readAttribute("min",var.min(),false);
          props.readAttribute("max",var.max(),false);
          props.readAttribute("start",var.start(),false);
          props.readAttribute("nominal",var.nominal(),false);
          props.readAttribute("free",var.free(),false);
          props.readAttribute("initialGuess",var.initialGuess(),false);
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
      bool has_der = false;
      SXElement bexpr = readExpr(beq[1][0],has_der,eliminate_dependent);
      casadi_assert(!has_der);
      
      // Add binding equation
      var.setBinding(bexpr);
      y.push_back(var); // legacy
      dep.append(bexpr); // legacy
    }
    
    // Resort the dependant parameters
    sortDependentParameters();
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
      bool has_der = false;
      SXElement de_new = readExpr(dnode[0],has_der,eliminate_dependent);
      if(has_der){
        ode.append(de_new);
      } else {
        alg.append(de_new);
      }
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
        bool has_der = false;
        initial.append(readExpr(inode[i],has_der,eliminate_dependent));
      }
    }
  }
  
  // **** Add optimization ****
  if(document[0].hasChild("opt:Optimization")){
    
    // Get a reference to the DynamicEquations node
    const XMLNode& opts = document[0]["opt:Optimization"];
    
    // Start time
    const XMLNode& intervalStartTime = opts["opt:IntervalStartTime"];
    if(intervalStartTime.hasChild("opt:Value"))
      intervalStartTime["opt:Value"].getText(t0);
    if(intervalStartTime.hasChild("opt:Free"))
      intervalStartTime["opt:Free"].getText(t0_free);
    if(intervalStartTime.hasChild("opt:InitialGuess"))
      intervalStartTime["opt:InitialGuess"].getText(t0_guess);

    // Terminal time
    const XMLNode& IntervalFinalTime = opts["opt:IntervalFinalTime"];
    if(IntervalFinalTime.hasChild("opt:Value"))
      IntervalFinalTime["opt:Value"].getText(tf);
    if(IntervalFinalTime.hasChild("opt:Free"))
      IntervalFinalTime["opt:Free"].getText(tf_free);
    if(IntervalFinalTime.hasChild("opt:InitialGuess"))
      IntervalFinalTime["opt:InitialGuess"].getText(tf_guess);

    // Time points
    const XMLNode& tpnode = opts["opt:TimePoints"];
    tp.resize(tpnode.size());
    for(int i=0; i<tp.size(); ++i){
      // Get index
      int index;
      tpnode[i].readAttribute("index",index);
      
      // Get value
      double value;
      tpnode[i].readAttribute("value",value);
      tp[i] = value;
      
      // Allocate all the timed variables
      for(int k=0; k<tpnode[i].size(); ++k){
        string qn = qualifiedName(tpnode[i][k]);
        variable(qn).atTime(value,true);
      }
    }
    
    for(int i=0; i<opts.size(); ++i){
      
      // Get a reference to the node
      const XMLNode& onode = opts[i];

      // Get the type
      if(onode.checkName("opt:ObjectiveFunction")){ // mayer term
        try{
          // Add components
          for(int i=0; i<onode.size(); ++i){
            const XMLNode& var = onode[i];
            
            // If string literal, ignore
            if(var.checkName("exp:StringLiteral"))
              continue;
            
            // Read expression
            bool has_der = false;
            SXElement v = readExpr(var,has_der,eliminate_dependent);
            casadi_assert(!has_der);
            mterm.append(v);
          }
        } catch(exception& ex){
          if(verbose){
            throw CasadiException(std::string("addObjectiveFunction failed: ") + ex.what());
          }
        }
      } else if(onode.checkName("opt:IntegrandObjectiveFunction")){
        try{
          for(int i=0; i<onode.size(); ++i){
            const XMLNode& var = onode[i];

            // If string literal, ignore
            if(var.checkName("exp:StringLiteral"))
              continue;
            
            // Read expression
            bool has_der = false;
            SXElement v = readExpr(var,has_der,eliminate_dependent);
            lterm.append(v);
          }
        } catch(exception& ex){
          throw CasadiException(std::string("addIntegrandObjectiveFunction failed: ") + ex.what());
        }
      } else if(onode.checkName("opt:IntervalStartTime")) {
        // Ignore, treated above
      } else if(onode.checkName("opt:IntervalFinalTime")) {
        // Ignore, treated above
      } else if(onode.checkName("opt:TimePoints")) {
        // Ignore, treated above
      } else if(onode.checkName("opt:PointConstraints")) {
        bool has_der = false; // Should we check that this remains false, i.e. should state derivatives be allowed in constraints?

        for(int i=0; i<onode.size(); ++i){
          const XMLNode& constr_i = onode[i];
          if(constr_i.checkName("opt:ConstraintLeq")){
            SXElement ex = readExpr(constr_i[0],has_der,eliminate_dependent);
            SXElement ub = readExpr(constr_i[1],has_der,eliminate_dependent);
            point.append(ex-ub);
            point_min.append(-numeric_limits<double>::infinity());
            point_max.append(0.);
          } else if(constr_i.checkName("opt:ConstraintGeq")){
            SXElement ex = readExpr(constr_i[0],has_der,eliminate_dependent);
            SXElement lb = readExpr(constr_i[1],has_der,eliminate_dependent);
            point.append(ex-lb);
            point_min.append(0.);
            point_max.append(numeric_limits<double>::infinity());
          } else if(constr_i.checkName("opt:ConstraintEq")){
            SXElement ex = readExpr(constr_i[0],has_der,eliminate_dependent);
            SXElement eq = readExpr(constr_i[1],has_der,eliminate_dependent);
            point.append(ex-eq);
            point_min.append(0.);
            point_max.append(0.);
          } else {
            cerr << "unknown constraint type" << constr_i.getName() << endl;
            throw CasadiException("SymbolicOCP::addConstraints");
          }
        }
        
      } else if(onode.checkName("opt:Constraints") || onode.checkName("opt:PathConstraints")) {
        
        bool has_der = false; // Should we check that this remains false, i.e. should state derivatives be allowed in constraints?
        for(int i=0; i<onode.size(); ++i){
          const XMLNode& constr_i = onode[i];
          if(constr_i.checkName("opt:ConstraintLeq")){
            SXElement ex = readExpr(constr_i[0],has_der,eliminate_dependent);
            SXElement ub = readExpr(constr_i[1],has_der,eliminate_dependent);
            path.append(ex-ub);
            path_min.append(-numeric_limits<double>::infinity());
            path_max.append(0.);
          } else if(constr_i.checkName("opt:ConstraintGeq")){
            SXElement ex = readExpr(constr_i[0],has_der,eliminate_dependent);
            SXElement lb = readExpr(constr_i[1],has_der,eliminate_dependent);
            path.append(ex-lb);
            path_min.append(0.);
            path_max.append(numeric_limits<double>::infinity());
          } else if(constr_i.checkName("opt:ConstraintEq")){
            SXElement ex = readExpr(constr_i[0],has_der,eliminate_dependent);
            SXElement eq = readExpr(constr_i[1],has_der,eliminate_dependent);
            path.append(ex-eq);
            path_min.append(0.);
            path_max.append(0.);
          } else {
            cerr << "unknown constraint type" << constr_i.getName() << endl;
            throw CasadiException("SymbolicOCP::addConstraints");
          }
        }
        
      } else throw CasadiException(string("SymbolicOCP::addOptimization: Unknown node ")+onode.getName());
    }
  }
  
  // Make sure that the dimensions are consistent at this point
  casadi_assert_warning(x.size()==ode.size(),"The number of differential equations (equations involving differentiated variables) does not match the number of differential states.");
  casadi_assert_warning(z.size()==alg.size(),"The number of algebraic equations (equations not involving differentiated variables) does not match the number of algebraic variables.");
  casadi_assert(q.size()==quad.size());
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
    sortODE();
    sortALG();
  }

  // Scale the equations
  if(scale_equations)
    scaleEquations();
    
  // Transform to semi-explicit form
  if(make_explicit){
    casadi_assert(eliminate_dependent);
    
    // Save the old number of dependent states
    int ny = y.size();
    
    // Make the ODE explicit
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

SXElement SymbolicOCP::readExpr(const XMLNode& node, bool& has_der, bool elim_binding){
  const string& fullname = node.getName();
  if (fullname.find("exp:")== string::npos) {
    casadi_error("SymbolicOCP::readExpr: unknown - expression is supposed to start with 'exp:' , got " << fullname);
  }
  
  // Chop the 'exp:'
  string name = fullname.substr(4);

  // The switch below is alphabetical, and can be thus made more efficient, for example by using a switch statement of the first three letters, if it would ever become a bottleneck
  if(name.compare("Add")==0){
    return readExpr(node[0],has_der,elim_binding) + readExpr(node[1],has_der,elim_binding);
  } else if(name.compare("Acos")==0){
    return acos(readExpr(node[0],has_der,elim_binding));
  } else if(name.compare("Asin")==0){
    return asin(readExpr(node[0],has_der,elim_binding));
  } else if(name.compare("Atan")==0){
    return atan(readExpr(node[0],has_der,elim_binding));
  } else if(name.compare("Cos")==0){
    return cos(readExpr(node[0],has_der,elim_binding));
  } else if(name.compare("Der")==0){
    Variable v = readVariable(node[0]);
    v.setDifferential(true);
    has_der = true;
    return v.der();
  } else if(name.compare("Div")==0){
    return readExpr(node[0],has_der,elim_binding) / readExpr(node[1],has_der,elim_binding);
  } else if(name.compare("Exp")==0){
    return exp(readExpr(node[0],has_der,elim_binding));
  } else if(name.compare("Identifier")==0){
    return readVariable(node).var();
  } else if(name.compare("IntegerLiteral")==0){
    int val;
    node.getText(val);
    return val;
  } else if(name.compare("Instant")==0){
    double val;
    node.getText(val);
    return val;
  } else if(name.compare("Log")==0){
    return log(readExpr(node[0],has_der,elim_binding));
  } else if(name.compare("LogLt")==0){ // Logical less than
    return readExpr(node[0],has_der,elim_binding) < readExpr(node[1],has_der,elim_binding);
  } else if(name.compare("LogGt")==0){ // Logical less than
    return readExpr(node[0],has_der,elim_binding) > readExpr(node[1],has_der,elim_binding);
  } else if(name.compare("Mul")==0){ // Multiplication
    return readExpr(node[0],has_der,elim_binding) * readExpr(node[1],has_der,elim_binding);
  } else if(name.compare("Neg")==0){
    return -readExpr(node[0],has_der,elim_binding);
  } else if(name.compare("NoEvent")==0) {
    // NOTE: This is a workaround, we assume that whenever NoEvent occurs, what is meant is a switch
    int n = node.size();
    
    // Default-expression
    SXElement ex = readExpr(node[n-1],has_der,elim_binding);
    
    // Evaluate ifs
    for(int i=n-3; i>=0; i -= 2) ex = if_else(readExpr(node[i],has_der,elim_binding),readExpr(node[i+1],has_der,elim_binding),ex);
    
    return ex;
  } else if(name.compare("Pow")==0){
    return pow(readExpr(node[0],has_der,elim_binding),readExpr(node[1],has_der,elim_binding));
  } else if(name.compare("RealLiteral")==0){
    double val;
    node.getText(val);
    return val;
  } else if(name.compare("Sin")==0){
    return sin(readExpr(node[0],has_der,elim_binding));
  } else if(name.compare("Sqrt")==0){
    return sqrt(readExpr(node[0],has_der,elim_binding));
  } else if(name.compare("StringLiteral")==0){
    throw CasadiException(node.getText());
  } else if(name.compare("Sub")==0){
    return readExpr(node[0],has_der,elim_binding) - readExpr(node[1],has_der,elim_binding);
  } else if(name.compare("Tan")==0){
    return tan(readExpr(node[0],has_der,elim_binding));
  } else if(name.compare("Time")==0){
    return t.toScalar();
  } else if(name.compare("TimedVariable")==0){
    // Get the index of the time point
    int index;
    node.readAttribute("timePointIndex",index);
    return readVariable(node[0]).atTime(tp[index]);
  }

  // throw error if reached this point
  throw CasadiException(string("SymbolicOCP::readExpr: Unknown node: ") + name);
  
}

void SymbolicOCP::repr(std::ostream &stream) const{
  stream << "Flat OCP";
}

void SymbolicOCP::print(ostream &stream) const{
  stream << "Dimensions: "; 
  stream << "#x = " << x.size() << ", ";
  stream << "#z = " << z.size() << ", ";
  stream << "#q = " << q.size() << ", ";
  stream << "#y = " << y.size() << ", ";
  stream << "#pi = " << pi.size() << ", ";
  stream << "#pd = " << pd.size() << ", ";
  stream << "#pf = " << pf.size() << ", ";
  stream << "#ci =  " << ci.size() << ", ";
  stream << "#cd =  " << cd.size() << ", ";
  stream << "#u = " << u.size() << ", ";
  stream << endl << endl;

  // Variables in the class hierarchy
  stream << "Variables" << endl;

  // Print the variables
  stream << "{" << endl;
  stream << "  t = " << t << endl;
  stream << "  x = " << x << endl;
  stream << "  z =  " << z << endl;
  stream << "  q =  " << q << endl;
  stream << "  y =  " << y << endl;
  stream << "  pi =  " << pi << endl;
  stream << "  pd =  " << pd << endl;
  stream << "  pf =  " << pf << endl;
  stream << "  ci =  " << ci << endl;
  stream << "  cd =  " << cd << endl;
  stream << "  u =  " << u << endl;
  stream << "}" << endl;
  
  stream << "Differential equations" << endl;
  for(int k=0; k<x.size(); ++k){
    stream << "0 == " << ode.at(k) << endl;
  }
  stream << endl;

  stream << "Algebraic equations" << endl;
  for(int k=0; k<z.size(); ++k){
    stream << "0 == " << alg.at(k) << endl;
  }
  stream << endl;
  
  stream << "Quadrature equations" << endl;
  for(int k=0; k<q.size(); ++k){
    stream << q.at(k).der() << " == " << quad.at(k) << endl;
  }
  stream << endl;

  stream << "Initial equations" << endl;
  for(vector<SXElement>::const_iterator it=initial.begin(); it!=initial.end(); it++){
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
  
  // Path constraint functions
  stream << "Path constraint functions" << endl;
  for(int i=0; i<path.size(); ++i)
    stream << path_min.at(i) << " <= " << path.at(i) << " <= " << path_max.at(i) << endl;
  stream << endl;
  
  // Point constraint functions
  stream << "Point constraint functions" << endl;
  for(int i=0; i<point.size(); ++i)
    stream << point_min.at(i) << " <= " << point.at(i) << " <= " << point_max.at(i) << endl;
  stream << endl;
  
  // Constraint functions
  stream << "Time horizon" << endl;
  stream << "t0 = " << t0 << endl;
  stream << "tf = " << tf << endl;
  stream << "tp = " << tp << endl;
}

void SymbolicOCP::eliminateInterdependencies(){
  substituteInPlace(var(y),dep,false);
  
  // Make sure that the dependent variables have been properly eliminated from the dependent expressions
  casadi_assert(!dependsOn(dep,var(y)));
}

vector<SX> SymbolicOCP::substituteDependents(const vector<SX>& x) const{
  return substitute(x,vector<SX>(1,var(y)),vector<SX>(1,dep));
}

void SymbolicOCP::eliminateDependent(bool eliminate_dependents_with_bounds){
  // All the functions to be replaced
  vector<SX> fcn(7);
  fcn[0] = ode;
  fcn[1] = alg;
  fcn[2] = quad;
  fcn[3] = initial;
  fcn[4] = path;
  fcn[5] = mterm;
  fcn[6] = lterm;
  
  // Replace all at once
  vector<SX> fcn_new = substituteDependents(fcn);
  
  // Save the new expressions
  ode = fcn_new[0];
  alg = fcn_new[1];
  quad = fcn_new[2];
  initial = fcn_new[3];
  path    = fcn_new[4];
  mterm   = fcn_new[5];
  lterm   = fcn_new[6];
}

void SymbolicOCP::eliminateLagrangeTerms(){
  // Index for the names
  int ind = 0;
  // For every integral term in the objective function
  for(vector<SXElement>::iterator it=lterm.begin(); it!=lterm.end(); ++it){
    
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
    q.push_back(qv);

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
  x.insert(x.end(),q.begin(),q.end());
  q.clear();
  
  // Move the equations to the list of ODEs
  ode.append(quad);
  quad.clear();
}

void SymbolicOCP::scaleVariables(){
  cout << "Scaling variables ..." << endl;
  double time1 = clock();
  
  // Variables
  Matrix<SXElement> _x = var(x);
  Matrix<SXElement> _xdot = der(x);
  Matrix<SXElement> _z = var(z);
  Matrix<SXElement> _pi = var(pi);
  Matrix<SXElement> _pf = var(pf);
  Matrix<SXElement> _u = var(u);
  
  // Collect all the variables
  Matrix<SXElement> v;
  v.append(t);
  v.append(_x);
  v.append(_xdot);
  v.append(_z);
  v.append(_pi);
  v.append(_pf);
  v.append(_u);
  
  // Nominal values
  Matrix<SXElement> t_n = 1.;
  Matrix<SXElement> x_n = getNominal(x);
  Matrix<SXElement> z_n = getNominal(z);
  Matrix<SXElement> pi_n = getNominal(pi);
  Matrix<SXElement> pf_n = getNominal(pf);
  Matrix<SXElement> u_n = getNominal(u);
  
  // Get all the old variables in expressed in the nominal ones
  Matrix<SXElement> v_old;
  v_old.append(t*t_n);
  v_old.append(_x*x_n);
  v_old.append(_xdot*x_n);
  v_old.append(_z*z_n);
  v_old.append(_pi*pi_n);
  v_old.append(_pf*pf_n);
  v_old.append(_u*u_n);
  
  // Temporary variable
  Matrix<SXElement> temp;

  // Substitute equations
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
  
  cout << "Scaling equations ..." << endl;
  double time1 = clock();

  // Variables
  enum Variables{T,X,XDOT,Z,PI,PF,U,NUM_VAR};
  vector<Matrix<SXElement> > v(NUM_VAR); // all variables
  v[T] = t;
  v[X] = var(x);
  v[XDOT] = der(x);
  v[Z] = var(z);
  v[PI] = var(pi);
  v[PF] = var(pf);
  v[U] = var(u);

  // Create the jacobian of the implicit equations with respect to [x,z,p,u] 
  Matrix<SXElement> xz;
  xz.append(v[X]);
  xz.append(v[Z]);
  xz.append(v[PI]);
  xz.append(v[PF]);
  xz.append(v[U]);
  SXFunction fcn = SXFunction(xz,ode);
  SXFunction J(v,fcn.jac());

  // Evaluate the Jacobian in the starting point
  J.init();
  J.setInput(0.0,T);
  J.setInput(getStart(x,true),X);
  J.input(XDOT).setAll(0.0);
  J.setInput(getStart(z,true),Z);
  J.setInput(getStart(pi,true),PI);
  J.setInput(getStart(pf,true),PF);
  J.setInput(getStart(u,true),U);
  J.evaluate();
  
  // Get the maximum of every row
  Matrix<double> &J0 = J.output();
  vector<double> scale(J0.size1(),0.0); // scaling factors
  for(int cc=0; cc<J0.size2(); ++cc){
    // Loop over non-zero entries of the column
    for(int el=J0.colind(cc); el<J0.colind(cc+1); ++el){
      // Row
      int rr=J0.row(el);
      
      // The scaling factor is the maximum norm, ignoring not-a-number entries
      if(!isnan(J0.at(el))){
        scale[rr] = max(scale[rr],fabs(J0.at(el)));
      }
    }
  }
  
  // Make sure nonzero factor found
  for(int rr=0; rr<J0.size1(); ++rr){
    if(scale[rr]==0){
      cout << "Warning: Could not generate a scaling factor for equation " << rr << "(0 == " << ode.at(rr) << "), selecting 1." << endl;
      scale[rr]=1.;
    }
  }
  
  // Scale the equations
  for(int i=0; i<ode.size(); ++i){
    ode[i] /= scale[i];
  }
  
  double time2 = clock();
  double dt = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "... equation scaling complete after " << dt << " seconds." << endl;
}

void SymbolicOCP::sortODE(){
  // Quick return if no differential states
  if(x.empty()) return;

  // Find out which differential equation depends on which differential state
  SXFunction f(der(x),ode);
  f.init();
  Sparsity sp = f.jacSparsity();
  
  // BLT transformation
  vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
  sp.dulmageMendelsohn(rowperm,colperm,rowblock,colblock,coarse_rowblock,coarse_colblock);

  // Permute equations
  vector<SXElement> ode_new(ode.size());
  for(int i=0; i<ode.size(); ++i){
    ode_new[i] = ode.at(rowperm[i]);
  }
  ode_new.swap(ode.data());
  
  // Permute variables
  vector<Variable> x_new(x.size());
  for(int i=0; i<x.size(); ++i){
    x_new[i]= x[colperm[i]];
  }
  x_new.swap(x);
}

void SymbolicOCP::sortALG(){
  // Quick return if no algebraic states
  if(z.empty()) return;
  
  // Find out which algebraic equation depends on which algebraic state
  SXFunction f(var(z),alg);
  f.init();
  Sparsity sp = f.jacSparsity();
  
  // BLT transformation
  vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
  sp.dulmageMendelsohn(rowperm,colperm,rowblock,colblock,coarse_rowblock,coarse_colblock);

  // Permute equations
  vector<SXElement> alg_new(alg.size());
  for(int i=0; i<alg.size(); ++i){
    alg_new[i] = alg.at(rowperm[i]);
  }
  alg_new.swap(alg.data());
  
  // Permute variables
  vector<Variable> z_new(z.size());
  for(int i=0; i<z.size(); ++i){
    z_new[i]= z[colperm[i]];
  }
  z_new.swap(z);
}

void SymbolicOCP::sortDependentParameters(){
  // Quick return if no algebraic states
  if(pd.empty()) return;
  
  // Find out which dependent parameter depends on which binding equation
  SX v = var(pd);
  SXFunction f(v,v-binding(pd));
  f.init();
  Sparsity sp = f.jacSparsity();
  
  // BLT transformation
  vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
  sp.dulmageMendelsohn(rowperm,colperm,rowblock,colblock,coarse_rowblock,coarse_colblock);

  // Permute variables
  vector<Variable> pd_new(pd.size());
  for(int i=0; i<pd.size(); ++i){
    pd_new[i]= pd[colperm[i]];
  }
  pd_new.swap(pd);  
}

void SymbolicOCP::makeExplicit(){
  // Quick return if there are no differential states
  if(x.empty()) return;
  
  // Make sure that the ODE is not already explicit
  if(!dependsOn(ode,der(x))){
    casadi_warning("The ODE is already explicit");
    return;
  }
  
  // Write the ODE as a function of the state derivatives
  SXFunction f(der(x),ode);
  f.init();

  // Get the sparsity of the Jacobian which can be used to determine which variable can be calculated from which other
  Sparsity sp = f.jacSparsity();

  // BLT transformation
  vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
  int nb = sp.dulmageMendelsohn(rowperm,colperm,rowblock,colblock,coarse_rowblock,coarse_colblock);

  // Permute equations
  vector<SXElement> ode_new(ode.size());
  for(int i=0; i<ode.size(); ++i){
    ode_new[i] = ode.at(rowperm[i]);
  }
  ode_new.swap(ode.data());
  ode_new.clear();
  
  // Permute variables
  vector<Variable> x_new(x.size());
  for(int i=0; i<x.size(); ++i){
    x_new[i]= x[colperm[i]];
  }
  x_new.swap(x);
  x_new.clear();

  // Now write the sorted ODE as a function of the state derivatives
  f = SXFunction(der(x),ode);
  f.init();

  // Get the Jacobian
  SX J = f.jac();
  
  // Block variables and equations
  vector<Variable> xb, xdb, xab;
  vector<SXElement> fb;

  // Explicit ODE
  vector<SXElement> ode_exp;
  
  // Loop over blocks
  for(int b=0; b<nb; ++b){
    
    // Block size
    int bs = rowblock[b+1] - rowblock[b];
    
    // Get variables in the block
    xb.clear();
    for(int i=colblock[b]; i<colblock[b+1]; ++i){
      xb.push_back(x[i]);
    }

    // Get equations in the block
    fb.clear();
    for(int i=rowblock[b]; i<rowblock[b+1]; ++i)
      fb.push_back(ode.at(i));

    // Get local Jacobian
    SX Jb = J(range(rowblock[b],rowblock[b+1]),range(colblock[b],colblock[b+1]));

    // If Jb depends on xb, then the state derivative does not enter linearly in the ODE and we cannot solve for the state derivative
    casadi_assert_message(!dependsOn(Jb,der(xb)),"Cannot find an explicit expression for variable(s) " << xb);
      
    // Divide fb into a part which depends on vb and a part which doesn't according to "fb == mul(Jb,vb) + fb_res"
    SX fb_res = substitute(fb,der(xb),SX::zeros(xb.size())).data();
    SX fb_exp;
      
    // Solve for vb
    if (bs <= 3){
      // Calculate inverse and multiply for very small matrices
      fb_exp = mul(inv(Jb),-fb_res);
    } else {
      // QR factorization
      fb_exp = solve(Jb,-fb_res);
    }

    // Add to explicitly determined equations and variables
    ode_exp.insert(ode_exp.end(),fb_exp.data().begin(),fb_exp.data().end());
  }
  
  // Eliminate inter-dependencies
  SX ode_expmat = ode_exp;
  substituteInPlace(der(x),ode_expmat,false);
  ode_exp = ode_expmat.data();
  ode_exp.swap(ode.data());
}

void SymbolicOCP::eliminateAlgebraic(){
  // Quick return if there are no algebraic states
  if(z.empty()) return;
  
  // Write the algebraic equations as a function of the algebraic states
  SXFunction f(var(z),alg);
  f.init();

  // Get the sparsity of the Jacobian which can be used to determine which variable can be calculated from which other
  Sparsity sp = f.jacSparsity();

  // BLT transformation
  vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
  int nb = sp.dulmageMendelsohn(rowperm,colperm,rowblock,colblock,coarse_rowblock,coarse_colblock);

  // Permute equations
  vector<SXElement> alg_new(alg.size());
  for(int i=0; i<alg.size(); ++i){
    alg_new[i] = alg.at(rowperm[i]);
  }
  alg_new.swap(alg.data());
  alg_new.clear();
  
  // Permute variables
  vector<Variable> z_new(z.size());
  for(int i=0; i<z.size(); ++i){
    z_new[i]= z[colperm[i]];
  }
  z_new.swap(z);
  z_new.clear();

  // Rewrite the sorted algebraic equations as a function of the algebraic states
  f = SXFunction(var(z),alg);
  f.init();

  // Get the Jacobian
  SX J = f.jac();
  
  // Block variables and equations
  vector<Variable> zb;
  vector<SXElement> fb;

  // Variables where we have found an explicit expression and where we haven't
  vector<Variable> z_exp, z_imp;
  
  // Explicit and implicit equations
  vector<SXElement> f_exp, f_imp;
  
  // Loop over blocks
  for(int b=0; b<nb; ++b){
    
    // Block size
    int bs = rowblock[b+1] - rowblock[b];
    
    // Get local variables
    zb.clear();
    for(int i=colblock[b]; i<colblock[b+1]; ++i){
      zb.push_back(z[i]);
    }

    // Get local equations
    fb.clear();
    for(int i=rowblock[b]; i<rowblock[b+1]; ++i)
      fb.push_back(alg.at(i));

    // Get local Jacobian
    SX Jb = J(range(rowblock[b],rowblock[b+1]),range(colblock[b],colblock[b+1]));

    // If Jb depends on zb, then we cannot (currently) solve for it explicitly
    if(dependsOn(Jb,var(zb))){
      
      // Add the equations to the new list of algebraic equations
      f_imp.insert(f_imp.end(),fb.begin(),fb.end());
        
      // ... and the variables accordingly
      z_imp.insert(z_imp.end(),zb.begin(),zb.end());
      
    } else { // The variables that we wish to determine enter linearly
      
      // Divide fb into a part which depends on vb and a part which doesn't according to "fb == mul(Jb,vb) + fb_res"
      SX fb_res = substitute(fb,var(zb),SX::zeros(zb.size())).data();
      SX fb_exp;
      
      // Solve for vb
      if (bs <= 3){
        // Calculate inverse and multiply for very small matrices
        fb_exp = mul(inv(Jb),-fb_res);
      } else {
        // QR factorization
        fb_exp = solve(Jb,-fb_res);
      }

      // Add to explicitly determined equations and variables
      z_exp.insert(z_exp.end(),zb.begin(),zb.end());
      f_exp.insert(f_exp.end(),fb_exp.data().begin(),fb_exp.data().end());
    }
  }
  
  // Eliminate inter-dependencies in fb_exp
  SX f_expmat = f_exp;
  substituteInPlace(var(z_exp),f_expmat,false);
  f_exp = f_expmat.data();

  // Add to the beginning of the dependent variables (since the other dependent variable might depend on them)
  y.insert(y.begin(),z_exp.begin(),z_exp.end());
  dep = vertcat(SX(f_exp),dep);
  
  // Save new algebraic equations
  z = z_imp;
  alg = f_imp;
  
  // Eliminate new dependent variables from the other equations
  eliminateDependent();
}

void SymbolicOCP::makeAlgebraic(const std::string& name){
  makeAlgebraic(variable(name));
}

void SymbolicOCP::makeAlgebraic(const Variable& v){
  casadi_assert(0);
#if 0
  // Find variable among the explicit variables
  for(int k=0; k<x.size(); ++k){
    if(x[k].get()==v.get()){
      
      // Add to list of algebraic variables and to the list of algebraic equations
      z.push_back(v);
      alg.append(ode.at(k));
      
      // Remove from list of differential variables and the list of differential equations
      x.erase(x.begin()+k);
      vector<SXElement> ode_ = ode.data();
      ode_.erase(ode_.begin()+k);
      ode = ode_;

      // Successfull return
      return;
    }
  }
  
  // Find the variable among the implicit variables
  for(int k=0; k<xz.size(); ++k){
    if(xz[k].get()==v.get()){
      
      // Substitute the state derivative with zero
      dae = substitute(dae,xz[k].der(),0.0);

      // Remove the highest state derivative expression from the variable
      xz[k].setDifferential(false);

      // Successfull return
      return;
    }
  }
  
  // Error if this point reached
  throw CasadiException("v not a differential state");
#endif
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
  
  // Add to the map of all variables
  varmap_[name] = var;
  
  // Sort by category
  switch(var.getCategory()){
    case CAT_DERIVATIVE:
      // Skip derivatives
      break;
    case CAT_STATE:
      x.push_back(var);
      break;
    case CAT_DEPENDENT_CONSTANT:
      cd.push_back(var);
      break;
    case CAT_INDEPENDENT_CONSTANT:
      ci.push_back(var);
      break;
    case CAT_DEPENDENT_PARAMETER:
      pd.push_back(var);
      break;
    case CAT_INDEPENDENT_PARAMETER:
      if(var.getFree()){
        pf.push_back(var);
      } else {
        pi.push_back(var);
      }
      break;
    case CAT_ALGEBRAIC:
      if(var.getCausality() == INTERNAL){
        z.push_back(var);
      } else if(var.getCausality() == INPUT){
        u.push_back(var);
      }
      break;
    default:
      casadi_assert_message(0,"Unknown category");
  }
}

std::string SymbolicOCP::qualifiedName(const XMLNode& nn){
  // Stringstream to assemble name
  stringstream qn;
  
  for(int i=0; i<nn.size(); ++i){
    // Add a dot
    if(i!=0) qn << ".";
    
    // Get the name part
    qn << nn[i].getAttribute("name");

    // Get the index, if any
    if(nn[i].size()>0){
      int ind;
      nn[i]["exp:ArraySubscripts"]["exp:IndexExpression"]["exp:IntegerLiteral"].getText(ind);
      qn << "[" << ind << "]";
    }
  }
  
  // Return the name
  return qn.str();
}

void SymbolicOCP::generateMuscodDatFile(const std::string& filename, const Dictionary& mc2_ops) const{
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
  vector<Variable> p;
  p.insert(p.end(),pi.begin(),pi.end());
  p.insert(p.end(),pf.begin(),pf.end());
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
  if(!x.empty()){
    datfile << "*  differential state start values, scale factors, and bounds" << endl;
    datfile << "sd(*,*)" << endl;
    for(int k=0; k<x.size(); ++k){
      datfile << k << ": " << x[k].getStart() << endl;
    }
    datfile << endl;
    
    datfile << "sd_sca(*,*)" << endl;
    for(int k=0; k<x.size(); ++k){
      datfile << k << ": " << x[k].getNominal() << endl;
    }
    datfile << endl;
    
    datfile << "sd_min(*,*)" << endl;
    for(int k=0; k<x.size(); ++k){
      datfile << k << ": " << x[k].getMin() << endl;
    }
    datfile << endl;
    
    datfile << "sd_max(*,*)" << endl;
    for(int k=0; k<x.size(); ++k){
      datfile << k << ": " << x[k].getMax() << endl;
    }
    datfile << endl;
    
    datfile << "sd_fix(*,*)" << endl;
    for(int k=0; k<x.size(); ++k){
      datfile << k << ": " << (x[k].getMin()==x[k].getMax()) << endl;
    }
    datfile << endl;

    datfile << "xd_name" << endl;
    for(int k=0; k<x.size(); ++k){
      datfile << k << ": " << x[k].getName() << endl;
    }
    datfile << endl;

    datfile << "xd_unit" << endl;
    for(int k=0; k<x.size(); ++k){
      datfile << k << ": " << x[k].getUnit() << endl;
    }
    datfile << endl;
  }
  
  // Algebraic state properties
  if(!z.empty()){
    datfile << "*  algebraic state start values, scale factors, and bounds" << endl;
    datfile << "sa(*,*)" << endl;
    for(int k=0; k<z.size(); ++k){
      datfile << k << ": " << z[k].getStart() << endl;
    }
    datfile << endl;
    
    datfile << "sa_sca(*,*)" << endl;
    for(int k=0; k<z.size(); ++k){
      datfile << k << ": " << z[k].getNominal() << endl;
    }
    datfile << endl;
    
    datfile << "sa_min(*,*)" << endl;
    for(int k=0; k<z.size(); ++k){
      datfile << k << ": " << z[k].getMin() << endl;
    }
    datfile << endl;
    
    datfile << "sa_max(*,*)" << endl;
    for(int k=0; k<z.size(); ++k){
      datfile << k << ": " << z[k].getMax() << endl;
    }
    datfile << endl;
    
    datfile << "sa_fix(*,*)" << endl;
    for(int k=0; k<z.size(); ++k){
      datfile << k << ": " << (z[k].getMin()==z[k].getMax()) << endl;
    }
    datfile << endl;

    datfile << "xa_name" << endl;
    for(int k=0; k<z.size(); ++k){
      datfile << k << ": " << z[k].getName() << endl;
    }
    datfile << endl;
    
    datfile << "xa_unit" << endl;
    for(int k=0; k<z.size(); ++k){
      datfile << k << ": " << z[k].getUnit() << endl;
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

std::vector<Variable>& SymbolicOCP::variableByType(SymbolicOCPVariables type){
  switch(type){
    case VAR_X: return x;
    case VAR_Z: return z;
    case VAR_Q: return q;
    case VAR_CI: return ci;
    case VAR_CD: return cd;
    case VAR_PI: return pi;
    case VAR_PD: return pd;
    case VAR_PF: return pf;
    case VAR_Y: return y;
    case VAR_U: return u;
    case NUM_VAR: casadi_error("SymbolicOCPVariable " << type << " out of range."); return x; // avoid -Wreturn-type
  }
}
    
const std::vector<Variable>& SymbolicOCP::variableByType(SymbolicOCPVariables type) const{
  return const_cast<SymbolicOCP*>(this)->variableByType(type);
}


} // namespace CasADi

