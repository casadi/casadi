#include "optimica_ocp.hpp"
#include <algorithm>
#include <set>
#include <cassert>
#include "../casadi/expression_tools.hpp"
#include "../casadi/casadi_exception.hpp"
#include "../casadi/stl_vector_tools.hpp"

using namespace std;
namespace CasADi{
  namespace Modelica{

OCP::OCP(){
  variables = Variable("variables");
}

void OCP::print(ostream &stream) const{
  // Variables sorted by type
  stream << "Variables" << endl;
  stream << variables << endl;

  // Print the initial equations
  stream << "Dynamic equations" << endl;
  for(vector<SX>::const_iterator it=dyneq.begin(); it!=dyneq.end(); it++){
    stream << "0 == "<< *it << endl;
  }
  stream << endl;

  stream << "Time state(s):                  " << t << endl;
  stream << "Differential states (implicit): " << x << endl;
  stream << "State derivatives:              " << xdot << endl;
  stream << "Differential states (explicit): " << xd << endl;
  stream << "Algebraic states:               " << xa << endl;
  stream << "Controls:                       " << u << endl;
  stream << "Parameter:                      " << p << endl;
  stream << "Dependent variables:            " << d << endl;
  stream << endl;
  
  // Print the equations
  stream << "Differential-algebraic equations (implicit)" << endl;
  for(vector<SX>::const_iterator it=dyneq.begin(); it!=dyneq.end(); it++)
    stream << "0 == " << *it << endl;
  stream << endl;
  
  stream << "Differential equations (explicit)" << endl;
  for(int i=0; i<diffeq.size(); ++i)
    stream << "der{" << xd[i] << "}" << " == " << diffeq[i] << endl;
  stream << endl;
  
  stream << "Algebraic equations" << endl;
  for(int i=0; i<algeq.size(); ++i)
    stream << "0 == " << algeq[i] << endl;
  stream << endl;
  
  stream << "Initial equations" << endl;
  for(vector<SX>::const_iterator it=initeq.begin(); it!=initeq.end(); it++){
    stream << "0 == " << *it << endl;
  }
  stream << endl;
  
  // Dependent equations
  stream << "Dependent equations" << endl;
  for(int i=0; i<d.size(); ++i)
    stream << d[i] << " == " << depdef[i] << endl;
  stream << endl;

  // Mayer terms
  stream << "Mayer objective terms" << endl;
  for(int i=0; i<mterm.size(); ++i)
    stream << mterm[i] << " at time == " << mtp[i] << endl;
  stream << endl;
  
  // Constraint functions
  stream << "Constraint functions" << endl;
  for(int i=0; i<cfcn.size(); ++i)
    stream << cfcn_lb[i] << " <= " << cfcn[i] << " <= " << cfcn_ub[i] << endl;
  stream << endl;
  
  // Constraint functions
  stream << "Time horizon" << endl;
  stream << "t0 = " << t0 << endl;
  stream << "tf = " << tf << endl;
  
}

void OCP::sortVariables(){
  // Clear the existing variables
  t.clear();
  x.clear();
  xdot.clear();
  xd.clear();
  xa.clear();
  u.clear();
  p.clear();
  d.clear();

  // Get all the variables
  vector<Variable> v = variables;

  cout << "variables  111  = " << variables << endl;
  
  // Loop over variables
  for(vector<Variable>::iterator it=v.begin(); it!=v.end(); ++it){
    // Make sure that the variable is initialized
    switch((*it)->type){
      case TYPE_TIME:               t.push_back(*it);  break;
      case TYPE_STATE:              x.push_back(it->sx());  xdot.push_back(it->der()); break;
      case TYPE_ALGEBRAIC:          xa.push_back(*it); break;
      case TYPE_CONTROL:            u.push_back(*it);  break;
      case TYPE_PARAMETER:          p.push_back(*it);  break;
      case TYPE_DEPENDENT:          d.push_back(*it);  break;
      default:
        throw CasadiException("OCP::sortVariables: unknown type");
    }
  }
  
  // Collect binding equations
  depdef.clear();
  depdef.reserve(d.size());
  for(vector<Variable>::iterator it=d.begin(); it!=d.end(); ++it)
    depdef.push_back(it->sx());
}

void OCP::makeExplicit(){
  // Dynamic equation
  SXMatrix dae(dyneq);
  SXMatrix xdot1(xdot);
  
  // Take the Jacobian of the ode with respect to xdot
  SXMatrix J = jacobian(dae,xdot1);

  // Make sure that J is invertable
  if(dependsOn(J,xdot1)) throw CasadiException("OCP::makeExplicit:: Dynamic equation not affine in state derivative");
  if(isZero(det(J))) throw CasadiException("OCP::makeExplicit: Jacobian not invertable");

  // Write the differential equation in explicit form
  SXMatrix rhs = solve(J,prod(J,xdot1-dae));

  // Simplify the expression
  simplify(rhs);
  
  // Save as explicit state
  xd = x;
  x.clear();
  xdot.clear();
  diffeq.resize(rhs.numel());
  for(int i=0; i<diffeq.size(); ++i)
    diffeq[i] = rhs.getElement(i);
  dyneq.clear();
}


void OCP::makeSemiExplicit(){
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
}

  

  } // namespace Modelica
} // namespace CasADi

