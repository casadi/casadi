/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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

#include <map>
#include <string>
#include <sstream>
#include <ctime>

#include "../std_vector_tools.hpp"
#include "../casadi_exception.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../function/integrator.hpp"
#include "../casadi_calculus.hpp"
#include "xml_file.hpp"

using namespace std;
namespace casadi {

  SymbolicOCP::SymbolicOCP(bool ignore_timed_variables)
      : ignore_timed_variables_(ignore_timed_variables) {
    t = SX::sym("t");
    t0 = t0_guess = numeric_limits<double>::quiet_NaN();
    tf = tf_guess = numeric_limits<double>::quiet_NaN();
    t0_free = false;
    tf_free = false;

    // Start with vectors of zero length
    this->s=this->x=this->z=this->q=this->ci=this->cd=this->pi=this->pd=this->p=this->y=
        this->u=this->path=this->point = SX::zeros(0, 1);
  }

  void SymbolicOCP::parseFMI(const std::string& filename) {

    // Load
    XmlFile xml_file("tinyxml");
    XmlNode document = xml_file.parse(filename);

    // **** Add model variables ****
    {
      //if (verbose) cout << "Adding model variables." << endl;

      // Get a reference to the ModelVariables node
      const XmlNode& modvars = document[0]["ModelVariables"];

      // Add variables
      for (int i=0; i<modvars.size(); ++i) {

        // Get a reference to the variable
        const XmlNode& vnode = modvars[i];

        // Get the attributes
        string name        = vnode.getAttribute("name");
        int valueReference;
        vnode.readAttribute("valueReference", valueReference);
        string variability = vnode.getAttribute("variability");
        string causality   = vnode.getAttribute("causality");
        string alias       = vnode.getAttribute("alias");

        // Skip to the next variable if its an alias
        if (alias.compare("alias") == 0 || alias.compare("negatedAlias") == 0)
          continue;

        // Get the name
        const XmlNode& nn = vnode["QualifiedName"];
        string qn = qualifiedName(nn);

        // Add variable, if not already added
        if (varmap_.find(qn)==varmap_.end()) {

          // Create variable
          Variable var;
          var.setName(name);

          // Value reference
          var.valueReference = valueReference;

          // Variability
          if (variability.compare("constant")==0)
            var.variability = CONSTANT;
          else if (variability.compare("parameter")==0)
            var.variability = PARAMETER;
          else if (variability.compare("discrete")==0)
            var.variability = DISCRETE;
          else if (variability.compare("continuous")==0)
            var.variability = CONTINUOUS;
          else
            throw CasadiException("Unknown variability");

          // ODE is zero for non-continuous variables
          if (var.variability!=CONTINUOUS) var.ode = 0;

          // Causality
          if (causality.compare("input")==0)
            var.causality = INPUT;
          else if (causality.compare("output")==0)
            var.causality = OUTPUT;
          else if (causality.compare("internal")==0)
            var.causality = INTERNAL;
          else
            throw CasadiException("Unknown causality");

          // Alias
          if (alias.compare("noAlias")==0)
            var.alias = NO_ALIAS;
          else if (alias.compare("alias")==0)
            var.alias = ALIAS;
          else if (alias.compare("negatedAlias")==0)
            var.alias = NEGATED_ALIAS;
          else
            throw CasadiException("Unknown alias");

          // Other properties
          if (vnode.hasChild("Real")) {
            const XmlNode& props = vnode["Real"];
            props.readAttribute("unit", var.unit, false);
            props.readAttribute("displayUnit", var.displayUnit, false);
            double dmin = -numeric_limits<double>::infinity();
            props.readAttribute("min", dmin, false);
            var.min = dmin;
            double dmax =  numeric_limits<double>::infinity();
            props.readAttribute("max", dmax, false);
            var.max = dmax;
            double dinitialGuess = 0;
            props.readAttribute("initialGuess", dinitialGuess, false);
            var.initialGuess = dinitialGuess;
            props.readAttribute("start", var.start, false);
            props.readAttribute("nominal", var.nominal, false);
            props.readAttribute("free", var.free, false);
          }

          // Variable category
          if (vnode.hasChild("VariableCategory")) {
            string cat = vnode["VariableCategory"].getText();
            if (cat.compare("derivative")==0)
              var.category = CAT_DERIVATIVE;
            else if (cat.compare("state")==0)
              var.category = CAT_STATE;
            else if (cat.compare("dependentConstant")==0)
              var.category = CAT_DEPENDENT_CONSTANT;
            else if (cat.compare("independentConstant")==0)
              var.category = CAT_INDEPENDENT_CONSTANT;
            else if (cat.compare("dependentParameter")==0)
              var.category = CAT_DEPENDENT_PARAMETER;
            else if (cat.compare("independentParameter")==0)
              var.category = CAT_INDEPENDENT_PARAMETER;
            else if (cat.compare("algebraic")==0)
              var.category = CAT_ALGEBRAIC;
            else
              throw CasadiException("Unknown variable category: " + cat);
          }

          // Add to list of variables
          addVariable(qn, var);

          // Sort expression
          switch (var.category) {
          case CAT_DERIVATIVE:
            // Skip - meta information about time derivatives is
            //        kept together with its parent variable
            break;
          case CAT_STATE:
            this->s.append(var.v);
            break;
          case CAT_DEPENDENT_CONSTANT:
            this->cd.append(var.v);
            break;
          case CAT_INDEPENDENT_CONSTANT:
            this->ci.append(var.v);
            break;
          case CAT_DEPENDENT_PARAMETER:
            this->pd.append(var.v);
            break;
          case CAT_INDEPENDENT_PARAMETER:
            if (var.free) {
              this->p.append(var.v);
            } else {
              this->pi.append(var.v);
            }
            break;
          case CAT_ALGEBRAIC:
            if (var.causality == INTERNAL) {
              this->s.append(var.v);
            } else if (var.causality == INPUT) {
              this->u.append(var.v);
            }
            break;
          default:
            casadi_error("Unknown category");
          }
        }
      }
    }

    // **** Add binding equations ****
    {
      //if (verbose) cout << "Adding binding equations." << endl;

      // Get a reference to the BindingEquations node
      const XmlNode& bindeqs = document[0]["equ:BindingEquations"];

      for (int i=0; i<bindeqs.size(); ++i) {
        const XmlNode& beq = bindeqs[i];

        // Get the variable and binding expression
        Variable& var = readVariable(beq[0]);
        SX bexpr = readExpr(beq[1][0]);
        setBeq(var.v, bexpr.toScalar());
      }
    }

    // **** Add dynamic equations ****
    {
      // Get a reference to the DynamicEquations node
      const XmlNode& dyneqs = document[0]["equ:DynamicEquations"];

      // Add equations
      for (int i=0; i<dyneqs.size(); ++i) {

        // Get a reference to the variable
        const XmlNode& dnode = dyneqs[i];

        // Add the differential equation
        SX de_new = readExpr(dnode[0]);
        dae.append(de_new);
      }
    }

    // **** Add initial equations ****
    {
      // Get a reference to the DynamicEquations node
      const XmlNode& initeqs = document[0]["equ:InitialEquations"];

      // Add equations
      for (int i=0; i<initeqs.size(); ++i) {

        // Get a reference to the node
        const XmlNode& inode = initeqs[i];

        // Add the differential equations
        for (int i=0; i<inode.size(); ++i) {
          initial.append(readExpr(inode[i]));
        }
      }
    }

    // **** Add optimization ****
    if (document[0].hasChild("opt:Optimization")) {

      // Get a reference to the DynamicEquations node
      const XmlNode& opts = document[0]["opt:Optimization"];

      // Start time
      const XmlNode& intervalStartTime = opts["opt:IntervalStartTime"];
      if (intervalStartTime.hasChild("opt:Value"))
        intervalStartTime["opt:Value"].getText(t0);
      if (intervalStartTime.hasChild("opt:Free"))
        intervalStartTime["opt:Free"].getText(t0_free);
      if (intervalStartTime.hasChild("opt:InitialGuess"))
        intervalStartTime["opt:InitialGuess"].getText(t0_guess);

      // Terminal time
      const XmlNode& IntervalFinalTime = opts["opt:IntervalFinalTime"];
      if (IntervalFinalTime.hasChild("opt:Value"))
        IntervalFinalTime["opt:Value"].getText(tf);
      if (IntervalFinalTime.hasChild("opt:Free"))
        IntervalFinalTime["opt:Free"].getText(tf_free);
      if (IntervalFinalTime.hasChild("opt:InitialGuess"))
        IntervalFinalTime["opt:InitialGuess"].getText(tf_guess);

      // Time points
      const XmlNode& tpnode = opts["opt:TimePoints"];
      tp.resize(tpnode.size());
      for (int i=0; i<tp.size(); ++i) {
        // Get index
        int index;
        tpnode[i].readAttribute("index", index);

        // Get value
        double value;
        tpnode[i].readAttribute("value", value);
        tp[i] = value;

        if (!ignore_timed_variables_) {
          // Allocate all the timed variables
          for (int k=0; k<tpnode[i].size(); ++k) {
            string qn = qualifiedName(tpnode[i][k]);
            atTime(qn, value, true);
          }
        }
      }

      for (int i=0; i<opts.size(); ++i) {

        // Get a reference to the node
        const XmlNode& onode = opts[i];

        // Get the type
        if (onode.checkName("opt:ObjectiveFunction")) { // mayer term
          try {
            // Add components
            for (int i=0; i<onode.size(); ++i) {
              const XmlNode& var = onode[i];

              // If string literal, ignore
              if (var.checkName("exp:StringLiteral"))
                continue;

              // Read expression
              SX v = readExpr(var);
              mterm.append(v);
            }
          } catch(exception& ex) {
            throw CasadiException(std::string("addObjectiveFunction failed: ") + ex.what());
          }
        } else if (onode.checkName("opt:IntegrandObjectiveFunction")) {
          try {
            for (int i=0; i<onode.size(); ++i) {
              const XmlNode& var = onode[i];

              // If string literal, ignore
              if (var.checkName("exp:StringLiteral"))
                continue;

              // Read expression
              SX v = readExpr(var);
              lterm.append(v);
            }
          } catch(exception& ex) {
            throw CasadiException(std::string("addIntegrandObjectiveFunction failed: ")
                                  + ex.what());
          }
        } else if (onode.checkName("opt:IntervalStartTime")) {
          // Ignore, treated above
        } else if (onode.checkName("opt:IntervalFinalTime")) {
          // Ignore, treated above
        } else if (onode.checkName("opt:TimePoints")) {
          // Ignore, treated above
        } else if (onode.checkName("opt:PointConstraints")) {
          for (int i=0; i<onode.size(); ++i) {
            const XmlNode& constr_i = onode[i];

            // Create a new variable
            Variable v;
            stringstream ss;
            ss << "point_" << i;
            v.setName(ss.str());

            // Get the definition and bounds
            if (constr_i.checkName("opt:ConstraintLeq")) {
              v.beq = readExpr(constr_i[0]).toScalar();
              v.max = readExpr(constr_i[1]).toScalar();
            } else if (constr_i.checkName("opt:ConstraintGeq")) {
              v.beq = readExpr(constr_i[0]).toScalar();
              v.min = readExpr(constr_i[1]).toScalar();
            } else if (constr_i.checkName("opt:ConstraintEq")) {
              v.beq = readExpr(constr_i[0]).toScalar();
              v.max = v.min = readExpr(constr_i[1]).toScalar();
            } else {
              cerr << "unknown constraint type" << constr_i.getName() << endl;
              throw CasadiException("SymbolicOCP::addConstraints");
            }

            // Add to list of variables and outputs
            addVariable(ss.str(), v);
            this->point.append(v.v);
          }
        } else if (onode.checkName("opt:Constraints") || onode.checkName("opt:PathConstraints")) {
          for (int i=0; i<onode.size(); ++i) {
            const XmlNode& constr_i = onode[i];
            // Create a new variable
            Variable v;
            stringstream ss;
            ss << "path_" << i;
            v.setName(ss.str());

            // Get the definition and bounds
            if (constr_i.checkName("opt:ConstraintLeq")) {
              v.beq = readExpr(constr_i[0]).toScalar();
              v.max = readExpr(constr_i[1]).toScalar();
            } else if (constr_i.checkName("opt:ConstraintGeq")) {
              v.beq = readExpr(constr_i[0]).toScalar();
              v.min = readExpr(constr_i[1]).toScalar();
            } else if (constr_i.checkName("opt:ConstraintEq")) {
              v.beq = readExpr(constr_i[0]).toScalar();
              v.max = v.min = readExpr(constr_i[1]).toScalar();
            } else {
              cerr << "unknown constraint type" << constr_i.getName() << endl;
              throw CasadiException("SymbolicOCP::addConstraints");
            }

            // Add to list of variables and outputs
            addVariable(ss.str(), v);
            this->path.append(v.v);
          }
        } else {
          throw CasadiException(string("SymbolicOCP::addOptimization: Unknown node ")
                                +onode.getName());
        }
      }
    }

    // Make sure that the dimensions are consistent at this point
    casadi_assert_warning(this->s.size()==this->dae.size(),
                          "The number of differential-algebraic equations does not match "
                          "the number of implicitly defined states.");
    casadi_assert_warning(this->z.size()==this->alg.size(),
                          "The number of algebraic equations (equations not involving "
                          "differentiated variables) does not match the number of "
                          "algebraic variables.");
  }

  Variable& SymbolicOCP::readVariable(const XmlNode& node) {
    // Qualified name
    string qn = qualifiedName(node);

    // Find and return the variable
    return variable(qn);
  }

  SX SymbolicOCP::readExpr(const XmlNode& node) {
    const string& fullname = node.getName();
    if (fullname.find("exp:")== string::npos) {
      casadi_error("SymbolicOCP::readExpr: unknown - expression is supposed to "
                   "start with 'exp:' , got " << fullname);
    }

    // Chop the 'exp:'
    string name = fullname.substr(4);

    // The switch below is alphabetical, and can be thus made more efficient,
    // for example by using a switch statement of the first three letters,
    // if it would ever become a bottleneck
    if (name.compare("Add")==0) {
      return readExpr(node[0]) + readExpr(node[1]);
    } else if (name.compare("Acos")==0) {
      return acos(readExpr(node[0]));
    } else if (name.compare("Asin")==0) {
      return asin(readExpr(node[0]));
    } else if (name.compare("Atan")==0) {
      return atan(readExpr(node[0]));
    } else if (name.compare("Cos")==0) {
      return cos(readExpr(node[0]));
    } else if (name.compare("Der")==0) {
      const Variable& v = readVariable(node[0]);
      return v.d;
    } else if (name.compare("Div")==0) {
      return readExpr(node[0]) / readExpr(node[1]);
    } else if (name.compare("Exp")==0) {
      return exp(readExpr(node[0]));
    } else if (name.compare("Identifier")==0) {
      return readVariable(node).v;
    } else if (name.compare("IntegerLiteral")==0) {
      int val;
      node.getText(val);
      return val;
    } else if (name.compare("Instant")==0) {
      double val;
      node.getText(val);
      return val;
    } else if (name.compare("Log")==0) {
      return log(readExpr(node[0]));
    } else if (name.compare("LogLt")==0) { // Logical less than
      return readExpr(node[0]) < readExpr(node[1]);
    } else if (name.compare("LogGt")==0) { // Logical less than
      return readExpr(node[0]) > readExpr(node[1]);
    } else if (name.compare("Mul")==0) { // Multiplication
      return readExpr(node[0]) * readExpr(node[1]);
    } else if (name.compare("Neg")==0) {
      return -readExpr(node[0]);
    } else if (name.compare("NoEvent")==0) {
      // NOTE: This is a workaround, we assume that whenever NoEvent occurs,
      // what is meant is a switch
      int n = node.size();

      // Default-expression
      SX ex = readExpr(node[n-1]);

      // Evaluate ifs
      for (int i=n-3; i>=0; i -= 2) ex = if_else(readExpr(node[i]), readExpr(node[i+1]), ex);

      return ex;
    } else if (name.compare("Pow")==0) {
      return pow(readExpr(node[0]), readExpr(node[1]));
    } else if (name.compare("RealLiteral")==0) {
      double val;
      node.getText(val);
      return val;
    } else if (name.compare("Sin")==0) {
      return sin(readExpr(node[0]));
    } else if (name.compare("Sqrt")==0) {
      return sqrt(readExpr(node[0]));
    } else if (name.compare("StringLiteral")==0) {
      throw CasadiException(node.getText());
    } else if (name.compare("Sub")==0) {
      return readExpr(node[0]) - readExpr(node[1]);
    } else if (name.compare("Tan")==0) {
      return tan(readExpr(node[0]));
    } else if (name.compare("Time")==0) {
      return t.toScalar();
    } else if (name.compare("TimedVariable")==0) {
      if (ignore_timed_variables_) {
        return readVariable(node[0]).v;
      } else {
        // Get the index of the time point
        int index;
        node.readAttribute("timePointIndex", index);
        return readVariable(node[0]).atTime(tp[index]);
      }
    }

    // throw error if reached this point
    throw CasadiException(string("SymbolicOCP::readExpr: Unknown node: ") + name);

  }

  void SymbolicOCP::repr(std::ostream &stream, bool trailing_newline) const {
    stream << "Flat OCP";
    if (trailing_newline) stream << endl;
  }

  void SymbolicOCP::print(ostream &stream, bool trailing_newline) const {
    stream << "Dimensions: ";
    stream << "#s = " << this->s.size() << ", ";
    stream << "#x = " << this->x.size() << ", ";
    stream << "#z = " << this->z.size() << ", ";
    stream << "#q = " << this->q.size() << ", ";
    stream << "#y = " << this->y.size() << ", ";
    stream << "#pi = " << this->pi.size() << ", ";
    stream << "#pd = " << this->pd.size() << ", ";
    stream << "#pf = " << this->p.size() << ", ";
    stream << "#ci =  " << this->ci.size() << ", ";
    stream << "#cd =  " << this->cd.size() << ", ";
    stream << "#u = " << this->u.size() << ", ";
    stream << endl << endl;

    // Variables in the class hierarchy
    stream << "Variables" << endl;

    // Print the variables
    stream << "{" << endl;
    stream << "  t = " << str(this->t) << endl;
    stream << "  s = " << str(this->s) << endl;
    stream << "  x = " << str(this->x) << endl;
    stream << "  z =  " << str(this->z) << endl;
    stream << "  q =  " << str(this->q) << endl;
    stream << "  y =  " << str(this->y) << endl;
    stream << "  pi =  " << str(this->pi) << endl;
    stream << "  pd =  " << str(this->pd) << endl;
    stream << "  pf =  " << str(this->p) << endl;
    stream << "  ci =  " << str(this->ci) << endl;
    stream << "  cd =  " << str(this->cd) << endl;
    stream << "  u =  " << str(this->u) << endl;
    stream << "}" << endl;

    if (!this->dae.isEmpty()) {
      stream << "Fully-implicit differential-algebraic equations" << endl;
      for (int k=0; k<this->dae.size(); ++k) {
        stream << "0 == " << this->dae.at(k) << endl;
      }
      stream << endl;
    }

    if (!this->x.isEmpty()) {
      stream << "Differential equations" << endl;
      for (int k=0; k<this->x.size(); ++k) {
        stream << str(der(this->x[k])) << " == " << str(ode(this->x[k])) << endl;
      }
      stream << endl;
    }

    if (!this->alg.isEmpty()) {
      stream << "Algebraic equations" << endl;
      for (int k=0; k<this->z.size(); ++k) {
        stream << "0 == " << str(this->alg[k]) << endl;
      }
      stream << endl;
    }

    if (!this->q.isEmpty()) {
      stream << "Quadrature equations" << endl;
      for (int k=0; k<this->q.size(); ++k) {
        stream << str(der(this->q[k])) << " == " << str(ode(this->q[k])) << endl;
      }
      stream << endl;
    }

    if (!this->initial.isEmpty()) {
      stream << "Initial equations" << endl;
      for (SX::const_iterator it=this->initial.begin(); it!=this->initial.end(); it++) {
        stream << "0 == " << *it << endl;
      }
      stream << endl;
    }

    if (!this->pi.isEmpty()) {
      stream << "Independent parameters" << endl;
      for (int i=0; i<this->pi.size(); ++i)
        stream << this->pi.at(i) << " == " << str(beq(this->pi.at(i))) << endl;
      stream << endl;
    }

    if (!this->pd.isEmpty()) {
      stream << "Dependent parameters" << endl;
      for (int i=0; i<this->pd.size(); ++i)
        stream << this->pd.at(i) << " == " << str(beq(this->pd.at(i))) << endl;
      stream << endl;
    }

    if (!this->ci.isEmpty()) {
      stream << "Independent constants" << endl;
      for (int i=0; i<this->ci.size(); ++i)
        stream << this->ci.at(i) << " == " << str(beq(this->ci.at(i))) << endl;
      stream << endl;
    }

    if (!this->cd.isEmpty()) {
      stream << "Dependent constants" << endl;
      for (int i=0; i<this->cd.size(); ++i)
        stream << this->cd.at(i) << " == " << str(beq(this->cd.at(i))) << endl;
      stream << endl;
    }

    if (!this->y.isEmpty()) {
      stream << "Output variables" << endl;
      for (int i=0; i<this->y.size(); ++i)
        stream << this->y.at(i) << " == " << str(beq(this->y.at(i))) << endl;
      stream << endl;
    }

    if (!this->mterm.isEmpty()) {
      stream << "Mayer objective terms" << endl;
      for (int i=0; i<this->mterm.size(); ++i)
        stream << this->mterm.at(i) << endl;
      stream << endl;
    }

    if (!this->lterm.isEmpty()) {
      stream << "Lagrange objective terms" << endl;
      for (int i=0; i<this->lterm.size(); ++i)
        stream << this->lterm.at(i) << endl;
      stream << endl;
    }

    if (!this->path.isEmpty()) {
      stream << "Path constraints" << endl;
      for (int i=0; i<this->path.size(); ++i)
        stream << str(max(this->path.at(i))) << " <= "
               << this->path.at(i) << " <= " << str(max(this->path.at(i))) << endl;
      stream << endl;
    }

    if (!this->point.isEmpty()) {
      stream << "Point constraints" << endl;
      for (int i=0; i<this->point.size(); ++i)
        stream << str(max(this->point.at(i))) << " <= "
               << this->point.at(i) << " <= " << str(max(this->point.at(i))) << endl;
      stream << endl;
    }

    // Constraint functions
    stream << "Time horizon" << endl;
    stream << "t0 = " << this->t0 << endl;
    stream << "tf = " << this->tf << endl;
    stream << "tp = " << this->tp << endl;
    if (trailing_newline) stream << endl;
  }

  void SymbolicOCP::eliminateLagrangeTerms() {
    // Index for the names
    int ind = 0;
    // For every integral term in the objective function
    for (SX::iterator it=this->lterm.begin(); it!=this->lterm.end(); ++it) {

      // Give a name to the quadrature state
      stringstream q_name;
      q_name << "q_" << ind++;

      // Create a new quadrature state
      Variable qv;
      qv.setName(q_name.str());

      // Set attributes
      qv.variability = CONTINUOUS;
      qv.causality = INTERNAL;
      qv.start = 0.0;
      if (tf==tf) qv.nominal = this->tf; // if not not-a-number

      // Add to the list of variables
      addVariable(q_name.str(), qv);

      // Add to the quadrature states
      this->q.append(qv.v);

      // Add the Lagrange term to the list of quadratures
      setOde(qv.v, *it);

      // Add to the list of Mayer terms
      this->mterm.append(qv.v);
    }

    // Remove the Lagrange terms
    this->lterm.clear();
  }

  void SymbolicOCP::eliminateQuadratureStates() {
    // Move all the quadratures to the list of differential states
    this->x.append(this->q);
    this->q = SX::zeros(0, 1);
  }

  void SymbolicOCP::scaleVariables() {
    // Gather variables and expressions to replace
    vector<SXElement> v_id, v_rep, ex_rep;
    for (VarMap::iterator it=varmap_.begin(); it!=varmap_.end(); ++it) {
      if (it->second.nominal!=1) {
        Variable& v=it->second;
        casadi_assert(v.nominal!=0);
        v.beq *= v.nominal;
        v.ode /= v.nominal;
        v.min /= v.nominal;
        v.max /= v.nominal;
        v.start /= v.nominal;
        v.derivativeStart /= v.nominal;
        v.initialGuess /= v.nominal;
        v_id.push_back(v.v);
        v_rep.push_back(v.beq);
        v_id.push_back(v.d);
        v_rep.push_back(v.d * v.nominal);
        ex_rep.push_back(v.ode);
        ex_rep.push_back(v.beq);
      }
    }

    // Quick return if no expressions to substitute
    if (v_id.empty()) return;

    // Collect all expressions to be replaced
    vector<SX> ex;
    ex.push_back(ex_rep);
    ex.push_back(this->dae);
    ex.push_back(this->alg);
    ex.push_back(this->initial);
    ex.push_back(this->mterm);
    ex.push_back(this->lterm);

    // Substitute all at once (more efficient since they may have common subexpressions)
    ex = substitute(ex, vector<SX>(1, v_id), vector<SX>(1, v_rep));

    // Get the modified expressions
    vector<SX>::const_iterator it=ex.begin();
    ex_rep = (*it++).data();
    this->dae = *it++;
    this->alg = *it++;
    this->initial = *it++;
    this->mterm = *it++;
    this->lterm = *it++;
    casadi_assert(it==ex.end());

    // Save the substituted expressions
    vector<SXElement>::iterator ex_rep_it = ex_rep.begin();
    for (VarMap::iterator it=varmap_.begin(); it!=varmap_.end(); ++it) {
      Variable& v=it->second;
      if (v.nominal!=1) {
        v.ode = *ex_rep_it++;
        v.beq = *ex_rep_it++;
        v.nominal=1;
      }
    }
    casadi_assert(ex_rep_it==ex_rep.end());
  }

  void SymbolicOCP::eliminateIndependentParameters() {
    // Collect all expressions to be replaced
    vector<SX> ex;
    ex.push_back(this->dae);
    ex.push_back(ode(this->x));
    ex.push_back(this->alg);
    ex.push_back(ode(this->q));
    ex.push_back(beq(this->y));
    ex.push_back(this->initial);
    ex.push_back(this->mterm);
    ex.push_back(this->lterm);
    ex.push_back(beq(this->pd));

    // Substitute all at once (since they may have common subexpressions)
    ex = substitute(ex, vector<SX>(1, this->pi), vector<SX>(1, beq(this->pi)));

    // Get the modified expressions
    vector<SX>::const_iterator it=ex.begin();
    this->dae = *it++;
    setOde(this->x, *it++);
    this->alg = *it++;
    setOde(this->q, *it++);
    setBeq(this->y, *it++);
    this->initial = *it++;
    this->mterm = *it++;
    this->lterm = *it++;
    setBeq(this->pd, *it++);
    casadi_assert(it==ex.end());
  }

  void SymbolicOCP::sortDependentParameters() {
    // Quick return if no dependent parameters
    if (this->pd.isEmpty()) return;

    // Find out which dependent parameter depends on which binding equation
    SXFunction f(this->pd, this->pd - beq(this->pd));
    f.init();
    Sparsity sp = f.jacSparsity();

    // BLT transformation
    vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    sp.dulmageMendelsohn(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);

    // Permute variables
    this->pd = this->pd(colperm);
  }

  void SymbolicOCP::eliminateDependentParameterInterdependencies() {
    // Quick return if no dependent parameters
    if (this->pd.isEmpty()) return;

    // Begin by sorting the parameters
    sortDependentParameters();

    // Sort the equations by causality
    SX pd_def = beq(this->pd);
    substituteInPlace(this->pd, pd_def, false);

    // Make sure that the dependent variables have been properly eliminated from the definitions
    casadi_assert(!dependsOn(pd_def, this->pd));

    // Save new binding equations
    setBeq(this->pd, pd_def);
  }

  void SymbolicOCP::eliminateDependentParameters() {
    // Quick return if no dependent parameters
    if (this->pd.isEmpty()) return;

    // Remove interdependencies
    eliminateDependentParameterInterdependencies();

    // Collect all expressions to be replaced
    vector<SX> ex;
    ex.push_back(this->dae);
    ex.push_back(ode(this->x));
    ex.push_back(this->alg);
    ex.push_back(ode(this->q));
    ex.push_back(beq(this->y));
    ex.push_back(this->initial);
    ex.push_back(this->mterm);
    ex.push_back(this->lterm);

    // Substitute all at once (since they may have common subexpressions)
    ex = substitute(ex, vector<SX>(1, this->pd), vector<SX>(1, beq(this->pd)));

    // Get the modified expressions
    vector<SX>::const_iterator it=ex.begin();
    this->dae = *it++;
    setOde(this->x, *it++);
    this->alg = *it++;
    setOde(this->q, *it++);
    setBeq(this->y, *it++);
    this->initial = *it++;
    this->mterm = *it++;
    this->lterm = *it++;
    casadi_assert(it==ex.end());
  }

  void SymbolicOCP::sortOutputs() {
    // Quick return if no outputs
    if (this->y.isEmpty()) return;

    // Find out which dependent parameter depends on which binding equation
    SXFunction f(this->y, this->y - beq(this->y));
    f.init();
    Sparsity sp = f.jacSparsity();

    // BLT transformation
    vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    sp.dulmageMendelsohn(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);

    // Permute variables
    this->y = this->y(colperm);
  }

  void SymbolicOCP::eliminateOutputInterdependencies() {
    // Quick return if no outputs
    if (this->y.isEmpty()) return;

    // Begin by sorting the outputs
    sortOutputs();

    // Sort the equations by causality
    SX y_def = beq(this->y);
    substituteInPlace(this->y, y_def, false);

    // Make sure that the outputs have been properly eliminated from the output definitions
    casadi_assert(!dependsOn(y_def, this->y));

    // Save new binding equations
    setBeq(this->y, y_def);
  }

  void SymbolicOCP::eliminateOutputs() {
    // Quick return if no dependent parameters
    if (this->y.isEmpty()) return;

    // Remove interdependencies
    eliminateOutputInterdependencies();

    // Collect all expressions to be replaced
    vector<SX> ex;
    ex.push_back(this->dae);
    ex.push_back(ode(this->x));
    ex.push_back(this->alg);
    ex.push_back(ode(this->q));
    ex.push_back(beq(this->y));
    ex.push_back(this->initial);
    ex.push_back(this->mterm);
    ex.push_back(this->lterm);

    // Substitute all at once (since they may have common subexpressions)
    ex = substitute(ex, vector<SX>(1, this->y), vector<SX>(1, beq(this->y)));

    // Get the modified expressions
    vector<SX>::const_iterator it=ex.begin();
    this->dae = *it++;
    setOde(this->x, *it++);
    this->alg = *it++;
    setOde(this->q, *it++);
    setBeq(this->y, *it++);
    this->initial = *it++;
    this->mterm = *it++;
    this->lterm = *it++;
    casadi_assert(it==ex.end());
  }

  void SymbolicOCP::scaleEquations() {
    casadi_error("SymbolicOCP::scaleEquations broken");

    cout << "Scaling equations ..." << endl;
    double time1 = clock();

    // Variables
    enum Variables {T, X, XDOT, Z, PI, PF, U, NUM_VAR};
    vector<SX > v(NUM_VAR); // all variables
    v[T] = this->t;
    v[X] = this->x;
    v[XDOT] = der(this->x); // BUG!!!
    v[Z] = this->z;
    v[PI] = this->pi;
    v[PF] = this->p;
    v[U] = this->u;

    // Create the jacobian of the implicit equations with respect to [x, z, p, u]
    SX xz;
    xz.append(v[X]);
    xz.append(v[Z]);
    xz.append(v[PI]);
    xz.append(v[PF]);
    xz.append(v[U]);
    SXFunction fcn = SXFunction(xz, ode(this->x));
    SXFunction J(v, fcn.jac());

    // Evaluate the Jacobian in the starting point
    J.init();
    J.setInput(0.0, T);
    J.setInput(start(this->x, true), X);
    J.input(XDOT).setAll(0.0);
    J.setInput(start(this->z, true), Z);
    J.setInput(start(this->pi, true), PI);
    J.setInput(start(this->p, true), PF);
    J.setInput(start(this->u, true), U);
    J.evaluate();

    // Get the maximum of every row
    Matrix<double> &J0 = J.output();
    vector<double> scale(J0.size1(), 0.0); // scaling factors
    for (int cc=0; cc<J0.size2(); ++cc) {
      // Loop over non-zero entries of the column
      for (int el=J0.colind(cc); el<J0.colind(cc+1); ++el) {
        // Row
        int rr=J0.row(el);

        // The scaling factor is the maximum norm, ignoring not-a-number entries
        if (!isnan(J0.at(el))) {
          scale[rr] = std::max(scale[rr], fabs(J0.at(el)));
        }
      }
    }

    // Make sure nonzero factor found
    for (int rr=0; rr<J0.size1(); ++rr) {
      if (scale[rr]==0) {
        cout << "Warning: Could not generate a scaling factor for equation " << rr;
        scale[rr]=1.;
      }
    }

    // Scale the equations
    setOde(this->x, ode(this->x)/scale);

    double time2 = clock();
    double dt = (time2-time1)/CLOCKS_PER_SEC;
    cout << "... equation scaling complete after " << dt << " seconds." << endl;
  }

  void SymbolicOCP::sortDAE() {
    // Quick return if no differential states
    if (this->x.isEmpty()) return;

    // Find out which differential equation depends on which differential state
    SXFunction f(der(this->s), this->dae);
    f.init();
    Sparsity sp = f.jacSparsity();

    // BLT transformation
    vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    sp.dulmageMendelsohn(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);

    // Permute equations
    this->dae = this->dae(rowperm);

    // Permute variables
    this->s = this->s(colperm);
  }

  void SymbolicOCP::sortALG() {
    // Quick return if no algebraic states
    if (this->z.isEmpty()) return;

    // Find out which algebraic equation depends on which algebraic state
    SXFunction f(this->z, this->alg);
    f.init();
    Sparsity sp = f.jacSparsity();

    // BLT transformation
    vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    sp.dulmageMendelsohn(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);

    // Permute equations
    this->alg = this->alg(rowperm);

    // Permute variables
    this->z = this->z(colperm);
  }

  void SymbolicOCP::makeSemiExplicit() {
    // Separate the algebraic variables and equations
    separateAlgebraic();

    // Quick return if there are no implicitly defined states
    if (this->s.isEmpty()) return;

    // Write the ODE as a function of the state derivatives
    SXFunction f(der(this->s), this->dae);
    f.init();

    // Get the sparsity of the Jacobian which can be used to determine which
    // variable can be calculated from which other
    Sparsity sp = f.jacSparsity();

    // BLT transformation
    vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    int nb = sp.dulmageMendelsohn(rowperm, colperm, rowblock, colblock,
                                  coarse_rowblock, coarse_colblock);

    // Permute equations
    this->dae = this->dae(rowperm);

    // Permute variables
    this->s = this->s(colperm);

    // Now write the sorted ODE as a function of the state derivatives
    f = SXFunction(der(this->s), this->dae);
    f.init();

    // Get the Jacobian
    SX J = f.jac();

    // Explicit ODE
    SX new_ode;

    // Loop over blocks
    for (int b=0; b<nb; ++b) {

      // Block size
      int bs = rowblock[b+1] - rowblock[b];

      // Get variables in the block
      SX xb = this->s(Slice(colblock[b], colblock[b+1]));

      // Get equations in the block
      SX fb = this->dae(Slice(rowblock[b], rowblock[b+1]));

      // Get local Jacobian
      SX Jb = J(Slice(rowblock[b], rowblock[b+1]), Slice(colblock[b], colblock[b+1]));

      // If Jb depends on xb, then the state derivative does not enter linearly
      // in the ODE and we cannot solve for the state derivative
      casadi_assert_message(!dependsOn(Jb, der(xb)),
                            "Cannot find an explicit expression for variable(s) " << xb);

      // Divide fb into a part which depends on vb and a part which doesn't according to
      // "fb == mul(Jb, vb) + fb_res"
      SX fb_res = substitute(fb, der(xb), SX::zeros(xb.sparsity()));
      SX fb_exp;

      // Solve for vb
      if (bs <= 3) {
        // Calculate inverse and multiply for very small matrices
        fb_exp = mul(inv(Jb), -fb_res);
      } else {
        // QR factorization
        fb_exp = solve(Jb, -fb_res);
      }

      // Add to explicitly determined equations and variables
      new_ode.append(fb_exp);
    }

    // Eliminate inter-dependencies
    substituteInPlace(der(this->s), new_ode, false);

    // Add to explicit differential states and ODE
    setOde(this->s, new_ode);
    this->x.append(this->s);
    this->dae = this->s = SX::zeros(0, 1);
  }

  void SymbolicOCP::eliminateAlgebraic() {
    // Quick return if there are no algebraic states
    if (this->z.isEmpty()) return;

    // Write the algebraic equations as a function of the algebraic states
    SXFunction f(this->z, alg);
    f.init();

    // Get the sparsity of the Jacobian which can be used to determine which
    // variable can be calculated from which other
    Sparsity sp = f.jacSparsity();

    // BLT transformation
    vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    int nb = sp.dulmageMendelsohn(rowperm, colperm, rowblock, colblock,
                                  coarse_rowblock, coarse_colblock);

    // Permute equations
    this->alg = this->alg(rowperm);

    // Permute variables
    this->z = this->z(colperm);

    // Rewrite the sorted algebraic equations as a function of the algebraic states
    f = SXFunction(this->z, this->alg);
    f.init();

    // Get the Jacobian
    SX J = f.jac();

    // Variables where we have found an explicit expression and where we haven't
    SX z_exp, z_imp;

    // Explicit and implicit equations
    SX f_exp, f_imp;

    // Loop over blocks
    for (int b=0; b<nb; ++b) {

      // Block size
      int bs = rowblock[b+1] - rowblock[b];

      // Get local variables
      SX zb = this->z(Slice(colblock[b], colblock[b+1]));

      // Get local equations
      SX fb = this->alg(Slice(rowblock[b], rowblock[b+1]));

      // Get local Jacobian
      SX Jb = J(Slice(rowblock[b], rowblock[b+1]), Slice(colblock[b], colblock[b+1]));

      // If Jb depends on zb, then we cannot (currently) solve for it explicitly
      if (dependsOn(Jb, zb)) {

        // Add the equations to the new list of algebraic equations
        f_imp.append(fb);

        // ... and the variables accordingly
        z_imp.append(zb);

      } else { // The variables that we wish to determine enter linearly

        // Divide fb into a part which depends on vb and a part which doesn't
        // according to "fb == mul(Jb, vb) + fb_res"
        SX fb_res = substitute(fb, zb, SX::zeros(zb.sparsity()));

        // Solve for vb
        SX fb_exp;
        if (bs <= 3) {
          // Calculate inverse and multiply for very small matrices
          fb_exp = mul(inv(Jb), -fb_res);
        } else {
          // QR factorization
          fb_exp = solve(Jb, -fb_res);
        }

        // Add to explicitly determined equations and variables
        z_exp.append(zb);
        f_exp.append(fb_exp);
      }
    }

    // Eliminate inter-dependencies in fb_exp
    substituteInPlace(z_exp, f_exp, false);

    // Add to the beginning of the dependent variables
    // (since the other dependent variable might depend on them)
    this->y = vertcat(z_exp, this->y);
    setBeq(z_exp, f_exp);

    // Save new algebraic equations
    this->z = z_imp;
    this->alg = f_imp;

    // Eliminate new dependent variables from the other equations
    eliminateOutputs();
  }

  void SymbolicOCP::makeExplicit() {
    // Start by transforming to semi-explicit form
    makeSemiExplicit();

    // Then eliminate the algebraic variables
    eliminateAlgebraic();

    // Error if still algebraic variables
    casadi_assert_message(this->z.isEmpty(), "Failed to eliminate algebraic variables");
  }

  const Variable& SymbolicOCP::variable(const std::string& name) const {
    return const_cast<SymbolicOCP*>(this)->variable(name);
  }

  Variable& SymbolicOCP::variable(const std::string& name) {
    // Find the variable
    VarMap::iterator it = varmap_.find(name);
    if (it==varmap_.end()) {
      casadi_error("No such variable: \"" << name << "\".");
    }

    // Return the variable
    return it->second;
  }

  void SymbolicOCP::addVariable(const std::string& name, const Variable& var) {
    // Try to find the component
    if (varmap_.find(name)!=varmap_.end()) {
      stringstream ss;
      ss << "Variable \"" << name << "\" has already been added.";
      throw CasadiException(ss.str());
    }

    // Add to the map of all variables
    varmap_[name] = var;
  }

  std::string SymbolicOCP::qualifiedName(const XmlNode& nn) {
    // Stringstream to assemble name
    stringstream qn;

    for (int i=0; i<nn.size(); ++i) {
      // Add a dot
      if (i!=0) qn << ".";

      // Get the name part
      qn << nn[i].getAttribute("name");

      // Get the index, if any
      if (nn[i].size()>0) {
        int ind;
        nn[i]["exp:ArraySubscripts"]["exp:IndexExpression"]["exp:IntegerLiteral"].getText(ind);
        qn << "[" << ind << "]";
      }
    }

    // Return the name
    return qn.str();
  }

  void SymbolicOCP::generateMuscodDatFile(const std::string& filename,
                                          const Dictionary& mc2_ops) const {
    // Print
    cout << "Generating: " << filename << endl;

    // Create the datfile
    ofstream datfile;
    datfile.open(filename.c_str());
    datfile.precision(numeric_limits<double>::digits10+2);
    datfile << scientific; // This is really only to force a decimal dot,
                           // would be better if it can be avoided

    // Print header
    datfile << "* This function was automatically generated by CasADi" << endl;
    datfile << endl;

    // User set options
    for (Dictionary::const_iterator it=mc2_ops.begin(); it!=mc2_ops.end(); ++it) {
      // Print the name
      datfile << it->first << endl;

      // Get the value
      const GenericType& val = it->second;

      // Print value
      if (val.isInt()) {
        datfile << static_cast<int>(val) << endl;
      } else if (val.isDouble()) {
        datfile << static_cast<double>(val) << endl;
      } else if (val.isString()) {
        datfile << string(val) << endl;
      } else if (val.isIntVector()) {
        vector<int> valv = val;
        for (int k=0; k<valv.size(); ++k) {
          datfile << k << ": " << valv[k] << endl;
        }
      } else if (val.isDoubleVector()) {
        vector<double> valv = val;
        for (int k=0; k<valv.size(); ++k) {
          datfile << k << ": " << valv[k] << endl;
        }
      } else if (val.isStringVector()) {
        vector<string> valv = val;
        for (int k=0; k<valv.size(); ++k) {
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
    if (!h_fix) {
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
    SX p = vertcat(this->pi, this->p);
    if (!p.isEmpty()) {
      datfile << "*  global model parameter start values, scale factors, and bounds" << endl;
      datfile << "p" << endl;
      for (int k=0; k<p.size(); ++k) {
        datfile << k << ": " << start(p[k]) << endl;
      }
      datfile << endl;

      datfile << "p_sca" << endl;
      for (int k=0; k<p.size(); ++k) {
        datfile << k << ": " << nominal(p[k]) << endl;
      }
      datfile << endl;

      datfile << "p_min" << endl;
      for (int k=0; k<p.size(); ++k) {
        datfile << k << ": " << min(p[k]) << endl;
      }
      datfile << endl;

      datfile << "p_max" << endl;
      for (int k=0; k<p.size(); ++k) {
        datfile << k << ": " << max(p[k]) << endl;
      }
      datfile << endl;

      datfile << "p_fix" << endl;
      for (int k=0; k<p.size(); ++k) {
        datfile << k << ": " << (min(p[k])==max(p[k])) << endl;
      }
      datfile << endl;

      datfile << "p_name" << endl;
      for (int k=0; k<p.size(); ++k) {
        datfile << k << ": " << p[k].getName() << endl;
      }
      datfile << endl;

      datfile << "p_unit" << endl;
      for (int k=0; k<p.size(); ++k) {
        datfile << k << ": " << unit(p[k]) << endl;
      }
      datfile << endl;
    }

    // Differential state properties
    if (!this->x.isEmpty()) {
      datfile << "*  differential state start values, scale factors, and bounds" << endl;
      datfile << "sd(*,*)" << endl;
      for (int k=0; k<this->x.size(); ++k) {
        datfile << k << ": " << start(this->x[k]) << endl;
      }
      datfile << endl;

      datfile << "sd_sca(*,*)" << endl;
      for (int k=0; k<this->x.size(); ++k) {
        datfile << k << ": " << nominal(this->x[k]) << endl;
      }
      datfile << endl;

      datfile << "sd_min(*,*)" << endl;
      for (int k=0; k<this->x.size(); ++k) {
        datfile << k << ": " << min(this->x[k]) << endl;
      }
      datfile << endl;

      datfile << "sd_max(*,*)" << endl;
      for (int k=0; k<this->x.size(); ++k) {
        datfile << k << ": " << max(this->x[k]) << endl;
      }
      datfile << endl;

      datfile << "sd_fix(*,*)" << endl;
      for (int k=0; k<this->x.size(); ++k) {
        datfile << k << ": " << (min(this->x[k])==max(this->x[k])) << endl;
      }
      datfile << endl;

      datfile << "xd_name" << endl;
      for (int k=0; k<this->x.size(); ++k) {
        datfile << k << ": " << this->x[k].getName() << endl;
      }
      datfile << endl;

      datfile << "xd_unit" << endl;
      for (int k=0; k<this->x.size(); ++k) {
        datfile << k << ": " << unit(this->x[k]) << endl;
      }
      datfile << endl;
    }

    // Algebraic state properties
    if (!this->z.isEmpty()) {
      datfile << "*  algebraic state start values, scale factors, and bounds" << endl;
      datfile << "sa(*,*)" << endl;
      for (int k=0; k<this->z.size(); ++k) {
        datfile << k << ": " << start(this->z[k]) << endl;
      }
      datfile << endl;

      datfile << "sa_sca(*,*)" << endl;
      for (int k=0; k<this->z.size(); ++k) {
        datfile << k << ": " << nominal(this->z[k]) << endl;
      }
      datfile << endl;

      datfile << "sa_min(*,*)" << endl;
      for (int k=0; k<this->z.size(); ++k) {
        datfile << k << ": " << min(this->z[k]) << endl;
      }
      datfile << endl;

      datfile << "sa_max(*,*)" << endl;
      for (int k=0; k<this->z.size(); ++k) {
        datfile << k << ": " << max(this->z[k]) << endl;
      }
      datfile << endl;

      datfile << "sa_fix(*,*)" << endl;
      for (int k=0; k<this->z.size(); ++k) {
        datfile << k << ": " << (min(this->z[k])==max(this->z[k])) << endl;
      }
      datfile << endl;

      datfile << "xa_name" << endl;
      for (int k=0; k<this->z.size(); ++k) {
        datfile << k << ": " << this->z[k].getName() << endl;
      }
      datfile << endl;

      datfile << "xa_unit" << endl;
      for (int k=0; k<this->z.size(); ++k) {
        datfile << k << ": " << unit(this->z[k]) << endl;
      }
      datfile << endl;
    }

    // Control properties
    if (!this->u.isEmpty()) {
      datfile << "* control start values, scale factors, and bounds" << endl;
      datfile << "u(*,*)" << endl;
      for (int k=0; k<this->u.size(); ++k) {
        datfile << k << ": " << start(this->u[k]) << endl;
      }
      datfile << endl;

      datfile << "u_sca(*,*)" << endl;
      for (int k=0; k<this->u.size(); ++k) {
        datfile << k << ": " << nominal(this->u[k]) << endl;
      }
      datfile << endl;

      datfile << "u_min(*,*)" << endl;
      for (int k=0; k<this->u.size(); ++k) {
        datfile << k << ": " << min(this->u[k]) << endl;
      }
      datfile << endl;

      datfile << "u_max(*,*)" << endl;
      for (int k=0; k<this->u.size(); ++k) {
        datfile << k << ": " << max(this->u[k]) << endl;
      }
      datfile << endl;

      datfile << "u_fix(*,*)" << endl;
      for (int k=0; k<this->u.size(); ++k) {
        datfile << k << ": " << (min(this->u[k])==max(this->u[k])) << endl;
      }
      datfile << endl;

      datfile << "u_name" << endl;
      for (int k=0; k<this->u.size(); ++k) {
        datfile << k << ": " << this->u[k].getName() << endl;
      }
      datfile << endl;

      datfile << "u_unit" << endl;
      for (int k=0; k<this->u.size(); ++k) {
        datfile << k << ": " << unit(this->u[k]) << endl;
      }
      datfile << endl;
    }

    // Close the datfile
    datfile.close();
  }

  SX SymbolicOCP::operator()(const std::string& name) const {
    return variable(name).v;
  }

  SX SymbolicOCP::der(const std::string& name) const {
    return variable(name).d;
  }

  SX SymbolicOCP::der(const SX& var) const {
    casadi_assert(var.isVector() && var.isSymbolic());
    SX ret = SX::zeros(var.sparsity());
    for (int i=0; i<ret.size(); ++i) {
      ret[i] = der(var.at(i).getName());
    }
    return ret;
  }

  SX SymbolicOCP::beq(const std::string& name) const {
    return variable(name).beq;
  }

  SX SymbolicOCP::beq(const SX& var) const {
    return attribute(&SymbolicOCP::beq, var);
  }

  void SymbolicOCP::setBeq(const std::string& name, const SX& val) {
    variable(name).beq = val.toScalar();
  }

  void SymbolicOCP::setBeq(const SX& var, const SX& val) {
    setAttribute(&SymbolicOCP::setBeq, var, val);
  }

  SX SymbolicOCP::ode(const std::string& name) const {
    return variable(name).ode;
  }

  SX SymbolicOCP::ode(const SX& var) const {
    casadi_assert(var.isVector() && var.isSymbolic());
    SX ret = SX::zeros(var.sparsity());
    for (int i=0; i<ret.size(); ++i) {
      ret[i] = ode(var.at(i).getName());
    }
    return ret;
  }

  void SymbolicOCP::setOde(const std::string& name, const SX& val) {
    variable(name).ode = val.toScalar();
  }

  void SymbolicOCP::setOde(const SX& var, const SX& val) {
    casadi_assert(var.isVector() && var.isSymbolic());
    casadi_assert(var.sparsity()==val.sparsity());
    for (int i=0; i<var.size(); ++i) {
      setOde(var.at(i).getName(), val.at(i));
    }
  }

  SX SymbolicOCP::atTime(const std::string& name, double t, bool allocate) const {
    return variable(name).atTime(t, allocate);
  }

  SX SymbolicOCP::atTime(const std::string& name, double t, bool allocate) {
    return variable(name).atTime(t, allocate);
  }

  void SymbolicOCP::separateAlgebraic() {
    // Quick return if no s
    if (this->s.isEmpty()) return;

    // We investigate the interdependencies in sdot -> dae
    vector<SX> f_in;
    f_in.push_back(der(this->s));
    SXFunction f(f_in, this->dae);
    f.init();

    // Number of s
    int ns = f.input().size();
    casadi_assert(f.output().size()==ns);

    // Input/output arrays
    bvec_t* f_sdot = reinterpret_cast<bvec_t*>(f.input().ptr());
    bvec_t* f_dae = reinterpret_cast<bvec_t*>(f.output().ptr());

    // First find out which equations depend on sdot
    f.spInit(true);

    // Seed all inputs
    std::fill(f_sdot, f_sdot+ns, bvec_t(1));

    // Propagate to f_dae
    std::fill(f_dae, f_dae+ns, bvec_t(0));
    f.spEvaluate(true);

    // Get the new differential and algebraic equations
    SX new_dae, new_alg;
    for (int i=0; i<ns; ++i) {
      if (f_dae[i]==bvec_t(1)) {
        new_dae.append(this->dae[i]);
      } else {
        casadi_assert(f_dae[i]==bvec_t(0));
        new_alg.append(this->dae[i]);
      }
    }

    // Now find out what sdot enter in the equations
    f.spInit(false);

    // Seed all outputs
    std::fill(f_dae, f_dae+ns, bvec_t(1));

    // Propagate to f_sdot
    std::fill(f_sdot, f_sdot+ns, bvec_t(0));
    f.spEvaluate(false);

    // Get the new algebraic variables and new states
    SX new_s, new_z;
    for (int i=0; i<ns; ++i) {
      if (f_sdot[i]==bvec_t(1)) {
        new_s.append(this->s[i]);
      } else {
        casadi_assert(f_sdot[i]==bvec_t(0));
        new_z.append(this->s[i]);
      }
    }

    // Make sure split was successful
    casadi_assert(new_dae.size()==new_s.size());

    // Divide up the s and dae
    this->dae = new_dae;
    this->s = new_s;
    this->alg.append(new_alg);
    this->z.append(new_z);
  }

  std::string SymbolicOCP::unit(const std::string& name) const {
    return variable(name).unit;
  }

  std::string SymbolicOCP::unit(const SX& var) const {
    casadi_assert_message(!var.isVector() && var.isSymbolic(),
                          "SymbolicOCP::unit: Argument must be a symbolic vector");
    if (var.isEmpty()) {
      return "n/a";
    } else {
      string ret = unit(var.at(0).getName());
      for (int i=1; i<var.size(); ++i) {
        casadi_assert_message(ret == unit(var.at(i).getName()),
                              "SymbolicOCP::unit: Argument has mixed units");
      }
      return ret;
    }
  }

  void SymbolicOCP::setUnit(const std::string& name, const std::string& val) {
    variable(name).unit = val;
  }

  double SymbolicOCP::nominal(const std::string& name) const {
    return variable(name).nominal;
  }

  void SymbolicOCP::setNominal(const std::string& name, double val) {
    variable(name).nominal = val;
  }

  std::vector<double> SymbolicOCP::nominal(const SX& var) const {
    casadi_assert_message(var.isVector() && var.isSymbolic(),
                          "SymbolicOCP::nominal: Argument must be a symbolic vector");
    std::vector<double> ret(var.size());
    for (int i=0; i<ret.size(); ++i) {
      ret[i] = nominal(var.at(i).getName());
    }
    return ret;
  }

  void SymbolicOCP::setNominal(const SX& var, const std::vector<double>& val) {
    casadi_assert_message(var.isVector() && var.isSymbolic(),
                          "SymbolicOCP::nominal: Argument must be a symbolic vector");
    casadi_assert_message(var.size()==var.size(), "SymbolicOCP::nominal: Dimension mismatch");
    for (int i=0; i<val.size(); ++i) {
      setNominal(var.at(i).getName(), val.at(i));
    }
  }

  std::vector<double> SymbolicOCP::attribute(getAtt f, const SX& var, bool normalized) const {
    casadi_assert_message(var.isVector() && var.isSymbolic(),
                          "SymbolicOCP::attribute: Argument must be a symbolic vector");
    std::vector<double> ret(var.size());
    for (int i=0; i<ret.size(); ++i) {
      ret[i] = (this->*f)(var.at(i).getName(), normalized);
    }
    return ret;
  }

  SX SymbolicOCP::attribute(getAttS f, const SX& var) const {
    casadi_assert_message(var.isVector() && var.isSymbolic(),
                          "SymbolicOCP::attribute: Argument must be a symbolic vector");
    SX ret = SX::zeros(var.sparsity());
    for (int i=0; i<ret.size(); ++i) {
      ret[i] = (this->*f)(var.at(i).getName());
    }
    return ret;
  }

  void SymbolicOCP::setAttribute(setAtt f, const SX& var, const std::vector<double>& val,
                                 bool normalized) {
    casadi_assert_message(var.isVector() && var.isSymbolic(),
                          "SymbolicOCP::setAttribute: Argument must be a symbolic vector");
    casadi_assert_message(var.size()==val.size(), "SymbolicOCP::setAttribute: Dimension mismatch");
    for (int i=0; i<val.size(); ++i) {
      (this->*f)(var.at(i).getName(), val.at(i), normalized);
    }
  }

  void SymbolicOCP::setAttribute(setAttS f, const SX& var, const SX& val) {
    casadi_assert_message(var.isVector() && var.isSymbolic(),
                          "SymbolicOCP::setAttribute: Argument must be a symbolic vector");
    casadi_assert_message(var.sparsity()==val.sparsity(),
                          "SymbolicOCP::setAttribute: Sparsity mismatch");
    for (int i=0; i<val.size(); ++i) {
      (this->*f)(var.at(i).getName(), val.at(i));
    }
  }

  SX SymbolicOCP::min(const std::string& name) const {
    return variable(name).min;
  }

  SX SymbolicOCP::min(const SX& var) const {
    return attribute(&SymbolicOCP::min, var);
  }

  void SymbolicOCP::setMin(const std::string& name, const SX& val) {
    variable(name).min = val.toScalar();
  }

  void SymbolicOCP::setMin(const SX& var, const SX& val) {
    setAttribute(&SymbolicOCP::setMin, var, val);
  }

  SX SymbolicOCP::max(const std::string& name) const {
    return variable(name).max;
  }

  SX SymbolicOCP::max(const SX& var) const {
    return attribute(&SymbolicOCP::max, var);
  }

  void SymbolicOCP::setMax(const std::string& name, const SX& val) {
    variable(name).max = val.toScalar();
  }

  void SymbolicOCP::setMax(const SX& var, const SX& val) {
    setAttribute(&SymbolicOCP::setMax, var, val);
  }

  SX SymbolicOCP::initialGuess(const std::string& name) const {
    return variable(name).initialGuess;
  }

  SX SymbolicOCP::initialGuess(const SX& var) const {
    return attribute(&SymbolicOCP::initialGuess, var);
  }

  void SymbolicOCP::setInitialGuess(const std::string& name, const SX& val) {
    variable(name).initialGuess = val.toScalar();
  }

  void SymbolicOCP::setInitialGuess(const SX& var, const SX& val) {
    setAttribute(&SymbolicOCP::setInitialGuess, var, val);
  }

  double SymbolicOCP::start(const std::string& name, bool normalized) const {
    const Variable& v = variable(name);
    return normalized ? v.start / v.nominal : v.start;
  }

  std::vector<double> SymbolicOCP::start(const SX& var, bool normalized) const {
    return attribute(&SymbolicOCP::start, var, normalized);
  }

  void SymbolicOCP::setStart(const std::string& name, double val, bool normalized) {
    Variable& v = variable(name);
    v.start = normalized ? val*v.nominal : val;
  }

  void SymbolicOCP::setStart(const SX& var, const std::vector<double>& val, bool normalized) {
    setAttribute(&SymbolicOCP::setStart, var, val, normalized);
  }

  double SymbolicOCP::derivativeStart(const std::string& name, bool normalized) const {
    const Variable& v = variable(name);
    return normalized ? v.derivativeStart / v.nominal : v.derivativeStart;
  }

  std::vector<double> SymbolicOCP::derivativeStart(const SX& var, bool normalized) const {
    return attribute(&SymbolicOCP::derivativeStart, var, normalized);
  }

  void SymbolicOCP::setDerivativeStart(const std::string& name, double val, bool normalized) {
    Variable& v = variable(name);
    v.derivativeStart = normalized ? val*v.nominal : val;
  }

  void SymbolicOCP::setDerivativeStart(const SX& var, const std::vector<double>& val,
                                       bool normalized) {
    setAttribute(&SymbolicOCP::setDerivativeStart, var, val, normalized);
  }

} // namespace casadi

