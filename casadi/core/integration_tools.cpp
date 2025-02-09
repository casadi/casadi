/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#include "integration_tools.hpp"
#include "integrator.hpp"
#include "rootfinder.hpp"
#include "polynomial.hpp"
#include "nlpsol.hpp"

namespace casadi {

const long double legendre_points1[] = { 0.50000000000000000000 };
const long double legendre_points2[] =
    { 0.21132486540518713447, 0.78867513459481286553 };
const long double legendre_points3[] =
  { 0.11270166537925824235, 0.50000000000000000000,
    0.88729833462074170214 };
const long double legendre_points4[] =
  { 0.06943184420297354720, 0.33000947820757187134,
    0.66999052179242823968, 0.93056815579702623076 };
const long double legendre_points5[] =
  { 0.04691007703066807366, 0.23076534494715861268,
    0.49999999999999994449, 0.76923465505284149835, 0.95308992296933192634 };
const long double legendre_points6[] =
  { 0.03376524289842347537, 0.16939530676686742616,
    0.38069040695840172805, 0.61930959304159849399, 0.83060469323313235179,
    0.96623475710157580298 };
const long double legendre_points7[] =
  { 0.02544604382862047931, 0.12923440720030288098,
    0.29707742431130129690, 0.50000000000000000000, 0.70292257568869853657,
    0.87076559279969734106, 0.97455395617137896558 };
const long double legendre_points8[] =
  { 0.01985507175123157886, 0.10166676129318691357,
    0.23723379504183561561, 0.40828267875217505445, 0.59171732124782483453,
    0.76276620495816449541, 0.89833323870681347501, 0.98014492824876797705 };
const long double legendre_points9[] =
  { 0.01591988024618706810, 0.08198444633668211523,
    0.19331428364970504319, 0.33787328829809543107, 0.49999999999999988898,
    0.66212671170190451342, 0.80668571635029517886, 0.91801555366331766272,
    0.98408011975381259884 };
const long double* legendre_points[] =
  { nullptr, legendre_points1, legendre_points2, legendre_points3, legendre_points4,
    legendre_points5, legendre_points6, legendre_points7, legendre_points8, legendre_points9};

// Radau collocation points
const long double radau_points1[] =
  { 1.00000000000000000000 };
const long double radau_points2[] =
  { 0.33333333333333337034, 1.00000000000000000000 };
const long double radau_points3[] =
  { 0.15505102572168222297, 0.64494897427831787695,
    1.00000000000000000000 };
const long double radau_points4[] =
  { 0.08858795951270420632, 0.40946686444073465694,
    0.78765946176084700170, 1.00000000000000000000 };
const long double radau_points5[] =
  { 0.05710419611451822419, 0.27684301363812369168,
    0.58359043236891683382, 0.86024013565621926247, 1.00000000000000000000 };
const long double radau_points6[] =
  { 0.03980985705146905529, 0.19801341787360787761,
    0.43797481024738621480, 0.69546427335363603106, 0.90146491420117347282,
    1.00000000000000000000 };
const long double radau_points7[] =
  { 0.02931642715978521885, 0.14807859966848435640,
    0.33698469028115418666, 0.55867151877155019069, 0.76923386203005450490,
    0.92694567131974103802, 1.00000000000000000000 };
const long double radau_points8[] =
  { 0.02247938643871305597, 0.11467905316090415413,
    0.26578982278458951338, 0.45284637366944457959, 0.64737528288683043876,
    0.81975930826310761113, 0.94373743946307731001, 1.00000000000000000000 };
const long double radau_points9[] =
  { 0.01777991514736393386, 0.09132360789979432347,
    0.21430847939563035798, 0.37193216458327238438, 0.54518668480342658000,
    0.71317524285556954666, 0.85563374295785443735, 0.95536604471003006012,
    1.00000000000000000000 };
const long double* radau_points[] =
  { nullptr, radau_points1, radau_points2, radau_points3, radau_points4, radau_points5,
    radau_points6, radau_points7, radau_points8, radau_points9};

template<typename RealT>
std::vector<RealT> collocation_pointsGen(casadi_int order, const std::string& scheme) {
  if (scheme=="radau") {
    casadi_assert(order>0 && order<10,
      "Error in collocationPoints(order, scheme): "
      "only order up to 9 supported for scheme 'radau', but got " + str(order) + ".");
    return std::vector<RealT>(radau_points[order], radau_points[order]+order);
  } else if (scheme=="legendre") {
    casadi_assert(order>0 && order<10,
      "Error in collocationPoints(order, scheme): "
      "only order up to 9 supported for scheme 'legendre', but got " + str(order) + ".");
    return std::vector<RealT>(legendre_points[order], legendre_points[order]+order);
  } else {
    casadi_error("Error in collocationPoints(order, scheme): unknown scheme '"
                  + scheme + "'. Select one of 'radau', 'legendre'.");
  }
}

std::vector<double> collocation_points(casadi_int order, const std::string& scheme) {
  return collocation_pointsGen<double>(order, scheme);
}

std::vector<long double> collocation_pointsL(casadi_int order, const std::string& scheme) {
  return collocation_pointsGen<long double>(order, scheme);
}

Function simpleRK(Function f, casadi_int N, casadi_int order) {
  // Consistency check
  casadi_assert(N>=1,
    "Parameter N (number of steps) must be at least 1, but got " + str(N) + ".");
  casadi_assert(order==4, "Only RK order 4 is supported now.");
  casadi_assert(f.n_in()==2, "Function must have two inputs: x and p");
  casadi_assert(f.n_out()==1, "Function must have one outputs: dot(x)");

  MX x0 = MX::sym("x0", f.sparsity_in(0));
  MX p = MX::sym("p", f.sparsity_in(1));
  MX h = MX::sym("h");

  std::vector<double> b(order);
  b[0]=1.0/6;
  b[1]=1.0/3;
  b[2]=1.0/3;
  b[3]=1.0/6;

  std::vector<double> c(order);
  c[0]=0;
  c[1]=1.0/2;
  c[2]=1.0/2;
  c[3]=1;

  std::vector< std::vector<double> > A(order-1);
  A[0].resize(1);
  A[0][0]=1.0/2;
  A[1].resize(2);
  A[1][0]=0;A[1][1]=1.0/2;
  A[2].resize(3);
  A[2][0]=0;
  A[2][1]=0;A[2][2]=1;

  // Time step
  MX dt = h/N;

  std::vector<MX> k(order);
  std::vector<MX> f_arg(2);

  // Integrate
  MX xf = x0;
  for (casadi_int i=0; i<N; ++i) {
    for (casadi_int j=0; j<order; ++j) {
      MX xL = 0;
      for (casadi_int jj=0; jj<j; ++jj) {
        xL += k.at(jj)*A.at(j-1).at(jj);
      }
      f_arg[0] = xf+xL;
      f_arg[1] = p;
      k[j] = dt*f(f_arg).at(0);
    }

    for (casadi_int j=0; j<order; ++j) {
      xf += b.at(j)*k.at(j);
    }
  }

  // Form discrete-time dynamics
  return Function("F", {x0, p, h}, {xf}, {"x0", "p", "h"}, {"xf"});
}

void collocation_interpolators(const std::vector<double> & tau,
                              std::vector< std::vector<double> > &C, std::vector< double > &D) {
  // Find the degree of the interpolation
  casadi_int deg = tau.size();

  // Include zero
  std::vector<double> etau_root = tau;
  etau_root.insert(etau_root.begin(), 0);

  // Allocate storage space for resulting coefficients
  C.resize(deg+1);
  for (casadi_int i=0;i<deg+1;++i) {
    C[i].resize(deg+1);
  }
  D.resize(deg+1);

  // Collocation point
  SX tau_sym = SX::sym("tau");

  // For all collocation points
  for (casadi_int j=0; j<deg+1; ++j) {
    // Construct Lagrange polynomials to get the polynomial basis at the collocation point
    SX L = 1;
    for (casadi_int j2=0; j2<deg+1; ++j2) {
      if (j2 != j) {
        L *= (tau_sym-etau_root[j2])/(etau_root[j]-etau_root[j2]);
      }
    }

    Function lfcn("lfcn", {tau_sym}, {L});

    // Evaluate the polynomial at the final time to get the
    // coefficients of the continuity equation
    D[j] = lfcn(std::vector<DM>{1.}).at(0)->front();

    // Evaluate the time derivative of the polynomial at all collocation points to
    // get the coefficients of the continuity equation
    Function tfcn("tfcn", {tau_sym}, {tangent(L, tau_sym)});
    for (casadi_int j2=0; j2<deg+1; ++j2) {
      C[j2][j] =  tfcn(std::vector<DM>{etau_root[j2]}).at(0)->front();
    }
  }
}

void collocation_coeff(const std::vector<double> & tau,
                              DM &C, DM &D, DM &B) {
  // Find the degree of the interpolation
  casadi_int deg = tau.size();

  // Include zero
  std::vector<double> etau_root = tau;
  etau_root.insert(etau_root.begin(), 0);

  // Check if tau ends with '1' (cfr. radau scheme)
  bool has_end = tau.back()==1;

  // Coefficients of the collocation equation
  std::vector<std::vector<double> > C_(deg+1, std::vector<double>(deg+1, 0));

  // Coefficients of the continuity equation
  std::vector<double> D_(deg+1, 0);

  // Coefficients of the quadratures
  std::vector<double> B_(deg+1, 0);

  // For all collocation points
  for (casadi_int j=0; j<deg+1; ++j) {

    // Construct Lagrange polynomials to get the polynomial basis at the collocation point
    Polynomial p = 1;
    for (casadi_int r=0; r<deg+1; ++r) {
      if (r!=j) {
        p *= Polynomial(-etau_root[r], 1)/(etau_root[j]-etau_root[r]);
      }
    }

    // Evaluate the polynomial at the final time to get the
    // coefficients of the continuity equation
    if (has_end) {
      D_[j] = j==deg ? 1 : 0;
    } else {
      D_[j] = p(1.0);
    }
    // Evaluate the time derivative of the polynomial at all collocation points to
    // get the coefficients of the continuity equation
    Polynomial dp = p.derivative();
    for (casadi_int r=0; r<deg+1; ++r) {
      C_[j][r] = dp(etau_root[r]);
    }

    // Integrate polynomial to get the coefficients of the quadratures
    Polynomial ip = p.anti_derivative();
    B_[j] = ip(1.0);
  }
  C = DM(C_);
  C = C(Slice(), Slice(1, deg+1)); // NOLINT(cppcoreguidelines-slicing)
  D = DM(D_);
  B = DM(std::vector<double>(B_.begin()+1, B_.end()));
}

Function simpleIRK(Function f, casadi_int N, casadi_int order, const std::string& scheme,
                      const std::string& solver,
                      const Dict& solver_options) {
  // Consistency check
  casadi_assert(N>=1,
    "Parameter N (number of steps) must be at least 1, but got " + str(N) + ".");
  casadi_assert(f.n_in()==2, "Function must have two inputs: x and p");
  casadi_assert(f.n_out()==1, "Function must have one outputs: dot(x)");

  // Obtain collocation points
  std::vector<double> tau_root = collocation_points(order, scheme);

  // Retrieve collocation interpolating matrices
  std::vector < std::vector <double> > C;
  std::vector < double > D;
  collocation_interpolators(tau_root, C, D);

  // Inputs of constructed function
  MX x0 = MX::sym("x0", f.sparsity_in(0));
  MX p = MX::sym("p", f.sparsity_in(1));
  MX h = MX::sym("h");

  // Time step
  MX dt = h/N;

  // Implicitly defined variables
  MX v = MX::sym("v", repmat(x0.sparsity(), order));
  std::vector<MX> x = vertsplit(v, x0.size1());
  x.insert(x.begin(), x0);

  // Collect the equations that implicitly define v
  std::vector<MX> V_eq, f_in(2), f_out;
  for (casadi_int j=1; j<order+1; ++j) {
    // Expression for the state derivative at the collocation point
    MX xp_j = 0;
    for (casadi_int r=0; r<=order; ++r) xp_j+= C[j][r]*x[r];

    // Collocation equations
    f_in[0] = x[j];
    f_in[1] = p;
    f_out = f(f_in);
    V_eq.push_back(dt*f_out.at(0)-xp_j);
  }

  // Root-finding function
  Function rfp("rfp", {v, x0, p, h}, {vertcat(V_eq)});

  // Create a implicit function instance to solve the system of equations
  Function ifcn = rootfinder("ifcn", solver, rfp, solver_options);

  // Get state at end time
  MX xf = x0;
  for (casadi_int k=0; k<N; ++k) {
    std::vector<MX> ifcn_out = ifcn({repmat(xf, order), xf, p, h});
    x = vertsplit(ifcn_out[0], xf.size1());

    // State at end of step
    xf = D[0]*xf;
    for (casadi_int i=1; i<=order; ++i) {
      xf += D[i]*x[i-1];
    }
  }

  // Form discrete-time dynamics
  return Function("F", {x0, p, h}, {xf}, {"x0", "p", "h"}, {"xf"});
}

Function simpleIntegrator(Function f, const std::string& plugin,
                          const Dict& plugin_options) {
  // Consistency check
  casadi_assert(f.n_in()==2, "Function must have two inputs: x and p");
  casadi_assert(f.n_out()==1, "Function must have one outputs: dot(x)");

  // Sparsities
  Sparsity x_sp = f.sparsity_in(0);
  Sparsity p_sp = f.sparsity_in(1);

  // Wrapper function inputs
  MX x = MX::sym("x", x_sp);
  MX u = MX::sym("u", vertcat(Sparsity::scalar(), vec(p_sp))); // augment p with t

  // Normalized xdot
  casadi_int u_offset[] = {0, 1, 1+p_sp.size1()};
  std::vector<MX> pp = vertsplit(u, std::vector<casadi_int>(u_offset, u_offset+3));
  MX h = pp[0];
  MX p = reshape(pp[1], p_sp.size());
  MX f_in[] = {x, p};
  MX xdot = f(std::vector<MX>(f_in, f_in+2)).at(0);
  xdot *= h;

  // Form DAE function
  MXDict dae = {{"x", x}, {"p", u}, {"ode", xdot}};

  // Create integrator function with normalized time from 0 to 1
  Function ifcn = integrator("integrator", plugin, dae, plugin_options);

  // Inputs of constructed function
  MX x0 = MX::sym("x0", x_sp);
  p = MX::sym("p", p_sp);
  h = MX::sym("h");

  // State at end
  MX xf = ifcn(MXDict{{"x0", x0}, {"p", vertcat(h, vec(p))}}).at("xf");

  // Form discrete-time dynamics
  return Function("F", {x0, p, h}, {xf}, {"x0", "p", "h"}, {"xf"});
}


std::vector<casadi_int> invert_lookup(const std::vector<casadi_int>& lookup) {
  std::vector<casadi_int> ret(lookup.size(), -1);
  for (casadi_int i=0;i<lookup.size();++i) {
    casadi_int e = lookup[i];
    if (e>=0) {
      ret[e] = i;
    }
  }
  return ret;
}

namespace IndexReduction {

  struct EquationStruct;
  /*
  struct VariableStruct {
    std::vector<struct EquationStruct*> eqs;
    // Eligability to serve as candidate, one for each pass during phase 1
    std::vector<bool> eligible1;
    struct EquationStruct* assign;
    // Which variable is this variable's derivative?
    struct VariableStruct* der;
    // Which variable produces this variable by differentiation?
    struct VariableStruct* ;
    bool visited;
    // Eligability to serve as candidate during phase 2
    bool eligible2;
  };*/

  struct VariableStruct {
    std::vector<struct EquationStruct*> eqs;
    std::vector<struct EquationStruct*> eqs0;
    struct EquationStruct* assign = nullptr;
    // Which variable is this variable's derivative?
    struct VariableStruct* der = nullptr;
    // Which variable produces this variable by differentiation?
    struct VariableStruct* der_inv = nullptr;
    casadi_int i; // Position in Graph::V
    bool visited = false;
    bool deleted = false;
  };

  struct EquationStruct {
    std::vector<struct VariableStruct*> vars;
    std::vector<struct VariableStruct*> vars0;
    struct VariableStruct* assign = nullptr;
    // Which equation is this equations's derivative?
    struct EquationStruct* der = nullptr;
    // Which equations produces this equations by differentiation?
    struct EquationStruct* der_inv = nullptr;
    casadi_int i; // Position in Graph::E
    bool visited = false;
  };

  typedef struct VariableStruct Variable;
  typedef struct EquationStruct Equation;

  // Bipartite graph
  struct GraphStruct {
    std::vector<Variable*> V;
    std::vector<Equation*> E;
    casadi_int ncol_orig;
    casadi_int nrow_orig;
  };

  typedef struct GraphStruct Graph;

  void graph_add_der(Graph& G, Variable* v) {
    // Push new equation to the graph
    G.V.push_back(new Variable());
    Variable* v_new = G.V.back();

    // Get a position
    v_new->i = G.V.size()-1;

    // Add derivative relationship between new and old variable
    v_new->der_inv = v;
    v->der = v_new;
  }

  void add_variable(Equation* e, Variable* v) {
    auto it = std::find(e->vars0.begin(), e->vars0.end(), v);
    if (it==e->vars0.end()) {
      e->vars0.push_back(v);
      if (!v->deleted) {
        e->vars.push_back(v);
        v->eqs.push_back(e);
      }
    }
  }

  void graph_add_der(Graph& G, Equation* e, bool add_old=false) {
    // Push new equation to the graph
    G.E.push_back(new Equation());

    Equation* e_new = G.E.back();

    // Get a position
    e_new->i = G.E.size()-1;

    // Add derivative relationship between new and old equation
    e_new->der_inv = e;
    e->der = e_new;

    // Loop over all variables in old equation
    for (auto* v : e->vars0) {
      // We depend on the old variable
      if (add_old) {
        add_variable(e_new, v);
      }
      // And we depend on its derivative (create if not present)
      if (!v->der) graph_add_der(G, v);
      add_variable(e_new, v->der);
    }

  }

  bool dfs_match_pantelides(Equation* i) {

    // Pantelides alg 3.2: Colour i (1)
    i->visited = true;

    // Look for unassigned candidate variables
    // Pantelides alg 3.2: (2) If a V-node j exists such that edge (i-j) exists...
    for (auto* j : i->vars) {
      // Pantelides alg 3.2: ... and ASSIGN (j)= 0 then:
      if (j->assign==nullptr) {
        // Found such variable
        // Assign
        // Pantelides alg 3.2: (2b) Set ASSIGN(j)=i
        j->assign = i; i->assign = j;

        // Pantelides alg 3.2: (2a) set PATHFOUND = TRUE, (2c) return
        return true; // End successfully
      }
    }

    // Look for assigned candidate variables
    // Pantelides alg 3.2: (3) For every j such that edge (i-j) exists
    for (auto j : i->vars) {
      // Pantelides alg 3.2: ... and j is uncoloured do:
      if (!j->visited) {
        // Pantelides alg 3.2: (3a) Colour j
        j->visited = true;

        // Pantelides alg 3.2: (3b) Set k = ASSIGN (j)
        Equation* k = j->assign;

        // Pantelides alg 3.2: (3c) AUGMENTPATH (k, PATHFOUND)
        if (dfs_match_pantelides(k)) {
          // Re-assignment
          // Pantelides alg 3.2: (3d-1) Set ASSIGN (j) = i
          j->assign = i; i->assign = j;
          // Pantelides alg 3.2: (3d-2) Return
          return true;
        }

      }
    }

    // Exhausted all options; return unsuccessfully
    return false;
  }


  /** \brief Perform Pantelides algorithm for DAE index reduction
   * 
   * The algorithm works purely on structure: no symbolics equations are used.
   *
   * \param graph   Structural relation between equations (columns) and variables (rows)
   * \param var_map Indicate for each variable where its derivative can be found in the variable list
   *                if var_map[i]>=0: var[var_map[i]] == dot(var[i])
   *                Size will increase to accommodate new variables introduced by index reduction 
   * \param eq_map  Indicate for each equation if it should be differentiated
   *                if eq_map[i]>=0:  eq[eq_map[i]] == dot(eq[i])
   *                Size will increase over the original number of equations to accommodate extra equations introduced by index reduction
  */
  bool dae_struct_detect_pantelides(Graph& G,
    std::vector<casadi_int>& var_map, std::vector<casadi_int>& eq_map,
    casadi_int max_iter) {

    // Pantelides alg 4.1: Set N'=N do (2)
    casadi_int Np = G.E.size();
    // Pantelides alg 4.1: For k=1 to N' do (3)
    for (casadi_int k=0;k<Np;++k) {
      // Pantelides alg 4.1: Set i=k (3a)
      Equation* i = G.E[k];

      casadi_int iter = 0;
      // Pantelides alg 4.1: Repeat (3b)
      while (true) {
        if (iter++ > max_iter) return false;
        // Pantelides alg 4.1: delete all V-nodes with A!=0 and their incident edges
        //       * from the graph (3b-1)
        for (auto* v : G.V) {
          if (v->der) {
            v->deleted = true;
            v->eqs.clear();
            for (auto* e : G.E) {
              auto it = std::find(e->vars.begin(), e->vars.end(), v);
              if (it!=e->vars.end()) e->vars.erase(it);
            }
          }
        }

        // Pantelides alg 4.1: Designate all nodes as "uncoloured" (3b-2)
        for (auto* v : G.V) v->visited = false;
        for (auto* e : G.E) e->visited = false;

        // Pantelides alg 4.1: Set PATHFOUND= FALSE        (3b-3)
        // Pantelides alg 4.1: AUGMENTPATH (i, PATHFOUND)  (3b-4)
        // Pantelides alg 4.1: Until PATHFOUND  (3c)
        if (dfs_match_pantelides(i)) break; // out of while loop

        // Pantelides alg 4.1: If PATHFOUND FALSE then (3b-5)
        casadi_int n = G.V.size();
        // Pantelides alg 4.1: For every coloured V-node j do (3b-5-i)
        for (casadi_int jj=0;jj<n;++jj) {
          Variable* j = G.V[jj];
          if (j->visited && !j->deleted) {
            graph_add_der(G, j);
          }
        }

        // Pantelides alg 4.1: For every coloured E-node l do (3b-5-ii)
        n = G.E.size();
        for (casadi_int ll=0;ll<n;++ll) {
          Equation* l = G.E[ll];
          if (l->visited) graph_add_der(G, l, true);
        }

        for (casadi_int ll=n;ll<G.E.size();++ll) {
          Equation* l = G.E[ll];
          bool valid = false;
          for (auto* v : l->vars0) {
            if (!v->eqs.empty()) valid = true;
          }
          if (!valid) return false;
        }

        // Pantelides alg 4.1: For every coloured V-node j set
        // ASSIGN (A(j))= B(ASSIGN (j))  (3b-5-iii)
        for (auto* j : G.V) {
          if (j->visited && !j->deleted) {
            j->der->assign = j->assign->der;
            j->assign->der->assign = j->der;
          }
        }

        // Pantelides alg 4.1: Set i=B(i) (3b-5-iv)
        i = i->der;

      }
    }
    return true;
  }

  void dae_struct_detect(const std::string& algorithm,
    const Sparsity& graph, std::vector<casadi_int>& var_map,
    std::vector<casadi_int>& eq_map,
    casadi_int max_iter) {
    // Input sanitization
    casadi_assert(var_map.size()==graph.size2(),
      "var_map size must match graph columns.");

    // Graph structure
    Graph G;

    // Allocate space for node representation
    G.E.resize(graph.size1(), nullptr);
    G.V.resize(graph.size2(), nullptr);
    G.nrow_orig = graph.size1();
    G.ncol_orig = graph.size2();

    for (auto*& e : G.E) e = new Equation();
    for (auto*& v : G.V) v = new Variable();

    // Set positions
    casadi_int i;
    i=0; for (auto* e : G.E) e->i = i++;
    i=0; for (auto* v : G.V) v->i = i++;

    // Create edges using incidence sparsity
    const casadi_int* colind = graph.colind();
    const casadi_int* row = graph.row();

    // Loop over incidence columns
    for (casadi_int c=0;c<graph.size2();++c) {
      // Loop over nonzeros
      for (casadi_int j=colind[c]; j<colind[c+1]; ++j) {
        casadi_int r = row[j];
        G.V[c]->eqs.push_back(G.E[r]); // Non-monotone
        G.E[r]->vars.push_back(G.V[c]);
      }
    }

    // Process var_map
    for (casadi_int i=0;i<var_map.size();++i) {
      Variable* v = G.V[i];
      casadi_int i_der = var_map[i];
      if (i_der>=0) {
        Variable* e = G.V[i_der];
        v->der = e;
        e->der_inv = v;
      }
    }

    for (auto* v : G.V) v->eqs0 = v->eqs;
    for (auto* e : G.E) e->vars0 = e->vars;

    bool detect = false;
    if (algorithm=="pantelides") {
      detect = dae_struct_detect_pantelides(G, var_map, eq_map, max_iter);
    } else {
      // Clean up
      for (auto* v : G.V) delete v;
      for (auto* e : G.E) delete e;
      casadi_error("Algorithm '" + algorithm + "' not recognized.");
    }
    if (!detect) {
      for (auto* v : G.V) delete v;
      for (auto* e : G.E) delete e;
      casadi_error("Structural detection failed.");
    }

    // Prepare outputs
    var_map.resize(G.V.size());
    eq_map.resize(G.E.size());

    int k = 0;
    // Populate var_map
    for (auto* v : G.V) {
      if (v->der) {
        var_map[v->i] = v->der->i;
      } else  {
        var_map[v->i] = -1;
      }
      casadi_assert_dev(v->i==k);
      k++;
    }

    k = 0;
    // Populate eq_map
    for (auto* e : G.E) {
      if (e->der) {
        eq_map[e->i] = e->der->i;
      } else  {
        eq_map[e->i] = -1;
      }
      casadi_assert_dev(e->i==k);
      k++;
    }

    // Clean up
    for (auto* v : G.V) delete v;
    for (auto* e : G.E) delete e;

  }
}  // namespace IndexReduction


using namespace IndexReduction;

template <class X>
std::map<std::string, X> add_defaults(const std::map<std::string, X>& in,
    const std::vector<std::string>& keys) {
  std::map<std::string, X> ret;
  for (const std::string& k : keys) {
    auto it = in.find(k);
    if (it==in.end()) {
      ret[k] = X(0, 1);
    } else {
      ret[k] = it->second;
    }
  }
  return ret;
}


std::vector<casadi_int> get_orders(const std::vector<casadi_int>& map) {
  std::vector<casadi_int> ret(map.size(), 0);
  for (casadi_int i=0;i<ret.size();++i) {
    if (map[i]>=0) ret[map[i]] = ret[i]+1;
  }
  return ret;
}

std::vector<casadi_int> get_inverse(const std::vector<casadi_int>& map) {
  std::vector<casadi_int> ret(map.size(), -1);
  for (casadi_int i=0;i<ret.size();++i) {
    if (map[i]>=0) ret[map[i]] = i;
  }
  return ret;
}

std::vector<casadi_int> path(const std::vector<casadi_int>& map, casadi_int i_start) {
  std::vector<casadi_int> ret;
  casadi_int i = i_start;
  while (true) {
    casadi_int i_next = map[i];
    if (i_next==-1) break;
    i = i_next;
    ret.push_back(i);
  }
  return ret;
}

template <class X>
const std::map<std::string, X>
reduce_index_gen(const std::map<std::string, X>& dae, Dict& stats, const Dict& opts) {
  double baumgarte_pole_ = 0;
  std::string algorithm = "pantelides";
  casadi_int max_iter = 500;
  // Option parsing
  for (auto&& op : opts) {
    if (op.first=="baumgarte_pole") {
      baumgarte_pole_ = op.second;
    } else if (op.first=="algorithm") {
      algorithm = op.second.to_string();
    } else if (op.first=="max_iter") {
      max_iter = op.second;
    } else {
      casadi_error("Unknown option '" + op.first + "'.");
    }
  }

  // Dae input sanitization

  for (const auto& e : dae) {
    casadi_assert(e.second.is_column() && e.second.is_dense(),
      "Dense column vector expected for key '" + e.first + "'. "
      "Got " + e.second.dim(true) + " instead.");
  }

  X x   = get_from_dict(dae, "x",   X(0, 1));
  X ode = get_from_dict(dae, "ode", X(0, 1));
  casadi_assert(x.numel()==ode.numel(),
    "Size of explicit differential state vector (x: " + str(x.numel()) + ") must "
    "match size of ode right-hand-side (ode: " + str(ode.numel()) + ").");
  casadi_assert(x.size1()==0, "Explicit differential states not supported yet.");

  X x_impl  = get_from_dict(dae, "x_impl",  X(0, 1));
  X dx_impl = get_from_dict(dae, "dx_impl", X(0, 1));
  casadi_assert(x_impl.numel()==dx_impl.numel(),
    "Size of implicit differential state vector (x_impl: " + str(x_impl.numel()) + ") must "
    "match size of implicit differential state derivative vector (dx_impl: "
    + str(dx_impl.numel()) + ").");

  X z   = get_from_dict(dae, "z",   X(0, 1));
  X alg = get_from_dict(dae, "alg", X(0, 1));
  casadi_assert(z.numel()+x_impl.numel()==alg.numel(),
    "Size of algebraic state vector (z: " + str(z.numel()) + ") + "
    "size of implicit states (x_impl: " + str(x_impl.numel()) + ") must "
    "match size of algebraic equations (alg: " + str(alg.numel()) + ").");

  X t = get_from_dict(dae, "t", X::sym("t"));
  casadi_assert(t.is_scalar(), "Time must be scalar. Got " + t.dim() + " instead.");

  X p = get_from_dict(dae, "p", X(0, 1));

  // Check that dx_impl are classified correctly
  std::vector<bool> dep = X::which_depends(alg, dx_impl);

  if (!all(dep)) {
    std::vector<X> dx_impl_split = vertsplit(dx_impl);
    std::vector<std::string> prob;
    for (casadi_int i=0;i<dep.size();++i) {
      if (!dep[i]) prob.push_back(dx_impl_split[i].name());
    }
    casadi_error("Found dx_impl variables that do not appear in alg: " + join(prob, ",") +
      ". They should be classified as z instead.");
  }


  bool normal_order = true;

  // Determine graph structure of problem for structural index reduction
  X V;
  if (normal_order) {
    V = vertcat(x_impl, dx_impl, z);
  } else {
    V = vertcat(x_impl, z, dx_impl);
  }

  Function temp("temp", {V, p, t}, {alg});
  Sparsity G = temp.jac_sparsity(0, 0);

  // Populate var_map: a list that associates variables with their derivatives
  int nx_impl = x_impl.numel();
  int nz = z.numel();
  std::vector<casadi_int> var_map(V.numel(), -1);
  for (casadi_int i=0;i<nx_impl;++i) {
    var_map[i] = i+nx_impl+(normal_order? 0: nz);
  }

  // Allocate space for eq_map: a list that associates equations with their derivatives
  std::vector<casadi_int> eq_map;

  // Modifies var_map and eq_map
  dae_struct_detect(algorithm, G, var_map, eq_map, max_iter);

  // Variables should not be removed
  casadi_assert_dev(var_map.size()>=2*nx_impl+nz);

  // List of scalarized variables
  std::vector<X> var_ext = vertsplit(V);
  // Allocate extra space for extension due to index reduction
  var_ext.resize(var_map.size());

  // Populate extra space
  for (casadi_int i=0;i<var_map.size();++i) {
    casadi_int i_der = var_map[i];
    // Derivate index is in the extra space?
    if (i_der>=V.numel())
      // Consstruct a new symbol with derived name
      var_ext[i_der] = X::sym("d"+var_ext[i].name());
  }

  // Prepare arguments for jtimes
  // In essence, this is like to the novel x_impl and dx_impl,
  // but a subset may be promoted to x (explicit diff. state)
  std::vector<X> xs; xs.reserve(var_map.size());
  std::vector<X> dxs; dxs.reserve(var_map.size());

  for (casadi_int i=0;i<var_map.size();++i) {
    casadi_int i_der = var_map[i];
    if (i_der>=0) {
      xs.push_back(var_ext[i]);
      dxs.push_back(var_ext[i_der]);
    }
  }

  // Stack to be used as jtimes arguments
  X x_cat = vertcat(xs);
  X dx_cat = vertcat(dxs);

  // Don't forget to differentiate time itself
  X j1 = vertcat(x_cat, t);
  X j2 = vertcat(dx_cat, 1);

  // List of scalarized equations
  std::vector<X> eq_ext = vertsplit(alg);
  // Allocate extra space for extension due to index reduction
  eq_ext.resize(eq_map.size());

  // For each equation, figure out how many times it was differentiated (=order)
  std::vector<casadi_int> order = get_orders(eq_map);
  casadi_int max_order = *std::max_element(order.begin(), order.end());

  // We lump together all relevant equations before doing differentiation in a target group,
  // in order to make use of shared subexpressions.
  // order -> all equations that give rise to an equation of this order
  std::map<casadi_int, std::vector<casadi_int> > targets;
  for (casadi_int i=0;i<eq_map.size();++i) {
    if (eq_map[i]>=0) {
      targets[order[eq_map[i]]].push_back(i);
    }
  }

  // Non-active equations will be dropped from the final set
  std::vector<bool> active(eq_map.size(), true);

  // Equations that get dropped, but should still hold for initial guess
  std::vector<X> dropped_constraints;

  // Loop over equation elements in increasing order
  for (casadi_int ord=1;ord<=max_order;++ord) {
    const std::vector<casadi_int> & target = targets[ord];

    // Remember to-be-dropped equations
    for (casadi_int i : target) dropped_constraints.push_back(eq_ext[i]);

    // Gather target's equations
    X block = vertcat(vector_slice(eq_ext, target));
    // Perform differentiation
    std::vector<X> dblocks = vertsplit(jtimes(block, j1, j2, false));
    // Scatter computed derivatives into target resultant equation slots
    for (casadi_int i=0;i<target.size();++i) {
      casadi_int t = target[i];
      eq_ext[eq_map[t]] = dblocks[i];
      active[t] = false;
    }
  }

  std::vector<casadi_int> eq_map_inverse = get_inverse(eq_map);

  // Construct the new list of algebraic equations
  // In theory, just collect all active members of eq_ext
  // In practice, we add some Baumgarte stabilization to avoid constraint drifting
  std::vector<X> alg_new;
  for (casadi_int i=0;i<eq_map.size();++i) {
    if (active[i]) {

      // Start from the pure equation
      X new_eq = eq_ext[i];

      // Construct coefficients of the polynomial (x-gamma)^order
      Polynomial p(1);
      for (casadi_int k=0;k<order[i];++k) {
        p *= Polynomial(1, -baumgarte_pole_);
      }

      const std::vector<double>& coeff = p.coeff();

      // Find indices of originating equations of equation i
      // e.g. if eq_i = der(eq_j), eq_j = der(eq_k), obtain [j,k]
      std::vector<casadi_int> gen_i = path(eq_map_inverse, i);

      // Add inner product of coefficients with originating equations of new_eq
      for (casadi_int k=0;k<gen_i.size();++k) {
        new_eq += coeff[k+1]*eq_ext[gen_i[k]];
      }

      // Ready to add
      alg_new.push_back(new_eq);
    }
  }

  std::vector<casadi_int> var_order = get_orders(var_map);
  std::vector<casadi_int> var_map_inverse = get_inverse(var_map);


  // Which implicit states become explicit pure integrators?
  std::vector<casadi_int> impl_to_expl;
  for (casadi_int i=0;i<var_map.size();++i) {
    if (var_order[i]>=2) {
      // Find indices of originating variables of variable i
      // e.g. if var_i = der(var_j), var_j = der(var_k), obtain [j,k]
      std::vector<casadi_int> gen_i = path(var_map_inverse, i);

      // Anything too deep should become explicit integrator
      if (gen_i.size()>1)
        impl_to_expl.insert(impl_to_expl.end(), gen_i.begin()+1, gen_i.end());

    }
  }

  // Take complement
  std::vector<casadi_int> impl_to_expl_compl = complement(impl_to_expl, xs.size());
  std::map<std::string, X> dae_result;

  dae_result["x"] = vertcat(vector_slice(xs, impl_to_expl));
  dae_result["ode"] = vertcat(vector_slice(dxs, impl_to_expl));

  dae_result["x_impl"] = vertcat(vector_slice(xs, impl_to_expl_compl));
  dae_result["dx_impl"] = vertcat(vector_slice(dxs, impl_to_expl_compl));

  dae_result["t"] = t;
  dae_result["p"] = p;

  // Which algebraic variables were not promoted to differential states?
  std::vector<X> z_orig = vertsplit(z);
  std::vector<X> z_new;
  for (casadi_int i=0;i<var_map.size();++i) {
    // Algebraic
    if (var_order[i]==0 && var_map[i]==-1) {
      z_new.push_back(z_orig.at(i-(normal_order?2:1)*nx_impl));
    }
  }
  dae_result["z"] = vertcat(z_new);
  dae_result["alg"] = vertcat(alg_new);

  dae_result["I"] = vertcat(dropped_constraints);

  stats["index"] = max_order+1;

  return dae_result;
}

MXDict dae_reduce_index(const MXDict& dae, Dict& stats, const Dict& opts) {
  return reduce_index_gen(dae, stats, opts);
}

SXDict dae_reduce_index(const SXDict& dae, Dict& stats, const Dict& opts) {
  return reduce_index_gen(dae, stats, opts);
}

template<class X>
std::map<std::string, X> map_semi_expl(const std::map<std::string, X>& dae,
  const std::map<std::string, X>& dae_red,
  Function& state_to_orig, Function& phi) {

  std::map<std::string, X> dae_se;
  dae_se["x"] = vertcat(dae_red.at("x"), dae_red.at("x_impl"));
  dae_se["z"] = vertcat(dae_red.at("dx_impl"), dae_red.at("z"));
  dae_se["ode"] = vertcat(dae_red.at("ode"), dae_red.at("dx_impl"));
  dae_se["alg"] = dae_red.at("alg");
  dae_se["t"] = dae_red.at("t");
  dae_se["p"] = dae_red.at("p");

  X x_impl  = get_from_dict(dae, "x_impl",  X(0, 1));
  X dx_impl = get_from_dict(dae, "dx_impl", X(0, 1));
  X z       = get_from_dict(dae, "z", X(0, 1));

  Sparsity var_map_sp = jacobian(vertcat(dae_red.at("x"), dae_red.at("x_impl"),
    dae_red.at("dx_impl"), dae_red.at("z")), vertcat(x_impl, dx_impl, z)).sparsity();
  DM var_map(var_map_sp, 1.0);

  state_to_orig = Function("state_to_orig",
    {dae_se["x"], dae_se["z"]},
    {x_impl, dx_impl, z},
    {"xf", "zf"},
    {"x_impl", "dx_impl", "z"});

  phi = Function("phi",
    {dae_se["x"], dae_se["z"], dae_se["p"], dae_se["t"]},
    {dae_red.at("I")},
    {"x", "z", "p", "t"},
    {"I"});

  return dae_se;
}

template<class X>
Function init_gen(const std::map<std::string, X>& dae,
  const std::map<std::string, X>& dae_red,
  const std::string& init_solver, const DMDict& init_strength, const Dict& init_solver_options) {

  X x_impl  = get_from_dict(dae, "x_impl",  X(0, 1));
  X dx_impl = get_from_dict(dae, "dx_impl", X(0, 1));
  X z       = get_from_dict(dae, "z", X(0, 1));
  X p_orig       = get_from_dict(dae, "p", X(0, 1));

  Sparsity var_map_sp = jacobian(
    vertcat(dae_red.at("x"), dae_red.at("x_impl"), dae_red.at("dx_impl"), dae_red.at("z")),
    vertcat(x_impl, dx_impl, z)).sparsity();
  DM var_map(var_map_sp, 1.0);

  // NLP to obtain consistent initial guesses
  std::map<std::string, X> nlp;

  MX pmx = MX::sym("pmx", p_orig.sparsity());
  MX tmx = MX::sym("t");

  X init_x_impl_orig  = X::sym("x_impl_init",  x_impl.sparsity());
  X init_dx_impl_orig = X::sym("dx_impl_init", dx_impl.sparsity());
  X init_z_orig       = X::sym("z_init",       z.sparsity());

  MX init_x_impl_orig_mx  = MX::sym("x_impl_init",  x_impl.sparsity());
  MX init_dx_impl_orig_mx = MX::sym("dx_impl_init", dx_impl.sparsity());
  MX init_z_orig_mx       = MX::sym("z_init",       z.sparsity());

  X init_orig = vertcat(init_x_impl_orig, init_dx_impl_orig, init_z_orig);
  MX init_orig_mx = vertcat(init_x_impl_orig_mx, init_dx_impl_orig_mx, init_z_orig_mx);

  X init_sym_red = mtimes(X(var_map), init_orig);
  nlp["p"] = vertcat(init_orig, dae_red.at("p"), dae_red.at("t"));
  MX p = vertcat(init_orig_mx, pmx, tmx);

  // Dae input sanitization

  for (const auto& e : init_strength) {
    casadi_assert(e.second.is_column() && e.second.is_dense(),
      "Dense column vector expected for key '" + e.first + "' of init_strength. "
      "Got " + e.second.dim(true) + " instead.");
  }

  DM init_strength_x_impl = get_from_dict(init_strength, "x_impl", DM::zeros(x_impl.numel(), 1));
  casadi_assert(x_impl.numel()==init_strength_x_impl.numel(),
    "init_strength 'x_impl' entry should have length " + str(x_impl.numel()) + ". "
    "Got length " + str(init_strength_x_impl.numel()) + " instead.");

  DM init_strength_dx_impl = get_from_dict(init_strength, "dx_impl",
    DM::zeros(dx_impl.numel(), 1));
  casadi_assert(dx_impl.numel()==init_strength_dx_impl.numel(),
    "init_strength 'dx_impl' entry should have length " + str(dx_impl.numel()) + ". "
    "Got length " + str(init_strength_dx_impl.numel()) + " instead.");

  DM init_strength_z = get_from_dict(init_strength, "z", DM::zeros(z.numel(), 1));
  casadi_assert(z.numel()==init_strength_z.numel(),
    "init_strength 'z' entry should have length " + str(z.numel()) + ". "
    "Got length " + str(init_strength_z.numel()) + " instead.");

  DM init_strength_cat = vertcat(init_strength_x_impl, init_strength_dx_impl, init_strength_z);

  DM weights = mtimes(var_map, init_strength_cat);
  // Should be written more cleanly
  std::vector<casadi_int> ind = sparsify(weights>0).sparsity().T().get_col();

  std::vector<X> x_parts = {dae_red.at("x"), dae_red.at("x_impl"),
                      dae_red.at("dx_impl"), dae_red.at("z")};

  nlp["x"] = vertcat(x_parts);
  nlp["g"] = vertcat(dae_red.at("I"), dae_red.at("alg"));
  nlp["f"] = sumsqr((nlp["x"](ind)-init_sym_red(ind))*X(weights(ind)));


  MX x0 = mtimes(var_map, init_orig_mx);
  MX lbx = -MX::inf(nlp.at("x").numel());
  MX ubx = MX::inf(lbx.size());

  // Should be written more cleanly
  ind = sparsify(mtimes(var_map, init_strength_cat)==-1).sparsity().T().get_col();
  lbx(ind) = ubx(ind) = x0(ind);

  Function solver = nlpsol("init_solver", init_solver, nlp, init_solver_options);
  MXDict res = solver(MXDict{{"x0", x0}, {"lbg", 0}, {"ubg", 0},
    {"lbx", lbx}, {"ubx", ubx}, {"p", p}});

  std::vector<casadi_int> x_parts_offsets = {0};
  for (casadi_int i=0;i<x_parts.size();i+=2) {
    x_parts_offsets.push_back(x_parts_offsets.back() + x_parts[i].numel() + x_parts[i+1].numel());
  }

  std::vector<MX> parts = vertsplit(res["x"], x_parts_offsets);

  return Function("init_gen",
    {init_x_impl_orig_mx, init_dx_impl_orig_mx, init_z_orig_mx, pmx, tmx},
    parts,
    {"x_impl", "dx_impl", "z", "p", "t"},
    {"x0", "z0"});
}

MXDict dae_map_semi_expl(const MXDict& dae, const MXDict& dae_red,
  Function& state_to_orig, Function& phi) {
    return map_semi_expl(dae, dae_red, state_to_orig, phi);
}

SXDict dae_map_semi_expl(const SXDict& dae, const SXDict& dae_red,
  Function& state_to_orig, Function& phi) {
    return map_semi_expl(dae, dae_red, state_to_orig, phi);
}

Function dae_init_gen(const MXDict& dae, const MXDict& dae_red,
  const std::string& init_solver, const DMDict& init_strength, const Dict& init_solver_options) {
    return init_gen(dae, dae_red, init_solver, init_strength, init_solver_options);
}

Function dae_init_gen(const SXDict& dae, const SXDict& dae_red,
  const std::string& init_solver, const DMDict& init_strength, const Dict& init_solver_options) {
    return init_gen(dae, dae_red, init_solver, init_strength, init_solver_options);
}

} // namespace casadi
