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

#include "ocp_tools.hpp"
#include "casadi/sx/sx_tools.hpp"
#include "casadi/fx/sx_function.hpp"
#include "casadi/stl_vector_tools.hpp"

#include <algorithm>
#include <set>
#include <queue>
#include <stack>

#ifdef HAVE_LAPACK
extern "C" void dsterf_(int *n, double *d, double *e, int *info); // needed by LG_matrices
#endif // HAVE_LAPACK

using namespace std;

namespace CasADi{
  namespace OptimalControl{

// Get the coefficeints for the collocation and continuity equations
void get_collocation_coeff(int K, vector<vector<double> >& C, vector<double>& D, CollocationPoints cp){
  D.resize(K+1);
  C.resize(K+1);
  
  // Collocation point
  SX tau("tau");
  
  // Roots
  double *tau_root = collocation_points[cp][K];

  // Lagrange polynomials
  vector<SXFunction> l(K+1);
  for(int j=0; j<=K; ++j){
    SX L = 1;
    for(int k=0; k<=K; ++k)
      if(k != j)
        L *= (tau-tau_root[k])/(tau_root[j]-tau_root[k]);
  
    l[j] = SXFunction(tau,L);
    stringstream ss;
    ss << "l(" << j << ")";
    l[j].setOption("name",ss.str());    
    l[j].init();
    cout << l[j] << endl;
  }
  
  // Get the coefficients of the continuity equation
  for(int j=0; j<=K; ++j){
    l[j].setInput(1.0);
    l[j].evaluate();
    l[j].getOutput(&D[j]);
  }

  // Get the coefficients of the collocation equation
  for(int j=0; j<=K; ++j){
    C[j].resize(K+1);
    for(int k=0; k<=K; ++k){
      l[j].setInput(&tau_root[k]);
      l[j].setFwdSeed(1.0);
      l[j].evaluate(1,0);
      l[j].getFwdSens(&C[j][k]);
    }
  }
}

// Collocate a variable (one variable per finite element)
void collocate(const SXMatrix& var, vector< SXMatrix >& VAR, int N){
  VAR.resize(N);

  // Loop over the finite elements
  for(int i=0; i<N; ++i){
    VAR[i] = SXMatrix(var.size1(),var.size2());
    for(int r=0; r<var.numel(); ++r){
      // Create state with appropriate name
      stringstream ss;
      ss << var(r) << "_" << i;
      VAR[i](r) = SX(ss.str());
    }
  }
}

// Collocate a variable (K variables per finite element)
void collocate(const SXMatrix& var, vector< vector< SXMatrix > >& VAR, int N, int K){
  VAR.resize(N);
  
  // Loop over the finite elements
  for(int i=0; i<N; ++i){
    VAR[i].resize(K+1);
    for(int j=0; j<=K; ++j){
      VAR[i][j] = SXMatrix(var.size1(),var.size2());
      for(int r=0; r<var.numel(); ++r){
        
        // Create state with appropriate name
        stringstream ss;
        ss << var(r) << "_" << i << "_" << j;
        VAR[i][j](r) = SX(ss.str());
      }
    }
  }
}
  
// Collocate a variable (K variables per finite element)
void collocate_final(const SXMatrix& var, SXMatrix &VARF){
  VARF = SXMatrix(var.size1(),var.size2());
  for(int i=0; i<var.numel(); ++i){
    stringstream ss;
    ss << var(i) << "_f";
    VARF(i) = SX(ss.str());
  }
}


class EquationSorter{
  public:
    EquationSorter(const std::vector<SX>& x, const std::vector<SX>& xdot, const std::vector<SX>& z, const std::vector<SX>& dae);
    ~EquationSorter(){}
    
    // Tarjan's algorithm
    void tarjan();
    void strongconnect(int v, int i);
    
    // Depth-first index
    int index_;
    
    // Index stack
    std::stack<int> S_;
    
    // Elements in stack
    std::vector<bool> instack_[2];
    
    // Depth-first index for all vertices
    std::vector<int> number_[2];
    
    // Index to some node reachable from the node
    std::vector<int> lowlink_[2];
    
    // Index ordering
    std::vector<int> scc_order_[2];
    
    // Execute Cellier's version of Tarjan's algorithm (returns true if algorithm finished successfully)
    void tarjanCellier();
    
    // Color edge (i,j) red
    void colorRed(int i, int j, int eq_no);
  
    // Color edge (i,j) blue
    void colorBlue(int i, int j);
    
    // Color edge 
    void decrease(int v,int i);

  protected:
    std::vector<SX> x_;
    std::vector<SX> xdot_;
    std::vector<SX> z_;
    std::vector<SX> dae_;
    
    // Sparsity pattern
    CRSSparsity I_[2];
    
    // Transpose of sparsity pattern and mapping
    std::vector<int> mapping_;
    
    // Occurence of each equation and variable
    std::vector<int> count_[2];
    
    // Vertices with no uncolored edge
    std::queue<int> singular_[2];
    
    // Vertices with exactly one uncolored edge
    std::queue<int> causal_[2];

    // The order of the vertex (-1 if not assigned)
    std::vector<int> order_[2];
    
    // The maximum and minumum assigned equation
    int max_eq_, min_eq_;
    
    // The equation associated with each variable (-1 if not assigned)
    std::vector<int> assign_;
    
    // Number of equations not yet assigned
    int num_unassigned_;
    
    // Edges
    std::vector<int> edges_;
        
};

EquationSorter::EquationSorter(const std::vector<SX>& x, const std::vector<SX>& xdot, const std::vector<SX>& z, const std::vector<SX>& dae) : x_(x), xdot_(xdot), z_(z), dae_(dae){
  vector<SX> v;
  v.insert(v.end(),xdot_.begin(),xdot_.end());
  v.insert(v.end(),z_.begin(),z_.end());
  
  SXFunction fcn(v,dae_);
  fcn.init();
  SXMatrix jac = fcn.jac(0,0);
  
  jac.sparsityRef().strongly_connected_components();
  
  // Get sparsity pattern (incidence matrix)
  I_[0] = jac.sparsity();
  
  // Get transpose of sparsity pattern
  I_[1] = I_[0].transpose(mapping_);

  // No elements in stack
  instack_[0].resize(I_[0].size1(),false);
  instack_[1].resize(I_[0].size2(),false);
  
  // Initialize indices to -1 (unassigned)
  number_[0].resize(I_[0].size1(),-1);
  number_[1].resize(I_[0].size2(),-1);
  
  // Index to some node reachable from the node
  lowlink_[0].resize(I_[0].size1(),-1);
  lowlink_[1].resize(I_[0].size2(),-1);
  
  scc_order_[0].reserve(I_[0].size1());
  scc_order_[1].reserve(I_[0].size2());
  
  // Depth-first index
  index_ = 0;

  // Count the number of occurences of each equation and variable
  count_[0].resize(I_[0].size1(),0);
  count_[1].resize(I_[0].size2(),0);
  
  // Loop over rows
  for(int i=0; i<I_[0].size1(); ++i){
    // Loop over non-zero elements
    for(int el=I_[0].rowind(i); el<I_[0].rowind(i+1); ++el){
      // column
      int j = I_[0].col(el);
      
      // Increase counters
      count_[0][i]++;
      count_[1][j]++;
    }
  }

  // Add equations to queue
  for(int i=0; i<count_[0].size(); ++i){
    if(count_[0][i]==0)
      singular_[0].push(i);
    else if(count_[0][i]==1)
      causal_[0].push(i);
  }
    
  // Add variables to queue
  for(int j=0; j<count_[1].size(); ++j){
    if(count_[1][j]==0)
      singular_[1].push(j);
    else if(count_[1][j]==1){
      causal_[1].push(j);
    }
  }
  
  // Equation assignment vector
  assign_.resize(I_[0].size2(),-1);
  num_unassigned_ = assign_.size();

  // Allocate edges
  edges_.resize(I_[0].size(),0);
  
  // Number of equations
  int neq = I_[0].size1();
  
  // Number of variables
  int nvar = I_[0].size2();
  
  order_[0].resize(neq,-1);
  order_[1].resize(nvar,-1);
  min_eq_ = -1;
  max_eq_ = neq;
  
}

void EquationSorter::colorBlue(int i, int j){
  decrease(0,i);
  decrease(1,j);
}
  
void EquationSorter::decrease(int v, int i){
  count_[v][i]--;
  if(order_[v][i]==-1){
    if(count_[v][i]==0)
      singular_[v].push(i);
    else if(count_[v][i]==1)
      causal_[v].push(i);
  }
}
      
void EquationSorter::colorRed(int i, int j, int eq_no){
  cout << "coloring (" << i << "," << j << "), equation " << eq_no << endl;
  
  // Number the equation and variable
  order_[1][j] = order_[0][i] = eq_no;
  
  // Decrease number of unassigned
  num_unassigned_--;
  
  // Loop over the variables associated with equation i
  for(int el=I_[0].rowind(i); el<I_[0].rowind(i+1); ++el){
    // Variable
    int jj=I_[0].col(el);
    
    // If variable matches
    if(j==jj){
      // Make sure edge not colored
      casadi_assert(edges_[el]==0);
      
      // Color edge (i,j)
      edges_[el] = 1;
      
      // Decrease number of uncolored variables and equations
      count_[0][i]--;
      count_[1][jj]--;
      
      // If doesn't match and edge is uncolored
    } else if(edges_[el]==0){
      
      // Color edge
      edges_[el] = 2;
      
      // Decrease number of uncolored variables and equations
      decrease(0,i);
      decrease(1,jj);
    }
  }
  
  // Loop over the equations associated with variable j
  for(int el=I_[1].rowind(j); el<I_[1].rowind(j+1); ++el){
    
    // Variable
    int ii=I_[1].col(el);

    // If edge is uncolored
    if(edges_[mapping_[el]]==0){
      
      // Color edge
      edges_[mapping_[el]] = 2;
      
      // Decrease number of uncolored variables and equations
      decrease(0,ii);
      decrease(1,j);
    }
  }
  
  
  
/*  // Assign
  assign_[j]=i;
  num_unassigned_--;
  
  // Color edge
  edges_[el]=1;
  
  // Decrease counters
  count_[0][i]--;
  count_[1][j]--;
  
  // Terminate search
  break;*/
}
  
void EquationSorter::tarjan(){
  // Start depth-first searches for all equations that have not been visited
  for(int i=0; i<I_[0].size1(); ++i){
    if(number_[0][i]<0){
      strongconnect(0,i);
    }
  }
  
  // Start depth-first searches for all variables that have not been visited
  for(int i=0; i<I_[1].size1(); ++i){
    if(number_[1][i]<0){
      strongconnect(1,i);
    }
  }
  
  //
  cout << scc_order_[0] << endl;
  cout << scc_order_[1] << endl;
  
  
  vector<SX> v;
  v.insert(v.end(),xdot_.begin(),xdot_.end());
  v.insert(v.end(),z_.begin(),z_.end());
  
  vector<SX> v_sorted(v.size());
  for(int i=0; i<v.size(); ++i)
    v_sorted[i] = v[scc_order_[1][i]];
  
  
  vector<SX> dae_sorted(dae_.size());
  for(int i=0; i<dae_.size(); ++i)
    dae_sorted[i] = dae_[scc_order_[0][i]];
  
  
  SXFunction fcn(v_sorted,dae_sorted);
  fcn.init();
  SXMatrix jac = fcn.jac(0,0);
  for(int k=0; k<jac.size(); ++k)
    jac.at(k) = 1;
  jac.printDense();

  
}

void EquationSorter::strongconnect(int v, int i){
  number_[v][i] = index_;
  lowlink_[v][i] = index_;
  index_++;
  S_.push(v + 2*i);
  instack_[v][i] = true;
  
  // Loop over edges going out from equation i
  for(int el=I_[v].rowind(i); el<I_[v].rowind(i+1); ++el){
    // Variable
    int j=I_[v].col(el);
    
    // Add to stack if index undefined
    if(number_[!v][j]<0){
      strongconnect(!v,j);
      lowlink_[v][i] = std::min(lowlink_[v][i],lowlink_[!v][j]);
    } else if(instack_[!v][j] && number_[!v][j] < number_[v][i]){
      lowlink_[v][i] = std::min(lowlink_[v][i],number_[!v][j]);
    }
  }
   
  // If v is a root node, pop the stack and generate an SCC
  if(number_[v][i] == lowlink_[v][i]){
    cout << "start scc" << endl;
    while(true){
      int w = S_.top();
      int wv = w % 2;
      int wi = w / 2;
      
      // Break loop?
      if(number_[wv][wi] < number_[v][i]) break;
      
      // Delete from stack
      S_.pop();
      scc_order_[wv].push_back(wi);
      instack_[wv][wi] = false;
      
      cout << "added " << wv << ", " << wi << endl;
    }
    
    cout << "end scc" << endl;
    
  }
}
    
void EquationSorter::tarjanCellier(){
  cout << "A" << count_[0] <<  ", "<< count_[1] << endl;

  
  // Process the queues
  while(num_unassigned_>0){

    // Make sure that no equation has become singular
    casadi_assert(singular_[0].empty());

    // Make sure that no variable has become singular
    casadi_assert(singular_[1].empty());

    // Check if there are equations with only one uncolored edge
    if(!causal_[0].empty()){
      // Get the equation from queue
      int i=causal_[0].front();
      causal_[0].pop();
      
      // Loop over edges
      for(int el=I_[0].rowind(i); el<I_[0].rowind(i+1); ++el){
        
        // Variable
        int j=I_[0].col(el);
        
        // If the edge is uncolored
        if(edges_[el]==0){
          
          // Color edge (i,j) red
          colorRed(i,j, ++min_eq_);
          break;
        }
      }
      continue;
    }
    
    // Check if there are variables with only one uncolored edge
    if(!causal_[1].empty()){
      
      // Get the variable from queue
      int j=causal_[1].front();
      causal_[1].pop();
      
      // Loop over edges
      for(int el=I_[1].rowind(j); el<I_[1].rowind(j+1); ++el){
        
        // Equation
        int i=I_[1].col(el);
        
        // If the edge is uncolored
        if(edges_[mapping_[el]]==0){
          
          // Color edge (i,j)
          colorRed(i,j, --max_eq_);
          break;
        }
      }
      continue;
    }
    
    break;
  }

  cout << order_[0] << endl;
  cout << order_[1] << endl;
  cout << "B" << count_[0] <<  ", "<< count_[1] << endl;

  if(num_unassigned_==0){
  
    vector<SX> v;
    v.insert(v.end(),xdot_.begin(),xdot_.end());
    v.insert(v.end(),z_.begin(),z_.end());
    
    vector<SX> v_sorted(v.size());
    for(int i=0; i<v.size(); ++i)
      v_sorted[order_[1][i]] = v[i];
    
    
    vector<SX> dae_sorted(dae_.size());
    for(int i=0; i<dae_.size(); ++i)
      dae_sorted[order_[0][i]] = dae_[i];
    
    
    SXFunction fcn(v_sorted,dae_sorted);
    fcn.init();
    SXMatrix jac = fcn.jac(0,0);
    for(int k=0; k<jac.size(); ++k)
      jac.at(k) = 1;
    
    jac.printDense();
  }

}




void tarjan(const std::vector<SX>& x, const std::vector<SX>& xdot, const std::vector<SX>& z, const std::vector<SX>& dae){
  EquationSorter diffeq(x,xdot,z,dae);

  diffeq.tarjan();
  
  

  
  
  
  
  // Transpose of sparsity pattern
/*  vector<int> mapping;
  CRSSparsity IT = I_[0].transpose(mapping);*/
  
  // 
  
  
  
  
/*  I_[0].print();*/
/*  IT.print();*/
  
  
  
  
}

  #if 0
  
void makeExplicit(AcadoOCP &ocp){
  
  vector<ODE> f_new;

  for(vector<ODE>::const_iterator it=ocp.f.begin(); it!=ocp.f.end(); ++it){
    // Try to split the ode in two
    SX rhs, lhs;
    try{
      it->split(lhs,rhs);

      // Make sure that the left hand side is symbolic
      if(lhs->isNan() || !lhs->isSymbolic())
	throw CasadiException("Left hand side not a variable");
    } catch(...){
      // Cannot split
      f_new.push_back(*it);
      continue;
    }
    
//     // Try to locate the right hand side
//     map<std::string, Derivative>::const_iterator itd=ocp.der_map.find(lhs->getName());
//     if(itd!=ocp.der_map.end()){
//       // Add explicit ODE
//       ExplicitODE eode(itd->second,rhs);
//       ocp.f_ex.push_back(eode);
//       continue;
//     }

    map<std::string, int>::iterator itv=ocp.varname.find(lhs->getName());
    if(itv!=ocp.varname.end()){
      // Add dependent equation
      DependentEquation eode(ocp.vars[itv->second],rhs);
      ocp.f_al_ex.push_back(eode);
      
      // Mark the variable dependent
      ocp.vars[itv->second]->type = Variable::DEPENDENT;
      
      continue;
    }
    
    // Default
    f_new.push_back(*it);
    continue;
  }
  ocp.f = f_new;
}
#endif

#if 0
void makeImplicit(AcadoOCP &ocp){
  for(vector<Equation>::iterator it=ocp.eqs.begin(); it!= ocp.eqs.end(); ++it){
    it->makeImplicit();
  }
  
  ocp.sortEquations();
}
  
  
  
#endif
  
  
  
  
  
  
  
  
  
  
  
#if 0
void eliminateTimeDerivatives(OCP_old &ocp){
assert(ocp.numFcn(DE_G) != 0);
  for(int i=0; i<ocp.numFcn(DE_G); ++i)
    replaceDerivatives(ocp.function(DE_G,i).fcn,ocp.getVars(VAR_X),ocp.getVars(VAR_XDOT));

}

void makeSmooth(OCP_old& ocp){
  Matrix bvar, bexp;

  makeSmooth(ocp.function(DE_F).fcn,bvar,bexp);
  makeSmooth(ocp.function(DE_XDOT).fcn,bvar,bexp);

  ocp.addVariable(bvar, FREE, CONTINUOUS, BINARY);

  // Save zero-crossing functions
  ocp.addToFunction(SWITCH, bexp);
}

void eliminateLagrangeTerm(OCP_old& ocp){
  if(!ocp.hasFcn(LTERM)) return;
  
  // Introduce a lagrangian term
  Matrix l("lterm"), ldot("lterm_dot");

  // Create a new state
  ocp.addState(l,ldot);

  // Add differential equation
  ocp.addEquation(ldot, ocp.function(LTERM).fcn);

  // Add initial constraint
  ocp.addConstraint(ocp.function(TPOINT,0).fcn, l, 0, 0);
  
  // Create a Mayer term
  ocp.addMayerTerm(l);

  // Delete the lagrange term
  ocp.function(LTERM).fcn.clear();
}

void sortInitialConditions(OCP_old& ocp, bool make_explicit){
  // Quick return if no implicitly defined initial conditions
  if(!ocp.hasFcn(IC_IMPFUN) || ocp.function(IC_IMPFUN).fcn.empty()) return;

  // All state variables
  Matrix eivar = ocp.getVars(VAR_X);

  // Make sure that ivar is allocated by adding an empty expression
  ocp.addToFunction(IC_IMPVAR,Matrix());

  // Referances to the implicit variables
  Matrix& ifun = ocp.function(IC_IMPFUN).fcn;
  Matrix& ivar = ocp.function(IC_IMPVAR).fcn;

  // Eliminate explicitly known variables and calculate the ivar
  if(!ocp.hasFcn(IC_EXPFUN) || ocp.function(IC_EXPFUN).fcn.empty()){
    ivar = eivar;
  } else {
    // References to the function and variable
    Matrix& efun = ocp.function(IC_EXPFUN).fcn;
    Matrix& evar = ocp.function(IC_EXPVAR).fcn;

    // Eliminate the explicitly known initial conditions from the implicitly defined equations
    substitute(ifun,evar,efun);

    // nodes of the variables that are known explicitly
    set<SXNode*> expvars; 
    for(int i=0; i<evar.size(); ++i)
      expvars.insert(evar[i].get());
    
    // ivar are all other state variables
    ivar.clear();
    for(int i=0; i<eivar.size(); ++i)
      if(expvars.count(eivar[i].get()) == 0)
	ivar << eivar[i];
  }
  
  if(make_explicit){
    // Differentiate ifun with respect to ivar
    Matrix A = gradient(ifun,ivar);

    // Make sure that A is invertable
    if(dependsOn(A,ivar)) throw "OCP::sortInitialConditions failed: IC not affine in state";
    if(isZero(det(A))) throw "OCP::sortInitialConditions failed: A not invertable";

    // Write the differential equation in explicit form
    Matrix rhs = solve(A,A*ivar-ifun);
    simplify(rhs);

    // Add to the explicit ic's
    ocp.addToFunction(IC_EXPVAR,ivar);
    ocp.addToFunction(IC_EXPFUN,rhs);
  
    // Delete the implicit ic's
    ifun.clear();
    ivar.clear();
  }
}
// void discretizeControls(OCP_old& ocp, int n_disc){
//   int n = ocp.ntp();
//   assert(n==2);
//   Matrix t_disc = linspace(ocp.getFcn(TPOINT,0),ocp.getFcn(TPOINT,n-1),n_disc+1);
//   for(int i=1; i<n_disc; ++i)
//     ocp.addTimePoint(t_disc[i],i);
// 
//   // Make all controls discrete (may only change values at substages)
//   for(int i=0; i<variables.size(); ++i){
//     if(ocp.variables[i].id.type == FREE && variables[i].id.dyn == CONTINUOUS)
//       ocp.variables[i].id.dyn = DISCRETE;
//   }
// 
//   // Update the indices
//   vector<int> & cind = ocp.findVars(FREE,CONTINUOUS);
//   vector<int> & dind = ocp.findVars(FREE,DISCRETE);
//   dind.insert(dind.end(),cind.begin(),cind.end()); // append indices to the end
//   cind.clear(); // remove indices to continuous variables
// }


void makeImplicit(OCP_old& ocp){
  ocp.function(DE_F).fcn -= ocp.function(DE_XDOT);
  ocp.function(DE_XDOT).fcn = Matrix(ocp.function(DE_XDOT).size());
}

void makeSemiExplicit(OCP_old& ocp){
  // Start with a fully implicit DAE
  makeImplicit(ocp);

  Matrix xdot = ocp.getVars(VAR_XDOT);

  // Derivate with respect to xdot to get xdot dependence
  Matrix A = gradient(ocp.function(DE_F).fcn,xdot);
  if(dependsOn(A,ocp.getVars(VAR_XDOT))) throw "OCP::makeSemiExplicit failed. The differential equation is not linear in xdot";  
  Matrix Adep = trans(A)*xdot;

  // Subtract from both sides
  ocp.function(DE_XDOT).fcn -= Adep;
  ocp.function(DE_F).fcn -= Adep;

  // Simplify the expressions (remove terms that cancel out)
  simplify(ocp.function(DE_XDOT).fcn);
  simplify(ocp.function(DE_F).fcn);
assert(0);
}

bool isExplicit(const OCP_old &ocp){
  // Check dimensions
 if(ocp.function(DE_XDOT).size() != ocp.getVars(VAR_XDOT).size()) throw "OCP::isExplicit: differential equation is not balanced";

  Matrix xdot = ocp.getVars(VAR_XDOT);

  // Check if any variable is different
  for(int i=0; i<ocp.function(DE_XDOT).size(); ++i)
    if(!ocp.function(DE_XDOT).fcn[i].isEqual(xdot[i])){
      return false;
    }
  return true; 
}

void makeExplicit(OCP_old& ocp){
  Matrix xdot = ocp.getVars(VAR_XDOT);
  Matrix& g = ocp.function(DE_G);

  // Derivate the left hand side with respect to xdot
  Matrix A = jacobian(g,xdot);

  // Make sure that A is invertable
  if(isZero(det(A))) throw "OCP::makeExplicit failed. The differential equation is not index 0 and affine in xdot";

  // Write the differential equation in explicit form
  ocp.addToFunction(DE_XDOT,xdot);
  ocp.addToFunction(DE_F,inv(A)*(A*xdot-g));
  ocp.function(DE_G).fcn.clear();
  simplify(ocp.function(DE_F).fcn);
}

Matrix parametrizeControls(OCP_old& ocp, const Matrix &u, int n_disc, const Matrix& met){
  assert(met.empty());
  
  // Uniform grid
  Matrix t_disc = linspace(ocp.getInitialTime(),ocp.getFinalTime(),n_disc+1);

  // Get the interior values of t_disc
  Matrix t_int;
  getRow(t_int, t_disc, 1, n_disc-1);

  // The points where the guess is evaluated
  Matrix t_guess;
  getRow(t_guess, t_disc, 0, n_disc);

  // Add time points
  assert(ocp.nsub()==1);
  for(int i=1; i<n_disc; ++i)
    ocp.addTimePoint(t_disc[i],i);

  // Locate the controls
  vector<int> inds = ocp.getVarInds(u);

  // Return value
  Matrix u_disc_all;

  // Loop over controls
  for(vector<int>::iterator it = inds.begin(); it != inds.end(); it++){
    // Create a set of parameters
    Matrix u_disc(getName(ocp.variable(*it).var),n_disc);
    
    // Save to the return value
    u_disc_all << trans(u_disc);

    // Add to the ocp
    ocp.addParameter(u_disc);
  
    // Get a guess by evaluating the guess for u at different points
    Matrix u_guess = ocp.variable(*it).solution;
    SXFunction u_guess_rt(ocp.getTime(),u_guess);

    Matrix temp;
    u_guess_rt->eval_symbolic(t_guess,temp);
    ocp.guessSolution(u_disc,temp);

    // Make a piecewise constant approximation of u in time and save to the variable structure
    Matrix u_def = pw_const(ocp.getTime(), t_int, u_disc);

    // from now on, u is a dependent variable
    ocp.variable(*it).id = VAR_DEP; 
    
    // Save the definition
    ocp.variable(*it).solution = u_def;
  }

  // Indices of controls variables
  vector<int> & u_ind = ocp.findVars(VAR_U);

  // Indices of dependent variables
  vector<int> & d_ind = ocp.findVars(VAR_DEP);

  // Save the indices to the list of dependent variables
  d_ind.insert(d_ind.end(),inds.begin(),inds.end());


  // Remove the indices from the vector of control indices (INEFFICIENT - better create a new vector!)
  while(!inds.empty()){
    bool found = false;
    for(int i=0; i<u_ind.size(); ++i)
      if(inds.back() == u_ind[i]){
        u_ind.erase(u_ind.begin()+i);
        inds.pop_back();
        found = true;
        break;
      }
    if(!found) throw "parametrizeControls: Could not find control";
  }
  return trans(u_disc_all);
}

void eventIteration(OCP_old &ocp, bool forward){
  
  double t_;
  ocp.getArg(ocp.findVars(VAR_T),&t_);
  if(forward)  t_ += 1e-3;
  else         t_ -= 1e-3;

  ocp.setArg(ocp.findVars(VAR_T),&t_);

  int nv_ = ocp.numVar(VAR_V);
  if(nv_ == 0) return;

  vector<double> new_v(nv_);
  vector<double> old_v(nv_);


  const int maxiter = 100;
  for(int iter=0; iter<maxiter; ++iter){

    // Save the old value
    ocp.getArg(ocp.findVars(VAR_VD),vecptr(old_v));

    // solve for the new v
    ocp.function(SWITCH).eval(vecptr(new_v));
    for(int i=0; i<nv_; ++i)
      new_v[i] = new_v[i]>=0;
    ocp.setArg(ocp.findVars(VAR_V),vecptr(new_v));

        // Return if all variables are equal
    bool all_equal = true;
    for(int i=0; i<nv_; ++i)
      if(old_v[i] != new_v[i]){
	all_equal = false; 
	break;
      }
    if(all_equal){
      if(forward)  t_ -= 1e-3;
      else         t_ += 1e-3;
      ocp.setArg(ocp.findVars(VAR_T),&t_);
      return;
    }
  }

  throw "OCP::eventIteration: The algorithm failed to converge";
}



void generateLegendreMatrices(int n, Matrix &d_, Matrix &w_, Matrix& tau_){
  double d[n][n+1];
  double w[n];
  double tau[n+2];

  for(int i=0; i<n; ++i)
    tau[i+1] = 0;
  tau[0] = -1;
  tau[n+1] = 1;
  double e[n-1]; // off-diagonal elements
  for(int i=1; i<n; ++i)
    e[i-1] =  double(i)/sqrt(4*i*i-1);

   // find the roots of the Legendre polynomials
   int info;
#ifdef HAVE_LAPACK
   dsterf_(&n,tau+1,e,&info);
#else
  throw "no dsterf_";
#endif
  // calculate the corresponding weights
  for(int i=0; i<n; ++i){
    double x=tau[i+1];
    // evaluate the legendre polynomials up to order n
    double p[n+1];
    p[0] = 1;
    p[1] = x;
    for(int j=1; j<n; ++j)
      p[j+1] = ((2*j+1)*x*p[j] - j*p[j-1])/(j+1);

    // evaluate the p_prim[n]
    double p_prim = (x*p[n]-p[n-1])*n/(x*x-1);

    // evaluate the weights
    w[i] = 2/(1-x*x)/(p_prim*p_prim);
  }
  // Calculate the differentiation matrix D (should be made more efficient!)
  for(int k=1; k<n+1; ++k){
    for(int i=0; i<n+1; ++i){
      d[k-1][i] = 0;    
      double factor = 1;
      for(int j=0; j<n+1; ++j)
        if(j != i && j != k)
          factor *= tau[k] - tau[j];

      for(int j=0; j<n+1; ++j)
        if(j != i)
          factor /= tau[i] - tau[j];

      for(int l=0; l<n+1; ++l){
        if(k==i){
          if(l==i)
            d[k-1][i] += 1; // term forgot in gpops?
          else
            d[k-1][i] += 1/(tau[i]-tau[l]);
        } else if(l!=i && l==k){
          d[k-1][i] += factor;
        }
      }
      if(k==i)
        d[k-1][i] -= 1; // reproduce the bug(?) of gpops
    }
  }

  tau_= Matrix(&tau[0],n+2,1);
  tau_ = (ones(tau_.size())+tau_)/2;
  w_ = Matrix(&w[0],n,1);
  d_ = Matrix(&d[0][0],n,n+1);
  d_.resize(n,n+2); // include one extra column
  d_ = trans(d_);
}

void generateLagrangePolynomials(const Matrix &t, const Matrix &tau, Matrix &L){
  
  // Number of LG points
  int n = tau.size()-2;

  // Carry out the subtraction to get the factors
  Matrix f = t*ones(tau.size()) - tau;

  // Create the forward and backward cumulative products:
  // Forward:   {1, f[0]. f[0]*f[1], ... , f[0]*...*f[i-1]}
  // Backward:  {1, f[n]. f[n]*f[n-1], ... , f[n]*...*f[i+1]}
  Matrix fw = 1, bw = 1;
  fw.reserve(n+1);
  bw.reserve(n+1);
  for(int i=0; i<n; ++i){
    fw << fw[i] * f[i];
    bw << bw[i] * f[n-i];
  }

  // Calculate the constant term
  Matrix c = ones(n+1,1);
  for(int i=0; i<=n; ++i)
    for(int j=0; j<=n; ++j)
      if(j != i)
        c[i] *= tau[i]-tau[j];

  // Calculate the Lagrange polynomial
  L = Matrix(n+2,1); // last row is zero
  for(int i=0; i<=n; ++i)
    L[i] = fw[i]*bw[n-i]/c[i];
}



ostream& operator<<(ostream &stream, const OCP_old& ocp){

   for(int i=0; i<ocp.numVar(); ++i){
      stream << "Variable " << i << ": " << ocp.variable(i) << endl;
    }
    
      if(ocp.hasFcn(DE_F)){
    stream << "--- DYNAMIC EQUATIONS - EXPLICIT --- " << endl;
    const Matrix& de_lhs = ocp.function(DE_XDOT);
    const Matrix& de_rhs = ocp.function(DE_F);
    assert(de_lhs.size() == de_rhs.size());
  
    for(int i=0; i<de_lhs.size(); ++i)  
      stream << de_lhs[i] << " == " << de_rhs[i] << endl;

    // Print an empty line
    stream << endl;
  }

  if(ocp.hasFcn(DE_G)){
    stream << "--- DYNAMIC EQUATIONS - IMPLICIT --- " << endl;
    const Matrix& de_g = ocp.function(DE_G);
  
    for(int i=0; i<de_g.size(); ++i)  
      stream << "0 == " << de_g[i] << endl;

    // Print an empty line
    stream << endl;
  }

  if(ocp.hasFcn(IC_EXPFUN)){
    stream << "--- INITIAL EQUATIONS - EXPLICIT --- " << endl;
    const Matrix& icvar = ocp.function(IC_EXPVAR);
    const Matrix& icfun = ocp.function(IC_EXPFUN);
  
    for(int i=0; i<icfun.size(); ++i)  
      stream << icvar[i] << " == " << icfun[i] << endl;

    // Print an empty line
    stream << endl;
  }

  if(ocp.hasFcn(IC_IMPFUN)){
    stream << "--- INITIAL EQUATIONS - IMPLICIT --- " << endl;
    const Matrix& icfun = ocp.function(IC_IMPFUN);
  
    for(int i=0; i<icfun.size(); ++i)  
      stream << "0 == " << icfun[i] << endl;

    // Print an empty line
    stream << endl;
  }

  
  stream << "--- OBJECTIVE FUNCTION --- " << endl;
  if(ocp.hasFcn(MTERM)){
    stream << "Mayer terms:      " << ocp.function(MTERM).fcn  << endl;
    stream << "Mayer times:      " << ocp.function(MTP).fcn << endl;
  }
  if(ocp.hasFcn(LTERM)){
    stream << "Lagrange terms:   " << ocp.function(LTERM).fcn  << endl;
  }
  if(ocp.hasFcn(LSQ_LTERM)){
    stream << "Least squares Lagrange terms:   " << ocp.function(LSQ_LTERM).fcn  << endl;
    stream << "Least squares Lagrange ref:   " << ocp.function(LSQ_LREF).fcn  << endl;
    stream << "Least squares Lagrange matrix weight:   " << ocp.function(LSQ_LWEIGHT).fcn  << endl;
  }

  stream << endl;

  return stream;
}

void eliminateDependent(OCP_old& ocp){
  cout << "Warning: eliminating dependent variables from the ocp may harm the performance" << endl;

  // get the indices of the dependent variables
  const vector<int>& inds = ocp.findVars(VAR_DEP);

  // get the dependent variables
  Matrix dep = ocp.getVars(inds);
  
  // get the definitions
  Matrix depdef;
  for(vector<int>::const_iterator it = inds.begin(); it != inds.end(); ++it)
    depdef << ocp.variable(*it).solution;

  // eliminate from all the functions
  for(int i=0; i<ocp.numFcn(); ++i){
      substitute(ocp.function(i).fcn,dep,depdef);
  }
}
#endif

} // namespace OptimalControl
} // namespace CasADi

