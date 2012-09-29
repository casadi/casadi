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
#include "variable_tools.hpp"

#include <algorithm>
#include <set>
#include <queue>
#include <stack>

#ifdef HAVE_LAPACK
extern "C" void dsterf_(int *n, double *d, double *e, int *info); // needed by LG_matrices
#endif // HAVE_LAPACK

using namespace std;

namespace CasADi{
  
void updateDependent(SymbolicOCP& ocp){
  // Quick return if no dependent parameters
  if(ocp.pd.empty()) return;
  
  // Get the binding equations
  SXMatrix pd_def = binding(ocp.pd);
  
  // Get expressions for the variables
  SXMatrix pd = var(ocp.pd);
  
  // Sort out interdependencies
  bool reverse=false;
  substituteInPlace(pd, pd_def, reverse);

  // Create a function which evaluates the binding equations numerically
  SXMatrix ci = var(ocp.ci);
  SXMatrix pi = var(ocp.pi);
  vector<SXMatrix> f_in(2);
  f_in[0] = ci;
  f_in[1] = pi;
  SXFunction f(f_in,pd_def);
  f.init();
  
  // Evaluate the start attribute
  f.setInput(getStart(ocp.ci),0);
  f.setInput(getStart(ocp.pi),1);
  f.evaluate();
  const vector<double>& res = f.output().data();
  
  // Save to variables
  for(int k=0; k<ocp.pd.size(); ++k){
    ocp.pd[k].setStart(res[k]);
  }
}
  

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
  SXMatrix jac = fcn.jac();
  
  casadi_assert(0);
  //jac.sparsity().stronglyConnectedComponents();
  
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
  SXMatrix jac = fcn.jac();
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
    SXMatrix jac = fcn.jac();
    for(int k=0; k<jac.size(); ++k)
      jac.at(k) = 1;
    
    jac.printDense();
  }

}

} // namespace CasADi

