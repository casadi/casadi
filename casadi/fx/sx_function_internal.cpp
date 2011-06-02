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

#include "sx_function_internal.hpp"
#include <cassert>
#include <limits>
#include <stack>
#include <deque>
#include <fstream>
#include <sstream>
#include "../stl_vector_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../sx/sx_node.hpp"
#include "../mx/evaluation.hpp"

namespace CasADi{

using namespace std;


SXFunctionInternal::SXFunctionInternal(const vector<Matrix<SX> >& inputv_, const vector<Matrix<SX> >& outputv_) : inputv(inputv_),  outputv(outputv_) {
  addOption("ad_mode",OT_STRING,"reverse");
  addOption("symbolic_jacobian",OT_BOOLEAN,true); // generate jacobian symbolically by source code transformation
  setOption("name","unnamed_sx_function");
  
  casadi_assert(!outputv.empty());

  // Check that inputs are symbolic
  for(int i=0; i<inputv_.size(); ++i) {
    if (!isSymbolicSparse(inputv_[i])) {
      stringstream ss;
      ss << "SXFunctionInternal::SXFunctionInternal: SXfunction input arguments must be purely symbolic." << endl;
      ss << "Argument #" << i << " is not symbolic." << endl;
      throw CasadiException(ss.str());
    }
  }

  // Stack
  stack<SXNode*> s;

  // Add the inputs to the stack
  for(vector<Matrix<SX> >::const_iterator it = inputv.begin(); it != inputv.end(); ++it)
    for(vector<SX>::const_iterator itc = it->begin(); itc != it->end(); ++itc)
      s.push(itc->get());

  // Add the outputs to the stack
  for(vector<Matrix<SX> >::const_iterator it = outputv.begin(); it != outputv.end(); ++it)
    for(vector<SX>::const_iterator itc = it->begin(); itc != it->end(); ++itc)
      s.push(itc->get());

  // Order the nodes in the order of dependencies using a depth-first topological sorting
  nodes.clear();
  sort_depth_first(s,nodes);

  // Sort the nodes by type
  vector<SXNode*> bnodes, cnodes, snodes;
  for(vector<SXNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
    SXNode* t = *it;
    if(t->isConstant())
      cnodes.push_back(t);
    else if(t->isSymbolic())
      snodes.push_back(t);
    else
      bnodes.push_back(t);
  }
  
  // Resort the nodes in a more cache friendly order (Kahn 1962)
  resort_bredth_first(bnodes);

  // Set the temporary variables to be the corresponding place in the sorted graph
  for(int i=0; i<nodes.size(); ++i)
    nodes[i]->temp = i;

  // Place in the work vector for each of the nodes in the tree
  vector<int> place(nodes.size(),-1);

  worksize = nodes.size();
  for(int i=0; i<place.size(); ++i)
    place[i] = i;


  // Add the binary operations
  algorithm.clear();
  algorithm.resize(bnodes.size());
  pder1.resize(algorithm.size());
  pder2.resize(algorithm.size());
  for(int i=0; i<algorithm.size(); ++i){

    // Locate the element
    SXNode* bnode = bnodes[i];

    // Save the node index
    algorithm[i].ind = place[bnode->temp];
    
    // Save the indices of the children
    algorithm[i].ch[0] = place[bnode->dep(0).get()->temp];
    algorithm[i].ch[1] = place[bnode->dep(1).get()->temp];

    // Operation
    algorithm[i].op = bnode->getOp();
  }

  // Indices corresponding to the inputs
  input_.resize(inputv.size());
  input_ind.resize(inputv.size());
  for(int i=0; i<inputv.size(); ++i){
    // References
    const Matrix<SX>& ip = inputv[i];
    input(i) = Matrix<double>(ip.sparsity());

    // Allocate space for the indices
    vector<int>& ii = input_ind[i];
    ii.resize(ip.size());

    // save the indices
    for(int j=0; j<ip.size(); ++j){
      ii[j] = place[ip[j].get()->temp];
    }
  }

  // Indices corresponding to each non-zero outputs
  output_.resize(outputv.size());
  output_ind.resize(output_.size());
  for(int i=0; i<outputv.size(); ++i){
    // References
    const Matrix<SX>& op = outputv[i];
    output(i) = Matrix<double>(op.size1(),op.size2(),op.col(),op.rowind());
    
    // Allocate space for the indices
    vector<int>& oi = output_ind[i];  
    oi.resize(op.size());

    // save the indices
    for(int j=0; j<op.size(); ++j){
      oi[j] = place[op[j].get()->temp];
    }
  }

  // Allocate a work vector
  work.resize(worksize,numeric_limits<double>::quiet_NaN());

  // Save the constants to the work vector
  for(vector<SXNode*>::iterator it=cnodes.begin(); it!=cnodes.end(); ++it)
    work[(**it).temp] = (**it).getValue();

  // Reset the temporary variables
  for(int i=0; i<nodes.size(); ++i){
    nodes[i]->temp = 0;
  }

}

SXFunctionInternal::~SXFunctionInternal(){
}


void SXFunctionInternal::evaluate(int nfdir, int nadir){
  
  // Copy the function arguments to the work vector
  for(int ind=0; ind<getNumInputs(); ++ind){
    const Matrix<double> &arg = input(ind);
    for(int i=0; i<arg.size(); ++i){
      work[input_ind[ind][i]] = arg[i];
    }
  }
  
  // Evaluate the algorithm
  if(nfdir==0 && nadir==0){
    // without taping
    for(vector<AlgEl>::iterator it=algorithm.begin(); it<algorithm.end(); ++it){
      // Get the arguments
      double x = work[it->ch[0]];
      double y = work[it->ch[1]];
      casadi_math<double>::fun[it->op](x,y,work[it->ind]);
    }
  } else {
    // with taping
    vector<AlgEl>::iterator it = algorithm.begin();
    vector<AlgElData<1> >::iterator it1 = pder1.begin();
    for(; it<algorithm.end(); ++it, ++it1){
      // Get the arguments
      double x = work[it->ch[0]];
      double y = work[it->ch[1]];
      casadi_math<double>::fun[it->op](x,y,work[it->ind]);
      casadi_math<double>::der[it->op](x,y,work[it->ind],it1->d);
    }
  }
  
  // Get the results
  for(int ind=0; ind<getNumOutputs(); ++ind){
    Matrix<double> &res = output(ind);
    for(int i=0; i<res.size(); ++i){
      res[i] = work[output_ind[ind][i]];
    }
  }

  // Loop over all forward directions
  for(int dir=0; dir<nfdir; ++dir){

    // Clear the seeds (not necessary if constants and parameters have zero value!)
    fill(dwork.begin(),dwork.end(),0);
    
    // Copy the function arguments to the work vector
    for(int ind=0; ind<input_.size(); ++ind){
      const Matrix<double> &fseed = fwdSeed(ind,dir);
      for(int i=0; i<fseed.size(); ++i){
        dwork[input_ind[ind][i]] = fseed[i];
      }
    }
  
    // Evaluate the algorithm for the sensitivities
    vector<AlgEl>::const_iterator it = algorithm.begin();
    vector<AlgElData<1> >::const_iterator it2 = pder1.begin();
    for(; it!=algorithm.end(); ++it, ++it2){
      dwork[it->ind] = it2->d[0] * dwork[it->ch[0]] + it2->d[1] * dwork[it->ch[1]];
    }
  
    // Get the forward sensitivities
    for(int ind=0; ind<output_.size(); ++ind){
      Matrix<double> &fsens = fwdSens(ind,dir);
      for(int i=0; i<output_ind[ind].size(); ++i){
        fsens[i] = dwork[output_ind[ind][i]];
      }
    }
  }
  
  for(int dir=0; dir<nadir; ++dir){

  // Clear the seeds (should not be necessary)
  fill(dwork.begin(),dwork.end(),0);

    // Pass the output seeds
    for(int ind=0; ind<output_.size(); ++ind){
      const Matrix<double> &aseed = adjSeed(ind,dir);
      for(int i=0; i<output_ind[ind].size(); ++i){
        dwork[output_ind[ind][i]] += aseed[i];
      }
    }

  for(int i=algorithm.size()-1; i>=0; --i){
    const AlgEl& ae = algorithm[i];
    const AlgElData<1>& aed = pder1[i];
    
    // copy the seed and clear the cache entry
    double seed = dwork[ae.ind];
    dwork[ae.ind] = 0;

    dwork[ae.ch[0]] += aed.d[0] * seed;
    dwork[ae.ch[1]] += aed.d[1] * seed;
  }

  // Collect the adjoint sensitivities
  for(int ind=0; ind<getNumInputs(); ++ind){
    Matrix<double> &asens = adjSens(ind,dir);
    for(int i=0; i<input_ind[ind].size(); ++i){
      asens[i] = dwork[input_ind[ind][i]];
    }
  }

  // Should clean up any zero added to constants and parameters!

  }
}

Matrix<SX> SXFunctionInternal::hess(int iind, int oind){
  if(output(oind).numel() != 1)
    throw CasadiException("SXFunctionInternal::hess: function must be scalar");
  
  // Reverse mode to calculate gradient
  Matrix<SX> g = grad(iind,oind);
  
  // Create function
  SXFunction gfcn(inputv.at(iind),g);
  
  // Initialize
  gfcn.init();
  
  // Return jacobian of the gradient
  return gfcn.jac();
}

Matrix<SX> SXFunctionInternal::grad(int iind, int oind){
  return trans(jac(iind,oind));
}

Matrix<SX> SXFunctionInternal::jac(int iind, int oind){
  if(input_ind.at(iind).empty() || output_ind.at(oind).empty()) return Matrix<SX>(); // quick return
  assert(input(iind).size2()==1);
  assert(output(oind).size2()==1);

  // Calculate the partial derivatives     // The loop can be executed in parallel!
  vector<SX> der1, der2;
  der1.reserve(algorithm.size());
  der2.reserve(algorithm.size());
  SX tmp[2];
  for(vector<AlgEl>::const_iterator it = algorithm.begin(); it!=algorithm.end(); ++it){
      SX f = SX(nodes[it->ind]);
      SX ch[2];
      ch[0] = SX(nodes[it->ch[0]]);
      ch[1] = SX(nodes[it->ch[1]]);
      casadi_math<SX>::der[it->op](ch[0],ch[1],f,tmp);
      if(!ch[0]->isConstant())  der1.push_back(tmp[0]);
      else                    	der1.push_back(0);

      if(!ch[1]->isConstant())  der2.push_back(tmp[1]);
      else                    	der2.push_back(0);
  }

  // Gradient (this is also the working array)
  vector<SX> g(nodes.size(),casadi_limits<SX>::zero);

  // if reverse AD
//  if(getOption("ad_mode") == "reverse"){
  if(1){ // problem with the forward mode!
    
  // Jacobian
  Matrix<SX> ret(output(oind).numel(),input_ind.at(iind).size()); 
  ret.reserve(input_ind.at(iind).size()+output(oind).numel());

#if 0
  // Backward seed (symbolic direction)
  Matrix<SX> bdir("bdir",output_.at(oind).col.size());
  for(int i=0; i<bdir.size(); ++i)
    bdir[i]->temp2 = i+1;
  
  // Pass to work array
  for(int i=0; i<output_.at(oind).size1(); ++i) // loop over rows
    for(int el=output_.at(oind).rowind_[i]; el<output_.at(oind).rowind_[i+1]; ++el){ // loop over the non-zero elements
      // add a seed to the element corresponding to the output of the function
      g[output_ind.at(oind)[el]] = bdir[el];
    }

   // backward sweep
   for(int k = algorithm.size()-1; k>=0; --k){
     const AlgEl& ae = algorithm[k];
     if(!nodes[ae.ch[0]]->isConstant()) g[ae.ch[0]] += der1[k] * g[ae.ind];
     if(!nodes[ae.ch[1]]->isConstant()) g[ae.ch[1]] += der2[k] * g[ae.ind];
     
     // Mark the place in the algorithm
     nodes[ae.ind]->temp = k;
   }

   // A row of the Jacobian
   Matrix<SX> jac_row(1,input_ind.at(iind).size());
    for(int j=0; j<input_ind[iind].size(); j++)
      if(input_ind.at(iind)[j]>=0 && !g[input_ind.at(iind)[j]]->isZero())
        jac_row[j] = g[input_ind.at(iind)[j]];

   // Loop over rows
   vector<SXNode*> deps;
   for(int i=0; i<jac_row.size(); ++i){
      // Get dependent nodes
      stack<SXNode*> s;
      s.push(jac_row[i].get());
      deps.clear();
      sort_depth_first(s,deps);
      
      // Add the dependent outputs to a new stack
      stack<SXNode*> s2;
      for(vector<SXNode*>::const_iterator it=deps.begin(); it!=deps.end(); ++it){
        for(int c=0; c<2; ++c){
          int ind = (*it)->child[c]->temp2-1;
          s2.push(outputv[oind].comp[ind].get());
        }
      }
      
      // Get the dependent nodes
/*      deps.clear();
      sort_depth_first(s2,deps);*/
 
      // Print
/*      for(vector<SXNode*>::const_iterator it=deps.begin(); it!=deps.end(); ++it)
        cout << SX(*it) << ",";
      cout << endl;*/
     
          
          /*      while(!s2.empty()){
        cout << SX(s2.top()) << ", " << endl;
        s2.pop();
      }
      cout << endl;*/
      
      

/*      cout << i << ": " << depind << endl;*/
   }
      
      
      
      
   // Reset the marks
   for(vector<AlgEl>::const_iterator it=algorithm.begin(); it!=algorithm.end(); ++it)
     nodes[it->ind]->temp = 0;
  
      cout << "ok, 1" << endl;
#endif 
#if 0
      
      
      // Mark the variables corresponding to directions
  for(int i=0; i<bdir.size(); ++i)
    bdir[i]->temp = 1;

  // Loop over rows
   vector<SXNode*> deps;
   vector<Matrix<SX> > bdir_i(jac_row.size());
   for(int i=0; i<jac_row.size(); ++i){
      // Get dependent nodes
      stack<SXNode*> s;
      s.push(jac_row[i].get());
      deps.clear();
      sort_depth_first(s,deps);
      
      for(vector<SXNode*>::iterator ii = deps.begin(); ii!=deps.end(); ++ii){
        for(int c=0; c<2; ++c)
          if((*ii)->child[c]->temp>0){
            bdir_i[i] << (*ii)->child[c];
            (*ii)->child[c]->temp = -1;
          }
      }
      
      // reset mark
      for(int j=0; j<bdir_i[i].size(); ++j)
        bdir_i[i][j]->temp = 1;
   }
   
  // Remove marks
  for(int i=0; i<bdir.size(); ++i)
    bdir[i]->temp = 0;
   
  // Create functions
  vector<SXFunction> jac_row_fcn(jac_row.size());
  for(int i=0; i<jac_row.size(); ++i){
     jac_row_fcn[i] = SXFunction(bdir_i[i],jac_row[i]);
  }
  
   // Mark where an output goes
  for(int i=0; i<bdir.size(); ++i)
    bdir[i]->temp = i;
     
  // loop over rows
  for(int i=0; i<output_.at(oind).size1(); ++i){
    
    // Loop over non-zeros in the jacobian
    for(int j=0; j<bdir_i[i].size(); ++j){
      
      // input vector
      vector<Matrix<SX> > inp(1);
      inp[0] = Matrix<SX>(bdir_i[i].size());
      inp[0][j] = SX::one;
      
      // output vector
      vector<Matrix<SX> > outp;
      
      // Evaluate
      jac_row_fcn[i]->evaluate(inp,outp);
      
      // Save the result
      if(!outp[0][0]->isZero())
        ret(i,bdir_i[i][j]->temp) = outp[0][0];
    }
  }
      
   // Remove marks
  for(int i=0; i<bdir.size(); ++i)
    
    bdir[i]->temp = 0;

  /*      cout << "i = " << i << endl;
      cout << bdir_i << endl;*/
      
/*      assert(bdir_i.size() < 10);*/
 
cout << "ok1" << endl;
return ret;
#endif 

    // Get symbolic nodes
    vector<int> snodes;
    for(int i=0; i<nodes.size(); ++i){
      if(nodes[i]->isSymbolic())
        snodes.push_back(i);
    }
              
    for(int i=0; i<output(oind).size1(); ++i) // loop over rows of the output
      for(int el=output(oind).rowind(i); el<output(oind).rowind(i+1); ++el){ // loop over the non-zero elements
        assert(output(oind).col(el) == 0); // column

        // Clear seeds (from symbolic components)
        for(vector<int>::const_iterator ii=snodes.begin(); ii!=snodes.end(); ++ii)
          g[*ii] = casadi_limits<SX>::zero;
                
        // add a seed to the element corresponding to the output of the function
        g[output_ind.at(oind)[el]] = casadi_limits<SX>::one;
        
        // backward sweep
        for(int k = algorithm.size()-1; k>=0; --k){
          const AlgEl& ae = algorithm[k];
        
          // Skip if the seed is zero
          if(!g[ae.ind]->isZero()){
        
            if(!der1[k]->isZero())
              g[ae.ch[0]] += der1[k] * g[ae.ind];
          
            if(!der2[k]->isZero())
              g[ae.ch[1]] += der2[k] * g[ae.ind];
          
            // Remove the seed
            g[ae.ind] = casadi_limits<SX>::zero;
          }
        }

        // extract and save the result
        for(int j=0; j<input_ind[iind].size(); j++){
          int ind = input_ind.at(iind)[j];
          if(ind>=0 && !g[ind]->isZero())
            ret(i,j) = g[ind];
        }
    }
    
    return ret;
  } else if(getOption("ad_mode") == "forward"){
    // Gradient
    Matrix<SX> ret(input_ind.at(iind).size(),output(oind).numel());
    ret.reserve(input_ind.at(iind).size()+output(oind).numel());
    
    for(int i=0; i<input(iind).size1(); ++i) // loop over rows of the gradient
      for(int el=input(iind).rowind(i); el<input(iind).rowind(i+1); ++el){ // loop over the non-zero elements
        assert(input(iind).col(el) == 0); // column
     
        // set all components to zero (a bit quicker than to use fill)
        for(vector<SX>::iterator it=g.begin(); it!=g.end(); ++it)
          if(!(*it)->isZero())
            *it = casadi_limits<SX>::zero;
                
        // add a seed to the element corresponding to the output of the function
        g[input_ind.at(iind)[el]] = casadi_limits<SX>::one;

        // forward sweep
        for(int k = 0; k<algorithm.size(); ++k){
          const AlgEl& ae = algorithm[k];

          // First argument
          if(!g[ae.ch[0]]->isZero() &&  !der1[k]->isZero()){
            g[ae.ind] += der1[k] * g[ae.ch[0]];
          }
          
          // Second argument
          if(!g[ae.ch[1]]->isZero() &&  !der2[k]->isZero()){
            g[ae.ind] += der2[k] * g[ae.ch[1]];
          }
          
          k++;
        }

        // extract and save the result
        for(int j=0; j<output_ind[oind].size(); j++){
          int ind = output_ind.at(oind)[j];
          if(ind>=0 && !g[ind]->isZero()){
            ret(i,j) = g[ind];
          }
        }
      }    
    return trans(ret);
    
  } else {
    throw CasadiException("SXFunctionInternal::jac: Unknown ad_mode");
  }
}

bool SXFunctionInternal::isSmooth() const{
    // Go through all nodes and check if any node is non-smooth
    for(vector<AlgEl>::const_iterator it = algorithm.begin(); it!=algorithm.end(); ++it){
      if(it->op == STEP || it->op == FLOOR )
        return false;
    }
    return true;
}

void SXFunctionInternal::print(ostream &stream) const{
 for(vector<AlgEl>::const_iterator it = algorithm.begin(); it!=algorithm.end(); ++it){
    int op = it->op;
    stringstream s,s0,s1;
    s << "i_" << it->ind;

    int i0 = it->ch[0], i1 = it->ch[1];
#if 0
    if(nodes[i0]->hasDep())  s0 << "i_" << i0;
    else                      s0 << SX(nodes[i0]);
    if(nodes[i1]->hasDep())  s1 << "i_" << i1;
    else                      s1 << SX(nodes[i1]);
#else
    s0 << "i_" << i0;
    s1 << "i_" << i1;
#endif

    stream << s.str() << " = ";
    casadi_math<double>::print[op](stream,s0.str(),s1.str());
    stream << ";" << endl;
  }
}


void SXFunctionInternal::generateCode(const string& src_name) const{
   // Output
  cout << "Generating: "<< src_name << endl;
  // Create the c source file
  ofstream cfile;
  cfile.open (src_name.c_str());

  // Print header
  cfile << "/* Automatically generated function. */" << endl;
//   cfile << "#include <stdio.h>" << endl << endl;
  cfile << "#include <math.h>" << endl << endl;

  // Dimensions
  cfile << "int n_in = " << input_.size() << ";" << endl;
  cfile << "int n_out = " << output_.size() << ";" << endl;

  // Input sizes
  // rows
  cfile << "int in_nrow[] = {";
  if(!input_.empty()){
    cfile << input(0).size1();
    for(int i=1; i<input_.size(); ++i)
      cfile << "," << input(i).size1();
  }
  cfile << "};" << endl;
  // columns
  cfile << "int in_ncol_[] = {";
  if(!input_.empty()){
    cfile << input(0).size2();
    for(int i=1; i<input_.size(); ++i)
      cfile << "," << input(i).size2();
  }
  cfile << "};" << endl;

  // Output sizes
  // rows
  cfile << "int out_nrow[] = {";
  if(!output_.empty()){
    cfile << output(0).size1();
    for(int i=1; i<output_.size(); ++i)
      cfile << "," << output(i).size1();
  }
  cfile << "};" << endl;
  // columns
  cfile << "int out_ncol_[] = {";
  if(!output_.empty()){
    cfile << output(0).size2();
    for(int i=1; i<output_.size(); ++i)
      cfile << "," << output(i).size2();
  }
  cfile << "};" << endl;

  // Memory
/*  cfile << "double i0[" << work.size() << "];" << endl;
  cfile << "double i1[" << dwork.size() << "];" << endl;;*/
  
  // Initializer
  cfile << "int init(int *n_in_, int *n_out_){" << endl;
  cfile << "*n_in_ = n_in, *n_out_ = n_out;" << endl;
  cfile << "return 0;" << endl;
  cfile << "}" << endl << endl;

  // Input sizes
  cfile << "int getInputSize(int n_in, int *n_row, int *n_col){" << endl;
  cfile << "*n_row = in_nrow[n_in]; *n_col = in_ncol_[n_in];" << endl;
  cfile << "return 0;" << endl;
  cfile << "}" << endl << endl;

  // Output sizes
  cfile << "int getOutputSize(int n_out, int *n_row, int *n_col){" << endl;
  cfile << "*n_row = out_nrow[n_out]; *n_col = out_ncol_[n_out];" << endl;
  cfile << "return 0;" << endl;
  cfile << "}" << endl << endl;

  // Evaluate function
  cfile << "int evaluate(const double** x, double** r){" << endl;

  // Which variables have been declared
  vector<bool> declared(dwork.size(),false);
  
  // Copy the function arguments to the work vector
  for(int ind=0; ind<input_.size(); ++ind){
    for(int i=0; i<input_ind[ind].size(); ++i){
      int el = input_ind[ind][i];
      cfile << "double i_" << el << "=x[" << ind << "][" << i << "];" << endl;
      declared[el] = true;
    }
  }

 // Run the algorithm
 for(vector<AlgEl>::const_iterator it = algorithm.begin(); it!=algorithm.end(); ++it){
    stringstream s0,s1;
    int op = it->op;
    if(!declared[it->ind]){
      cfile << "double ";
      declared[it->ind]=true;
    }
    cfile << "i_" << it->ind << "=";
    if(nodes[it->ch[0]]->isConstant())  s0 << nodes[it->ch[0]]->getValue();
    else                               s0 << "i_" << it->ch[0];
    if(nodes[it->ch[1]]->isConstant())  s1 << nodes[it->ch[1]]->getValue();
    else                               s1 << "i_" << it->ch[1];
    casadi_math<double>::print[op](cfile ,s0.str(),s1.str());
    cfile  << ";" << endl;
  }
  
  // Get the results
  for(int ind=0; ind<output_.size(); ++ind){
    for(int i=0; i<output_ind[ind].size(); ++i){
      cfile << "r[" << ind << "][" << i << "]=";
      int el = output_ind[ind][i];
      if(nodes[el]->isConstant())
        cfile << nodes[el]->getValue() << ";" << endl;
      else
        cfile << "i_" << el << ";" << endl;
    }
  }

  cfile << "return 0;" << endl;
  cfile << "}" << endl << endl;
  
  // Close the results file
  cfile.close();
}

void SXFunctionInternal::init(){
  // Call the init function of the base class
  XFunctionInternal::init();

  // allocate a vector with the values at the nodes, the first vector also contains the partial derivatives
  
  if(nfdir_>0 || nadir_>0)
    dwork.resize(worksize,numeric_limits<double>::quiet_NaN());
}

FX SXFunctionInternal::jacobian(int iind, int oind){
  Matrix<SX> J = jac(iind,oind); // NOTE: Multiple input, multiple output
  return SXFunction(inputv,J);
}

FX SXFunctionInternal::hessian(int iind, int oind){
  Matrix<SX> H = hess(iind,oind); // NOTE: Multiple input, multiple output
  return SXFunction(inputv,H);
}

void SXFunctionInternal::evaluateSX(const vector<Matrix<SX> >& input_s, vector<Matrix<SX> >& output_s, bool eliminate_constants){
  // Assert input dimension
  assert(input_.size() == input_s.size());
  
  // Create a symbolic work vector if not existing
  if(swork.size() != nodes.size()){
    swork.resize(nodes.size());
    for(int i=0; i<nodes.size(); ++i){
      if(!nodes[i]->hasDep())
        swork[i] = SX(nodes[i]);
    }
  }
  
  // Copy the function arguments to the work vector
  for(int ind=0; ind<input_s.size(); ++ind){
    for(int i=0; i<input_ind[ind].size(); ++i){
      swork[input_ind[ind][i]] = input_s[ind][i];
    }
  }
  
  // Evaluate the algorithm
  for(vector<AlgEl>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it){
      // Get the arguments
    SX x = swork[it->ch[0]];
    SX y = swork[it->ch[1]];
    if(eliminate_constants && x.isConstant() && y.isConstant()){
      // Check if both arguments are constants
      double temp;
      casadi_math<double>::fun[it->op](x.getValue(),y.getValue(),temp);
      swork[it->ind] = temp;
    } else {
      casadi_math<SX>::fun[it->op](x,y,swork[it->ind]);
    }
  }

  // Get the results
  for(int ind=0; ind<output_.size(); ++ind){
    for(int i=0; i<output_ind[ind].size(); ++i)
      output_s[ind][i] = swork[output_ind[ind][i]];
  }
}

SXFunctionInternal* SXFunctionInternal::clone() const{
  return new SXFunctionInternal(*this);
}


void SXFunctionInternal::clearSymbolic(){
  inputv.clear();
  outputv.clear();
  swork.clear();
}

FX SXFunctionInternal::jacobian(const std::vector<std::pair<int,int> >& jblocks){
  // Jacobian blocks
  vector<SXMatrix> jac_out(jblocks.size());
  jac_out.reserve(jblocks.size());
  for(int el=0; el<jac_out.size(); ++el){
    if(jblocks[el].second==-1){
      // Undifferentiated function
      jac_out[el] = outputv[jblocks[el].first];
    } else {
      // Jacobian
      jac_out[el] = jac(jblocks[el].second,jblocks[el].first);
    }
  }
  
  // Return function
  return SXFunction(inputv,jac_out);
}


CRSSparsity SXFunctionInternal::getJacSparsity(int iind, int oind){
  // FIXME: Inefficient algorithm
  return jac(iind,oind).sparsity();
}

} // namespace CasADi

