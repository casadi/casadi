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

#include "snopt_internal.hpp"
#include "symbolic/stl_vector_tools.hpp"
#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/mx/mx_tools.hpp"
#include "symbolic/matrix/sparsity_tools.hpp"
#include "symbolic/fx/mx_function.hpp"
#include <ctime>

#include "wsnopt.hpp"

using namespace std;

namespace CasADi{

  SnoptInternal::SnoptInternal(const FX& nlp) : NLPSolverInternal(nlp){

  }

  SnoptInternal::~SnoptInternal(){

  }
  
  SnoptInternal* SnoptInternal::clone() const{ 
    // Use default copy routine
    SnoptInternal* node = new SnoptInternal(*this); 
    
    return node;
  }

  void SnoptInternal::init(){

    // Call the init method of the base class
    NLPSolverInternal::init();

  }

  void SnoptInternal::reset(){

  }

  void SnoptInternal::setQPOptions() {
    setOption("UserHM", true);
  }

  void SnoptInternal::passOptions() {

  }

  std::string SnoptInternal::formatStatus(int status) const {
    if (status_.find(status)==status_.end()) {
      std::stringstream ss;
      ss << "Unknown status: " << status;
      return ss.str();
    } else {
      return (*status_.find(status)).second;
    }
  }

  void SnoptInternal::evaluate(){
    log("SnoptInternal::evaluate");    
    snopt();    
  }

  bool SnoptInternal::eval_h(const double* x, double obj_factor, const double* lambda, double* values){
    try{
      log("eval_h started");
      double time1 = clock();

      // Make sure generated
      casadi_assert_warning(!hessLag_.isNull(),"Hessian function not pregenerated");

      // Get Hessian 
      FX& hessLag = this->hessLag();

      // Pass input
      hessLag.setInput(x,HESSLAG_X);
      hessLag.setInput(input(NLP_SOLVER_P),HESSLAG_P);
      hessLag.setInput(obj_factor,HESSLAG_LAM_F);
      hessLag.setInput(lambda,HESSLAG_LAM_G);

      // Evaluate
      hessLag.evaluate();

      // Get results
      const DMatrix& H = hessLag.output();
      const vector<int>& rowind = H.rowind();
      const vector<int>& col = H.col();
      const vector<double>& data = H.data();
      
      if(monitored("eval_h")){
        std::cout << "x = " <<  hessLag.input(HESSLAG_X) << std::endl;
        std::cout << "obj_factor= " << obj_factor << std::endl;
        std::cout << "lambda = " << hessLag.input(HESSLAG_LAM_G) << std::endl;
        std::cout << "H = " << hessLag.output(HESSLAG_HESS) << std::endl;
      }

      if (regularity_check_ && !isRegular(hessLag.output(HESSLAG_HESS).data())) casadi_error("SnoptInternal::eval_h: NaN or Inf detected.");
      
      double time2 = clock();
      t_eval_h_ += double(time2-time1)/CLOCKS_PER_SEC;
      log("eval_h ok");
      return true;
    } catch (exception& ex){
      cerr << "eval_h failed: " << ex.what() << endl;
      return false;
    }
  }

  bool SnoptInternal::eval_jac_g(const double* x,double* values){
    try{
      log("eval_jac_g started");
   
      // Make sure generated
      casadi_assert(!jacG_.isNull());
 
      // Get Jacobian
      FX& jacG = this->jacG();

      double time1 = clock();

      // Pass the argument to the function
      jacG.setInput(x,JACG_X);
      jacG.setInput(input(NLP_SOLVER_P),JACG_P);
    
      // Evaluate the function
      jacG.evaluate();

      // Transpose the result
      const DMatrix& J = jacG.output(JACG_JAC);
      const vector<double>& J_data = J.data();
      const vector<int>& J_rowind = J.rowind();
      const vector<int>& JT_col = spJacG_T_.col();
      copy(J_rowind.begin(),J_rowind.end(),jacG_tmp_.begin());
      for(vector<int>::const_iterator i=JT_col.begin(); i!=JT_col.end(); ++i){
        *values++ = J_data[jacG_tmp_[*i]++];
      }
    
      if(monitored("eval_jac_g")){
        cout << "x = " << jacG_.input().data() << endl;
        cout << "J = " << endl;
        jacG_.output().printSparse();
      }
    
      double time2 = clock();
      t_eval_jac_g_ += double(time2-time1)/CLOCKS_PER_SEC;
    
      log("eval_jac_g ok");
      return true;
    } catch (exception& ex){
      cerr << "eval_jac_g failed: " << ex.what() << endl;
      return false;
    }
  }

  bool SnoptInternal::eval_f(const double* x, double scale, double& obj_value){
    try {
      log("eval_f started");
    
      // Log time
      double time1 = clock();

      // Pass the argument to the function
      nlp_.setInput(x, NL_X);
      nlp_.setInput(input(NLP_SOLVER_P),NL_P);
      
      // Evaluate the function
      nlp_.evaluate();

      // Get the result
      nlp_.getOutput(obj_value,NL_F);

      // Printing
      if(monitored("eval_f")){
        cout << "x = " << nlp_.input(NL_X) << endl;
        cout << "obj_value = " << obj_value << endl;
      }
      obj_value *= scale;

      if (regularity_check_ && !isRegular(nlp_.output().data())) casadi_error("SnoptInternal::eval_f: NaN or Inf detected.");

      double time2 = clock();
      t_eval_f_ += double(time2-time1)/CLOCKS_PER_SEC;

      log("eval_f ok");
      return true;
    } catch (exception& ex){
      cerr << "eval_f failed: " << ex.what() << endl;
      return false;
    }
  }

  bool SnoptInternal::eval_g(const double* x, double* g)
  {
    try {
      log("eval_g started");
      double time1 = clock();


      // Pass the argument to the function
      nlp_.setInput(x,NL_X);
      nlp_.setInput(input(NLP_SOLVER_P),NL_P);

      // Evaluate the function and tape
      nlp_.evaluate();

      // Ge the result
      nlp_.getOutput(g,NL_G);

      // Printing
      if(monitored("eval_g")){
        cout << "x = " << nlp_.input(NL_X) << endl;
        cout << "g = " << nlp_.output(NL_G) << endl;
      }
      

      if (regularity_check_ && !isRegular(nlp_.output(NL_G).data())) casadi_error("SnoptInternal::eval_g: NaN or Inf detected.");
    
      double time2 = clock();
      t_eval_g_ += double(time2-time1)/CLOCKS_PER_SEC;
    
      log("eval_g ok");
      return true;
    } catch (exception& ex){
      cerr << "eval_g failed: " << ex.what() << endl;
      return false;
    }
  }

  bool SnoptInternal::eval_grad_f(const double* x, double scale , double* grad_f )
  {
    try {
      log("eval_grad_f started");
      double time1 = clock();
    
      // Pass the argument to the function
      gradF_.setInput(x,NL_X);
      gradF_.setInput(input(NLP_SOLVER_P),NL_P);
      
      // Evaluate, adjoint mode
      gradF_.evaluate();
      
      // Get the result
      gradF_.output().get(grad_f,DENSE);

      // Scale
      for(int i=0; i<nx_; ++i){
        grad_f[i] *= scale;
      }
      
      // Printing
      if(monitored("eval_grad_f")){
        cout << "grad_f = " << gradF_.output() << endl;
      }
      
      if (regularity_check_ && !isRegular(gradF_.output().data())) casadi_error("SnoptInternal::eval_grad_f: NaN or Inf detected.");
    
      double time2 = clock();
      t_eval_grad_f_ += double(time2-time1)/CLOCKS_PER_SEC;

      // Check the result for regularity
      for(int i=0; i<nx_; ++i){
        if(isnan(grad_f[i]) || isinf(grad_f[i])){
          log("eval_grad_f: result not regular");
          return false;
        }
      }

      log("eval_grad_f ok");
      return true;
    } catch (exception& ex){
      cerr << "eval_jac_f failed: " << ex.what() << endl;
      return false;
    }
  }

  void SnoptInternal::setOptionsFromFile(const std::string & file) {
    
  }


} // namespace CasADi

