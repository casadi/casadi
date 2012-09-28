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

#include "ipopt_internal.hpp"
#include "ipopt_nlp.hpp"
#include "casadi/stl_vector_tools.hpp"
#include <ctime>

using namespace std;
#include <IpIpoptApplication.hpp>

// Headers for sIPOPT
#ifdef WITH_SIPOPT
#include <SensApplication.hpp>
#include <IpPDSearchDirCalc.hpp>
#include <IpIpoptAlg.hpp>
#include <SensRegOp.hpp>
#endif // WITH_SIPOPT

namespace CasADi{

IpoptInternal::IpoptInternal(const FX& F, const FX& G, const FX& H, const FX& J, const FX& GF) : NLPSolverInternal(F,G,H,J), GF_(GF){
  addOption("pass_nonlinear_variables", OT_BOOLEAN, true);
  addOption("print_time",               OT_BOOLEAN, true, "print information about execution time");
  
  // Monitors
  addOption("monitor",                  OT_STRINGVECTOR, GenericType(),  "", "eval_f|eval_g|eval_jac_g|eval_grad_f", true);

  // For passing metadata to IPOPT
  addOption("var_string_md",            OT_DICTIONARY, GenericType(), "String metadata (a dictionary with lists of strings) about variables to be passed to IPOPT");
  addOption("var_integer_md",           OT_DICTIONARY, GenericType(), "Integer metadata (a dictionary with lists of integers) about variables to be passed to IPOPT");
  addOption("var_numeric_md",           OT_DICTIONARY, GenericType(), "Numeric metadata (a dictionary with lists of reals) about variables to be passed to IPOPT");
  addOption("con_string_md",            OT_DICTIONARY, GenericType(), "String metadata (a dictionary with lists of strings) about constraints to be passed to IPOPT");
  addOption("con_integer_md",           OT_DICTIONARY, GenericType(), "Integer metadata (a dictionary with lists of integers) about constraints to be passed to IPOPT");
  addOption("con_numeric_md",           OT_DICTIONARY, GenericType(), "Numeric metadata (a dictionary with lists of reals) about constraints to be passed to IPOPT");
  
  // Sensitivity calculations with sIPOPT
  addOption("run_sens",                 OT_BOOLEAN,    false,         "Run sIPOPT sensitivity calculations");
  addOption("compute_red_hessian",      OT_BOOLEAN,    false,         "Compute the (inverse of) the reduced Hessian");
  
  // Set pointers to zero
  app_ = 0;
  userclass_ = 0;
  #ifdef WITH_SIPOPT
  app_sens_ = 0;
  #endif

  // Start an application (temporarily)
  Ipopt::IpoptApplication temp_app;

  // Get all options available in IPOPT
  map<string, Ipopt::SmartPtr<Ipopt::RegisteredOption> > regops = temp_app.RegOptions()->RegisteredOptionsList();
  for(map<string, Ipopt::SmartPtr<Ipopt::RegisteredOption> >::const_iterator it=regops.begin(); it!=regops.end(); ++it){
    // Option identifier
    string opt_name = it->first;
    
    // Short description goes here, even though we do have a longer description
    string opt_desc = it->second->ShortDescription() + " (see IPOPT documentation)";
    
    // Get the type
    Ipopt::RegisteredOptionType ipopt_type = it->second->Type();
    opt_type casadi_type;
    
    // Map Ipopt option category to a CasADi options type
    switch(ipopt_type){
      case Ipopt::OT_Number:    casadi_type = OT_REAL;          break;
      case Ipopt::OT_Integer:   casadi_type = OT_INTEGER;       break;
      case Ipopt::OT_String:    casadi_type = OT_STRING;        break;
      case Ipopt::OT_Unknown:   continue; // NOTE: No mechanism to handle OT_Unknown options
      default:                  continue; // NOTE: Unknown Ipopt options category
    }
    
    addOption(opt_name, casadi_type, GenericType(), opt_desc);
    
    // Set default values of IPOPT options 
    if (casadi_type == OT_REAL) {
      setDefault(opt_name,it->second->DefaultNumber());
    } else if (casadi_type == OT_INTEGER) {
      setDefault(opt_name,it->second->DefaultInteger());
    } else if (casadi_type == OT_STRING) {
      setDefault(opt_name,it->second->DefaultString());
    };
    
    // Save to map containing IPOPT specific options
    ops_[opt_name] = casadi_type;
  }
}

void IpoptInternal::freeIpopt(){
  // Free Ipopt application instance
  if(app_){
    delete app_; 
    app_ = 0;
  }

  // Free Ipopt user class
  if(userclass_ != 0){
    Ipopt::SmartPtr<Ipopt::TNLP> *ucptr = static_cast<Ipopt::SmartPtr<Ipopt::TNLP>*>(userclass_);
    delete ucptr;
    userclass_ = 0;
  }
  
  // Free sensitivity application
  #ifdef WITH_SIPOPT
  if(app_sens_){
    delete app_sens_;
    app_sens_ = 0;
  }
  #endif  
}

IpoptInternal::~IpoptInternal(){
  freeIpopt();
}

void IpoptInternal::init(){
  // Free existing IPOPT instance
  freeIpopt();

  // Call the init method of the base class
  NLPSolverInternal::init();
    
  // Read user options
  run_sens_ = getOption("run_sens");
  compute_red_hessian_ = getOption("compute_red_hessian");
  #ifndef WITH_SIPOPT
  casadi_assert_message(!run_sens_ && !compute_red_hessian_,  "Sensitivity calculations with sIPOPT requires that CasADi has been compiled with sIPOPT support. See CasADi documentation.");
  #endif
  
  // Gradient of the objective function, remove?
  if(!GF_.isNull()) GF_.init();
  if(!GF_.isNull()) {
    if (parametric_) {
      casadi_assert_message(GF_.getNumInputs()==2, "Wrong number of input arguments to GF for parametric NLP. Must be 2, but got " << GF_.getNumInputs());
    } else {
      casadi_assert_message(GF_.getNumInputs()==1, "Wrong number of input arguments to GF for non-parametric NLP. Must be 1, but got " << GF_.getNumInputs() << " instead. Do you perhaps intend to use fixed parameters? Then use the 'parametric' option.");
    }
    casadi_assert_message(GF_.getNumOutputs()>=1, "Wrong number of output arguments to GF");
    casadi_assert_message(GF_.input().numel()==n_,"Inconsistent dimensions");
    casadi_assert_message((GF_.output().size1()==n_ && GF_.output().size2()==1) || (GF_.output().size1()==1 && GF_.output().size2()==n_),"Inconsistent dimensions");
  }

  // Start a new IPOPT application
  app_ = new Ipopt::IpoptApplication();
    
  // Create an Ipopt user class -- need to use Ipopts spart pointer class
  Ipopt::SmartPtr<Ipopt::TNLP> *ucptr = new Ipopt::SmartPtr<Ipopt::TNLP>();
  userclass_ = static_cast<void*>(ucptr);
  Ipopt::SmartPtr<Ipopt::TNLP> &uc = *ucptr;
  uc = new IpoptUserClass(this);
  
  // read options
  exact_hessian_ = !H_.isNull();
  if(hasSetOption("hessian_approximation")){
    if(getOption("hessian_approximation")=="limited-memory"){
      exact_hessian_ = false;
    } else {
      exact_hessian_ = true;
      casadi_assert_message(!H_.isNull(), "No hessian has been provided");
    }
  } else {
    if(!exact_hessian_){
      setOption("hessian_approximation","limited-memory");
    }
  }
  
  if(verbose_){
    cout << "There are " << n_ << " variables and " << m_ << " constraints." << endl;
    if(exact_hessian_) cout << "Using exact Hessian" << endl;
    else             cout << "Using limited memory Hessian approximation" << endl;
  }
 
  bool ret = true;
       
  // Pass all the options to ipopt
  for(map<string,opt_type>::const_iterator it=ops_.begin(); it!=ops_.end(); ++it)
    if(hasSetOption(it->first)){
      GenericType op = getOption(it->first);
      switch(it->second){
        case OT_REAL:
          ret &= app_->Options()->SetNumericValue(it->first,op.toDouble(),false);
          break;
        case OT_INTEGER:
          ret &= app_->Options()->SetIntegerValue(it->first,op.toInt(),false);
          break;
        case OT_STRING:
          ret &= app_->Options()->SetStringValue(it->first,op.toString(),false);
          break;
        default:
          throw CasadiException("Illegal type");
      }
    }
  
  if (!ret) casadi_error("IpoptInternal::Init: Invalid options were detected by Ipopt.");
  // Intialize the IpoptApplication and process the options
  Ipopt::ApplicationReturnStatus status = app_->Initialize();
  if (status != Solve_Succeeded) {
    throw "Error during initialization!\n";
  }
}

void IpoptInternal::evaluate(int nfdir, int nadir){
  casadi_assert(nfdir==0 && nadir==0);

  checkInitialBounds();

  // Set the static parameter
  if (parametric_) {
    if (!F_.isNull()) F_.setInput(input(NLP_P),F_.getNumInputs()-1);
    if (!G_.isNull()) G_.setInput(input(NLP_P),G_.getNumInputs()-1);
    if (!H_.isNull()) H_.setInput(input(NLP_P),H_.getNumInputs()-1);
    if (!J_.isNull()) J_.setInput(input(NLP_P),J_.getNumInputs()-1);
    if (!GF_.isNull()) GF_.setInput(input(NLP_P),GF_.getNumInputs()-1);
  }

  // Reset the counters
  t_eval_f_ = t_eval_grad_f_ = t_eval_g_ = t_eval_jac_g_ = t_eval_h_ = t_callback_fun_ = t_callback_prepare_ = 0;
  
  // Get back the smart pointer
  Ipopt::SmartPtr<Ipopt::TNLP> *ucptr = (Ipopt::SmartPtr<Ipopt::TNLP>*)userclass_;
  Ipopt::SmartPtr<Ipopt::TNLP> &uc = *ucptr;

  // Ask Ipopt to solve the problem
  Ipopt::ApplicationReturnStatus status = app_->OptimizeTNLP(uc);
  
  if (hasOption("print_time") && bool(getOption("print_time"))) {
    // Write timings
    cout << "time spent in eval_f: " << t_eval_f_ << " s." << endl;
    cout << "time spent in eval_grad_f: " << t_eval_grad_f_ << " s." << endl;
    cout << "time spent in eval_g: " << t_eval_g_ << " s." << endl;
    cout << "time spent in eval_jac_g: " << t_eval_jac_g_ << " s." << endl;
    cout << "time spent in eval_h: " << t_eval_h_ << " s." << endl;
    cout << "time spent in callback function: " << t_callback_fun_ << " s." << endl;
    cout << "time spent in callback preparation: " << t_callback_prepare_ << " s." << endl;
  }

  if (status == Solve_Succeeded)
    stats_["return_status"] = "Solve_Succeeded";
  if (status == Solved_To_Acceptable_Level)
    stats_["return_status"] = "Solved_To_Acceptable_Level";
  if (status == Infeasible_Problem_Detected)
    stats_["return_status"] = "Infeasible_Problem_Detected";
  if (status == Search_Direction_Becomes_Too_Small)
    stats_["return_status"] = "Search_Direction_Becomes_Too_Small";
  if (status == Diverging_Iterates)
    stats_["return_status"] = "Diverging_Iterates";
  if (status == User_Requested_Stop)
    stats_["return_status"] = "User_Requested_Stop";
  if (status == Maximum_Iterations_Exceeded)
    stats_["return_status"] = "Maximum_Iterations_Exceeded";
  if (status == Restoration_Failed)
    stats_["return_status"] = "Restoration_Failed";
  if (status == Error_In_Step_Computation)
    stats_["return_status"] = "Error_In_Step_Computation";
  if (status == Not_Enough_Degrees_Of_Freedom)
    stats_["return_status"] = "Not_Enough_Degrees_Of_Freedom";
  if (status == Invalid_Problem_Definition)
    stats_["return_status"] = "Invalid_Problem_Definition";
  if (status == Invalid_Option)
    stats_["return_status"] = "Invalid_Option";
  if (status == Invalid_Number_Detected)
    stats_["return_status"] = "Invalid_Number_Detected";
  if (status == Unrecoverable_Exception)
    stats_["return_status"] = "Unrecoverable_Exception";
  if (status == NonIpopt_Exception_Thrown)
    stats_["return_status"] = "NonIpopt_Exception_Thrown";
  if (status == Insufficient_Memory)
    stats_["return_status"] = "Insufficient_Memory";
  

}

bool IpoptInternal::intermediate_callback(const double* x, const double* z_L, const double* z_U, const double* g, const double* lambda, double obj_value, int iter, double inf_pr, double inf_du,double mu,double d_norm,double regularization_size,double alpha_du,double alpha_pr,int ls_trials) {
  try {
    log("intermediate_callback started");
    double time1 = clock();
    if (!callback_.isNull()) {
#ifdef WITH_IPOPT_CALLBACK 
      copy(x,x+n_,callback_.input(NLP_X_OPT).begin());
      
      vector<double>& lambda_x = callback_.input(NLP_LAMBDA_X).data();
      for(int i=0; i<lambda_x.size(); ++i){
        lambda_x[i] = z_U[i]-z_L[i];
      }
      copy(lambda,lambda+m_,callback_.input(NLP_LAMBDA_G).begin());
      copy(g,g+m_,callback_.input(NLP_G).begin());
      
#endif // WITH_IPOPT_CALLBACK 

#ifndef WITH_IPOPT_CALLBACK 
   if (iter==0) {
      cerr << "Warning: intermediate_callback is disfunctional in your installation. You will only be able to use getStats(). See https://sourceforge.net/apps/trac/casadi/wiki/enableIpoptCallback to enable it." << endl;
   }
#endif // WITH_IPOPT_CALLBACK 
      callback_.input(NLP_COST).at(0) = obj_value;
      callback_->stats_["iter"] = iter;
      callback_->stats_["inf_pr"] = inf_pr;
      callback_->stats_["inf_du"] = inf_du;
      callback_->stats_["mu"] = mu;
      callback_->stats_["d_norm"] = d_norm;
      callback_->stats_["regularization_size"] = regularization_size;
      callback_->stats_["alpha_pr"] = alpha_pr;
      callback_->stats_["alpha_du"] = alpha_du;
      callback_->stats_["ls_trials"] = ls_trials;
      callback_.evaluate();
      double time2 = clock();
      t_callback_fun_ += double(time2-time1)/CLOCKS_PER_SEC;
      return  !callback_.output(0).at(0);
    } else {
      return 1;
    }
  } catch (exception& ex){
    if (getOption("iteration_callback_ignore_errors")) {
      cerr << "intermediate_callback: " << ex.what() << endl;
    } else {
      throw ex;
    }
    return 1;
  }
}

void IpoptInternal::finalize_solution(const double* x, const double* z_L, const double* z_U, const double* g, const double* lambda, double obj_value, int iter_count){
  try {
    // Get primal solution
    copy(x,x+n_,output(NLP_X_OPT).begin());

    // Get optimal cost
    output(NLP_COST).at(0) = obj_value;

    // Get dual solution (simple bounds)
    vector<double>& lambda_x = output(NLP_LAMBDA_X).data();
    for(int i=0; i<lambda_x.size(); ++i){
      lambda_x[i] = z_U[i]-z_L[i];
    }

    // Get dual solution (nonlinear bounds)
    copy(lambda,lambda+m_,output(NLP_LAMBDA_G).begin());
    
    // Get the constraints
    copy(g,g+m_,output(NLP_G).begin());
    
    // Get statistics
    stats_["iter_count"] = iter_count;

  } catch (exception& ex){
    cerr << "finalize_solution failed: " << ex.what() << endl;
  }
}

bool IpoptInternal::eval_h(const double* x, bool new_x, double obj_factor, const double* lambda,bool new_lambda, int nele_hess, int* iRow,int* jCol, double* values){
  try{
    log("eval_h started");
    double time1 = clock();
    if (values == NULL) {
      int nz=0;
      vector<int> rowind,col;
      H_.output().sparsity().getSparsityCRS(rowind,col);
      for(int r=0; r<rowind.size()-1; ++r)
        for(int el=rowind[r]; el<rowind[r+1]; ++el){
         if(col[el]<=r){
            iRow[nz] = r;
            jCol[nz] = col[el];
            nz++;
         }
        }
    } else {
      // Number of inputs to the hessian
      int n_hess_in = H_.getNumInputs() - (parametric_ ? 1 : 0);
      
      // Pass input
      H_.setInput(x);
      if(n_hess_in>1){
        H_.setInput(lambda, n_hess_in==4? 2 : 1);
        H_.setInput(obj_factor, n_hess_in==4? 3 : 2);
      }

      // Evaluate
      H_.evaluate();

      // Scale objective
      if(n_hess_in==1 && obj_factor!=1.0){
        for(vector<double>::iterator it=H_.output().begin(); it!=H_.output().end(); ++it){
          *it *= obj_factor;
        }
      }

      // Get results
      H_.output().get(values,SPARSESYM);
      
    }
    double time2 = clock();
    t_eval_h_ += double(time2-time1)/CLOCKS_PER_SEC;
    log("eval_h ok");
    return true;
  } catch (exception& ex){
    cerr << "eval_h failed: " << ex.what() << endl;
    return false;
  }
}

bool IpoptInternal::eval_jac_g(int n, const double* x, bool new_x,int m, int nele_jac, int* iRow, int *jCol,double* values){
  try{
    log("eval_jac_g started");
    
    // Quich finish if no constraints
    if(m==0){
      log("eval_jac_g quick return (m==0)");
      return true;
    }
    
    double time1 = clock();
    if (values == NULL) {
      int nz=0;
      vector<int> rowind,col;
      J_.output().sparsity().getSparsityCRS(rowind,col);
      for(int r=0; r<rowind.size()-1; ++r)
        for(int el=rowind[r]; el<rowind[r+1]; ++el){
  //        if(col[el]>=r){
            iRow[nz] = r;
            jCol[nz] = col[el];
            nz++;
    //      }
        }
    } else {
      // Pass the argument to the function
      J_.setInput(x);
      
       // Evaluate the function
      J_.evaluate();

      // Get the output
      J_.getOutput(values);
      
      if(monitored("eval_jac_g")){
        cout << "x = " << J_.input().data() << endl;
        cout << "J = " << endl;
        J_.output().printSparse();
      }
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

bool IpoptInternal::eval_f(int n, const double* x, bool new_x, double& obj_value)
{
  try {
    log("eval_f started");
    
    // Log time
    double time1 = clock();
    casadi_assert(n == n_);

    // Pass the argument to the function
    F_.setInput(x);
      
    // Evaluate the function
    F_.evaluate();

    // Get the result
    F_.getOutput(obj_value);

    // Printing
    if(monitored("eval_f")){
      cout << "x = " << F_.input() << endl;
      cout << "obj_value = " << obj_value << endl;
    }

    double time2 = clock();
    t_eval_f_ += double(time2-time1)/CLOCKS_PER_SEC;

    log("eval_f ok");
    return true;
  } catch (exception& ex){
    cerr << "eval_f failed: " << ex.what() << endl;
    return false;
  }
    
}

bool IpoptInternal::eval_g(int n, const double* x, bool new_x, int m, double* g)
{
  try {
    log("eval_g started");
    double time1 = clock();

    if(m>0){
      // Pass the argument to the function
      G_.setInput(x);

      // Evaluate the function and tape
      G_.evaluate();

      // Ge the result
      G_.getOutput(g);

      // Printing
      if(monitored("eval_g")){
        cout << "x = " << G_.input() << endl;
        cout << "g = " << G_.output() << endl;
      }
    }
      
    double time2 = clock();
    t_eval_g_ += double(time2-time1)/CLOCKS_PER_SEC;
    
    log("eval_g ok");
    return true;
  } catch (exception& ex){
    cerr << "eval_g failed: " << ex.what() << endl;
    return false;
  }
}

bool IpoptInternal::eval_grad_f(int n, const double* x, bool new_x, double* grad_f)
{
  try {
    log("eval_grad_f started");
    double time1 = clock();
    casadi_assert(n == n_);
    
    // If no gradient function has been provided, use AD adjoint
    if(GF_.isNull()){
    
      // Pass the argument to the function
      F_.setInput(x);
      
      // Give a seed to the function
      F_.setAdjSeed(1.0);

      // Evaluate, adjoint mode
      F_.evaluate(0,1);

      // Get the result
      F_.getAdjSens(grad_f);

      // Printing
      if(monitored("eval_grad_f")){
        cout << "grad_f = " << F_.adjSens() << endl;
      }
      
    } else {
      
      // Pass the argument to the function
      GF_.setInput(x);
      
      // Evaluate, adjoint mode
      GF_.evaluate();

      // Get the result
      GF_.getOutput(grad_f);
      
      // Printing
      if(monitored("eval_grad_f")){
        cout << "grad_f = " << GF_.output() << endl;
      }
    }
    
    double time2 = clock();
    t_eval_grad_f_ += double(time2-time1)/CLOCKS_PER_SEC;

    // Check the result for regularity
    for(int i=0; i<n; ++i){
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

bool IpoptInternal::get_bounds_info(int n, double* x_l, double* x_u,
                                int m, double* g_l, double* g_u)
{
  try {
    casadi_assert(n == n_);
    casadi_assert(m == m_);
    input(NLP_LBX).getArray(x_l,n);
    input(NLP_UBX).getArray(x_u,n);
    input(NLP_LBG).getArray(g_l,m);
    input(NLP_UBG).getArray(g_u,m);
    return true;
  } catch (exception& ex){
    cerr << "get_bounds_info failed: " << ex.what() << endl;
    return false;
  }
}

bool IpoptInternal::get_starting_point(int n, bool init_x, double* x,
                                   bool init_z, double* z_L, double* z_U,
                                   int m, bool init_lambda,
                                   double* lambda)
{
  try {
    bool warmstart = hasSetOption("warm_start_init_point") && getOption("warm_start_init_point")=="yes";
    casadi_assert_warning(init_x,"Not initializing x");
    if (warmstart) {
      casadi_assert_warning(init_lambda,"Not initializing lambda");
      casadi_assert_warning(init_z,"Not initializing z");
    }
      
    if(init_x) 
      input(NLP_X_INIT).getArray(x,n);
    
    if (init_z) {
      // Get dual solution (simple bounds)
      vector<double>& lambda_x = output(NLP_LAMBDA_X).data();
      for(int i=0; i<lambda_x.size(); ++i){
        z_L[i] = max(0.,-lambda_x[i]);
        z_U[i] = max(0., lambda_x[i]);
      }
    }
    
    if (init_lambda)
      input(NLP_LAMBDA_INIT).getArray(lambda,m);
    
    return true;
  } catch (exception& ex){
    cerr << "get_starting_point failed: " << ex.what() << endl;
    return false;
  }
}

void IpoptInternal::get_nlp_info(int& n, int& m, int& nnz_jac_g,int& nnz_h_lag)
{
  try {
    n = n_;               // number of variables
    m = m_;               // number of constraints

    // Get Jacobian sparsity pattern
    if(G_.isNull())
      nnz_jac_g = 0;
    else
      nnz_jac_g = J_.output().size();

    // Get Hessian sparsity pattern
    if(exact_hessian_)
      nnz_h_lag = H_.output().sparsity().sizeL();
    else
      nnz_h_lag = 0;
  } catch (exception& ex){
    cerr << "get_nlp_info failed: " << ex.what() << endl;
  }
}

int IpoptInternal::get_number_of_nonlinear_variables() const{
  try {
    if(H_.isNull() || !bool(getOption("pass_nonlinear_variables"))){
      // No Hessian has been interfaced
      return -1;
    } else {
      // Number of variables that appear nonlinearily
      int nv = 0;
      
      // Loop over the rows
      for(int i=0; i<H_.output().size1(); ++i){
        // If the row contains any non-zeros, the corresponding variable appears nonlinearily
        if(H_.output().rowind(i)!=H_.output().rowind(i+1))
          nv++;
      }
      
      // Return the number
      return nv;
    }
  } catch (exception& ex){
    cerr << "get_number_of_nonlinear_variables failed: " << ex.what() << endl;
    return -1;
  }
}

bool IpoptInternal::get_list_of_nonlinear_variables(int num_nonlin_vars, int* pos_nonlin_vars) const{
  try {
    // Running index
    int el = 0;
    
    // Loop over the rows
    for(int i=0; i<H_.output().size1(); ++i){
      // If the row contains any non-zeros, the corresponding variable appears nonlinearily
      if(H_.output().rowind(i)!=H_.output().rowind(i+1)){
        pos_nonlin_vars[el++] = i;
      }
    }
    
    // Assert number and return
    casadi_assert(el==num_nonlin_vars);
    return true;
  } catch (exception& ex){
    cerr << "get_list_of_nonlinear_variables failed: " << ex.what() << endl;
    return false;
  }
}

bool IpoptInternal::get_var_con_metadata(int n,
                                         map<string,vector<string> >& var_string_md, 
                                         map<string,vector<int> >& var_integer_md,
                                         map<string,vector<double> >& var_numeric_md,
                                         int m,
                                         map<string,vector<string> >& con_string_md,
                                         map<string,vector<int> >& con_integer_md,
                                         map<string,vector<double> >& con_numeric_md){
  if(hasSetOption("var_string_md")){
    Dictionary dict = getOption("var_string_md");
    for(Dictionary::const_iterator it=dict.begin(); it!=dict.end(); ++it){
      string key = it->first; // Get the key
      vector<string> entry = it->second; // Get the entry
      casadi_assert_message(entry.size()==n, "Inconsistent length of IPOPT metadata."); // Check length for consistency
      var_string_md[key] = entry; // Save to IPOPT data structure
    }
  }

  if(hasSetOption("var_integer_md")){
    Dictionary dict = getOption("var_integer_md");
    for(Dictionary::const_iterator it=dict.begin(); it!=dict.end(); ++it){
      string key = it->first; // Get the key
      vector<int> entry = it->second; // Get the entry
      casadi_assert_message(entry.size()==n, "Inconsistent length of IPOPT metadata."); // Check length for consistency
      var_integer_md[key] = entry; // Save to IPOPT data structure
    }
  }

  if(hasSetOption("var_numeric_md")){
    Dictionary dict = getOption("var_numeric_md");
    for(Dictionary::const_iterator it=dict.begin(); it!=dict.end(); ++it){
      string key = it->first; // Get the key
      vector<double> entry = it->second; // Get the entry
      casadi_assert_message(entry.size()==n, "Inconsistent length of IPOPT metadata."); // Check length for consistency
      var_numeric_md[key] = entry; // Save to IPOPT data structure
    }
  }
  
  if(hasSetOption("con_string_md")){
    Dictionary dict = getOption("con_string_md");
    for(Dictionary::const_iterator it=dict.begin(); it!=dict.end(); ++it){
      string key = it->first; // Get the key
      vector<string> entry = it->second; // Get the entry
      casadi_assert_message(entry.size()==m, "Inconsistent length of IPOPT metadata."); // Check length for consistency
      con_string_md[key] = entry; // Save to IPOPT data structure
    }
  }

  if(hasSetOption("con_integer_md")){
    Dictionary dict = getOption("con_integer_md");
    for(Dictionary::const_iterator it=dict.begin(); it!=dict.end(); ++it){
      string key = it->first; // Get the key
      vector<int> entry = it->second; // Get the entry
      casadi_assert_message(entry.size()==m, "Inconsistent length of IPOPT metadata."); // Check length for consistency
      con_integer_md[key] = entry; // Save to IPOPT data structure
    }
  }

  if(hasSetOption("con_numeric_md")){
    Dictionary dict = getOption("con_numeric_md");
    for(Dictionary::const_iterator it=dict.begin(); it!=dict.end(); ++it){
      string key = it->first; // Get the key
      vector<double> entry = it->second; // Get the entry
      casadi_assert_message(entry.size()==m, "Inconsistent length of IPOPT metadata."); // Check length for consistency
      con_numeric_md[key] = entry; // Save to IPOPT data structure
    }
  }
  
  return true;
}

void IpoptInternal::finalize_metadata(int n,
                                      const std::map<std::string,std::vector<std::string> >& var_string_md, 
                                      const std::map<std::string,std::vector<int> >& var_integer_md,
                                      const std::map<std::string,std::vector<double> >& var_numeric_md,
                                      int m,
                                      const std::map<std::string,std::vector<std::string> >& con_string_md,
                                      const std::map<std::string,std::vector<int> >& con_integer_md,
                                      const std::map<std::string,std::vector<double> >& con_numeric_md){
  
}

} // namespace CasADi
