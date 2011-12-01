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

#include "nlp_solver_internal.hpp"
#include "mx_function.hpp"
#include "sx_function.hpp"
#include "../sx/sx_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../fx/fx_tools.hpp"

INPUTSCHEME(NLPInput)
OUTPUTSCHEME(NLPOutput)

using namespace std;
namespace CasADi{

NLPSolverInternal::NLPSolverInternal(const FX& F, const FX& G, const FX& H, const FX& J) : F_(F), G_(G), H_(H), J_(J){
  // set default options
  setOption("name",            "unnamed NLP solver"); // name of the function
  addOption("expand_f",         OT_BOOLEAN,     false,         "Expand the objective function in terms of scalar operations, i.e. MX->SX");
  addOption("expand_g",         OT_BOOLEAN,     false,         "Expand the constraint function in terms of scalar operations, i.e. MX->SX");
  addOption("generate_hessian", OT_BOOLEAN,     false,         "Generate an exact Hessian of the Lagrangian");
  addOption("iteration_callback", OT_FX,     FX(),            "A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.");
  addOption("ignore_check_vec", OT_BOOLEAN,     false,            "If set to true, the input shape of F will not be checked.");
  
  n_ = 0;
  m_ = 0;
  pn_ = 0;
  pm_ = 0;
}

NLPSolverInternal::~NLPSolverInternal(){
}

void NLPSolverInternal::init(){
  // Read the verbosity option already now
  verbose_ = getOption("verbose");
  
  // Initialize the functions
  casadi_assert_message(!F_.isNull(),"No objective function");
  if(!F_.isInit()){
    F_.init();
    log("Objective function initialized");
  }
  if(!G_.isNull() && !G_.isInit()){
    G_.init();
    log("Constraint function initialized");
  }

  // Get dimensions
  n_ = F_.input(0).numel();
  m_ = G_.isNull() ? 0 : G_.output(0).numel();

  // Basic sanity checks
  casadi_assert_message(F_.getNumInputs()==1 || F_.getNumInputs()==2, "Wrong number of input arguments to F. Must be 1 or 2");
  casadi_assert_message(getOption("ignore_check_vec") || F_.input().size2()==1,
     "To avoid confusion, the input argument to F must be vector. You supplied " << F_.input().dimString() << endl <<
     " We suggest you make the following changes:" << endl <<
     "   -  F is an SXFunction:  SXFunction([X],[rhs]) -> SXFunction([vec(X)],[rhs])" << endl <<
     "             or            F -                   ->  F = vec(F) " << 
     "   -  F is an MXFunction:  MXFunction([X],[rhs]) -> " <<  endl <<
     "                                     X_vec = MX(\"X\",vec(X.sparsity())) " << endl <<
     "                                     F_vec = MXFunction([X_flat],[F.call([X_flat.reshape(X.sparsity())])[0]]) " << endl <<
     "             or            F -                   ->  F = vec(F) " << 
     " You may ignore this warning by setting the 'ignore_check_vec' option to true." << endl
  );
  
  casadi_assert_message(F_.getNumOutputs()>=1, "Wrong number of output arguments to F");
  casadi_assert_message(F_.output().scalar() && F_.output().dense(), "Output argument of F not dense scalar.");
  casadi_assert_message(F_.input().dense(), "Input argument of F must be dense. You supplied " << F_.input().dimString());
  
  if(!G_.isNull()) {
    casadi_assert_message(G_.getNumInputs()>=1, "Wrong number of input arguments to G");
    casadi_assert_message(G_.getNumOutputs()>=1, "Wrong number of output arguments to G");
    casadi_assert_message(G_.input().numel()==n_, "Inconsistent dimensions");
    casadi_assert_message(G_.input().sparsity()==F_.input().sparsity(), "F and G input dimension must match. F " << F_.input().dimString() << ". G " << G_.input().dimString());
  }
  
  // Find out if we are to expand the objective function in terms of scalar operations
  bool expand_f = getOption("expand_f");
  if(expand_f){
    log("Expanding objective function");
    
    // Cast to MXFunction
    MXFunction F_mx = shared_cast<MXFunction>(F_);
    if(F_mx.isNull()){
      casadi_warning("Cannot expand objective function as it is not an MXFunction");
    } else {
      // Take use the input scheme of G if possible (it might be an SXFunction)
      vector<SXMatrix> inputv;
      if(!G_.isNull() && F_.getNumInputs()==G_.getNumInputs()){
        inputv = G_.symbolicInputSX();
      } else {
        inputv = F_.symbolicInputSX();
      }
      
      // Try to expand the MXFunction
      F_ = F_mx.expand(inputv);
      F_.setOption("number_of_fwd_dir",F_mx.getOption("number_of_fwd_dir"));
      F_.setOption("number_of_adj_dir",F_mx.getOption("number_of_adj_dir"));
      F_.init();
    }
  }
  
  
  // Find out if we are to expand the constraint function in terms of scalar operations
  bool expand_g = getOption("expand_g");
  if(expand_g){
    log("Expanding constraint function");
    
    // Cast to MXFunction
    MXFunction G_mx = shared_cast<MXFunction>(G_);
    if(G_mx.isNull()){
      casadi_warning("Cannot expand constraint function as it is not an MXFunction");
    } else {
      // Take use the input scheme of F if possible (it might be an SXFunction)
      vector<SXMatrix> inputv;
      if(F_.getNumInputs()==G_.getNumInputs()){
        inputv = F_.symbolicInputSX();
      } else {
        inputv = G_.symbolicInputSX();
      }
      
      // Try to expand the MXFunction
      G_ = G_mx.expand(inputv);
      G_.setOption("number_of_fwd_dir",G_mx.getOption("number_of_fwd_dir"));
      G_.setOption("number_of_adj_dir",G_mx.getOption("number_of_adj_dir"));
      G_.init();
    }
  }
  
  // Find out if we are to expand the constraint function in terms of scalar operations
  bool generate_hessian = getOption("generate_hessian");
  if(generate_hessian && H_.isNull()){
    log("generating hessian");
    
    // Simple if unconstrained
    if(G_.isNull()){
      // Create Hessian of the objective function
      FX HF = F_.hessian();
      HF.init();
      
      // Symbolic inputs of HF
      vector<MX> HF_in = F_.symbolicInput();
      
      // Lagrange multipliers
      MX lam("lam",0);
      
      // Objective function scaling
      MX sigma("sigma");
      
      // Inputs of the Hessian function
      vector<MX> H_in = HF_in;
      H_in.insert(H_in.begin()+1, lam);
      H_in.insert(H_in.begin()+2, sigma);
      
      // Get an expression for the Hessian of F
      MX hf = HF.call(HF_in).at(0);
      
      // Create the scaled Hessian function
      H_ = MXFunction(H_in, sigma*hf);
      log("Unconstrained Hessian function generated");
      
    } else { // G_.isNull()
      
      // Check if the functions are SXFunctions
      SXFunction F_sx = shared_cast<SXFunction>(F_);
      SXFunction G_sx = shared_cast<SXFunction>(G_);
      
      // Efficient if both functions are SXFunction
      if(!F_sx.isNull() && !G_sx.isNull()){
        // Parametric NLP not supported
        casadi_assert_message(G_.getNumInputs()==1 && F_.getNumInputs()==1, "exact hessian generation currently not supported for parametric NLP");

        // Expression for f and g
        SXMatrix f = F_sx.outputSX();
        SXMatrix g = G_sx.outputSX();
        
        // Numeric hessian
        bool f_num_hess = F_sx.getOption("numeric_hessian");
        bool g_num_hess = G_sx.getOption("numeric_hessian");
        
        // Number of derivative directions
        bool f_num_fwd = F_sx.getOption("number_of_fwd_dir");
        bool g_num_fwd = G_sx.getOption("number_of_fwd_dir");
        bool f_num_adj = F_sx.getOption("number_of_adj_dir");
        bool g_num_adj = G_sx.getOption("number_of_adj_dir");
        
        // Substitute symbolic variables in f if different input variables from g
        if(!isEqual(F_sx.inputSX(),G_sx.inputSX())){
          f = substitute(f,F_sx.inputSX(),G_sx.inputSX());
        }
        
        // Lagrange multipliers
        SXMatrix lam = ssym("lambda",g.size1());

        // Objective function scaling
        SXMatrix sigma = ssym("sigma");        
        
        // Lagrangian function
        vector<SXMatrix> lfcn_in(3);
        lfcn_in[0] = G_sx.inputSX();
        lfcn_in[1] = lam;
        lfcn_in[2] = sigma;
        SXFunction lfcn(lfcn_in, sigma*f + inner_prod(lam,g));
        lfcn.setOption("verbose",getOption("verbose"));
        lfcn.setOption("numeric_hessian",f_num_hess || g_num_hess);
        lfcn.setOption("number_of_fwd_dir",std::min(f_num_fwd,g_num_fwd));
        lfcn.setOption("number_of_adj_dir",std::min(f_num_adj,g_num_adj));
        lfcn.init();
        
        // Hessian of the Lagrangian
        H_ = static_cast<FX&>(lfcn).hessian();
        H_.setOption("verbose",getOption("verbose"));
        log("SX Hessian function generated");
        
      } else { // !F_sx.isNull() && !G_sx.isNull()
        // Check if the functions are SXFunctions
        MXFunction F_mx = shared_cast<MXFunction>(F_);
        MXFunction G_mx = shared_cast<MXFunction>(G_);
        
        // If they are, check if the arguments are the same
        if(!F_mx.isNull() && !G_mx.isNull() && isEqual(F_mx.inputMX(),G_mx.inputMX())){
          casadi_warning("Exact Hessian calculation for MX is still experimental");
          
          // Parametric NLP not supported
          casadi_assert_message(G_mx.getNumInputs()==1 && F_mx.getNumInputs()==1, "exact hessian generation currently not supported for parametric NLP");
          
          // Expression for f and g
          MX f = F_mx.outputMX();
          MX g = G_mx.outputMX();
          
          // Lagrange multipliers
          MX lam("lam",g.size1());
      
          // Objective function scaling
          MX sigma("sigma");

          // Inputs of the Lagrangian function
          vector<MX> lfcn_in(3);
          lfcn_in[0] = G_mx.inputMX();
          lfcn_in[1] = lam;
          lfcn_in[2] = sigma;

          // Lagrangian function
          MXFunction lfcn(lfcn_in,sigma*f+ inner_prod(lam,g));
          lfcn.init();
      
          bool adjoint_mode = true;
          if(adjoint_mode){
          
            // Gradient of the lagrangian
            MX gL = trans(lfcn.grad().at(0));
            
            MXFunction glfcn(lfcn_in,gL);
            glfcn.setOption("number_of_fwd_dir",n_);
            glfcn.init();
            
            // Hessian of the Lagrangian
            H_ = glfcn.jacobian();
          } else {

            // Hessian of the Lagrangian
            H_ = lfcn.hessian();
            
          }
          log("MX Hessian function generated");
          
        } else {
          casadi_assert_message(0, "Automatic calculation of exact Hessian currently only for F and G both SXFunction or MXFunction ");
        }
      } // !F_sx.isNull() && !G_sx.isNull()
    } // G_.isNull()
  } // generate_hessian && H_.isNull()
  if(!H_.isNull() && !H_.isInit()) {
    H_.init();
    log("Hessian function initialized");
  }

  // Create a Jacobian if it does not already exists
  if(!G_.isNull() && J_.isNull()){
    J_ = G_.jacobian();
    log("Jacobian function generated");
  }
  if(!J_.isNull() && !J_.isInit()){
    J_.init();
    log("Jacobian function initialized");
  }

  
  if(!H_.isNull()) {
    casadi_assert_message(H_.getNumInputs()>=1, "Wrong number of input arguments to H");
    casadi_assert_message(H_.getNumOutputs()>=1, "Wrong number of output arguments to H");
    casadi_assert_message(H_.input(0).numel()==n_,"Inconsistent dimensions");
    casadi_assert_message(H_.output().size1()==n_,"Inconsistent dimensions");
    casadi_assert_message(H_.output().size2()==n_,"Inconsistent dimensions");
  }

  if(!J_.isNull()){
    casadi_assert_message(J_.getNumInputs()>=1, "Wrong number of input arguments to J");
    casadi_assert_message(J_.getNumOutputs()>=1, "Wrong number of output arguments to J");
    casadi_assert_message(J_.input().numel()==n_,"Inconsistent dimensions");
    casadi_assert_message(J_.output().size2()==n_,"Inconsistent dimensions");
  }
  
  pn_=0;
  pm_=0;
  
  int pn=0;
  int pm=0;
    
  // Check if any of the functions have a second argument (i.e. to pass parameters)
  FX* ff[4] = {&F_, &G_, &H_, &J_};
  for (int k=0; k<4; ++k){
    if (ff[k]->isNull()) continue;
    if (ff[k]->getNumInputs()!=2) continue;
    pn = ff[k]->input(1).size1();
    pm = ff[k]->input(1).size2();
    
    if (pn==0 || pm==0)
     continue;
    
    casadi_assert_message((pn==pn_ && pm==pm_) || pn_==0 || pm_==0,
      "One of your supplied functions had a second input argument, which was interpreted as a parameter of shape (" << pn_ << "x" << pm_ << ")." << std::endl <<
      "However, another function had a second input argument of shape (" << pn << "x" << pm << ")." << std::endl <<
      "This is inconsistent."
    );

    pn_ = pn;
    pm_ = pm;
  }
  
  // Infinity
  double inf = numeric_limits<double>::infinity();
  
  // Allocate space for inputs
  input_.resize(NLP_NUM_IN);
  input(NLP_X_INIT)      = DMatrix(n_,1,0);
  input(NLP_LBX)         = DMatrix(n_,1,-inf);
  input(NLP_UBX)         = DMatrix(n_,1, inf);
  input(NLP_LBG)         = DMatrix(m_,1,-inf);
  input(NLP_UBG)         = DMatrix(m_,1, inf);
  input(NLP_LAMBDA_INIT) = DMatrix(m_,1,0);
  input(NLP_P)           = DMatrix(pn_,pm_,0);
  
  // Allocate space for outputs
  output_.resize(NLP_NUM_OUT);
  output(NLP_X_OPT)      = DMatrix(n_,1,0);
  output(NLP_COST)       = DMatrix(1,1,0);
  output(NLP_LAMBDA_OPT) = DMatrix(m_,1,0);
  output(NLP_LAMBDA_LBX) = DMatrix(n_,1,0);
  output(NLP_LAMBDA_UBX) = DMatrix(n_,1,0);
  
  
  if (hasSetOption("iteration_callback")) {
   callback_ = getOption("iteration_callback");
   if (!callback_.isNull()) {
     if (!callback_.isInit()) callback_.init();
     casadi_assert_message(callback_.getNumOutputs()==1, "Callback function should have one output, a scalar that indicates wether to break. 0 = continue");
     casadi_assert_message(callback_.output(0).size()==1, "Callback function should have one output, a scalar that indicates wether to break. 0 = continue");
     casadi_assert_message(callback_.getNumInputs()==NLP_NUM_OUT, "Callback function should have the output scheme of NLPSolver as input scheme. i.e. " <<NLP_NUM_OUT << " inputs instead of the " << callback_.getNumInputs() << " you provided." );
     for (int i=0;i<NLP_NUM_OUT;i++) {
       casadi_assert_message(callback_.input(i).sparsity()==output(i).sparsity(),
         "Callback function should have the output scheme of NLPSolver as input scheme. " << 
         "Input #" << i << " was found to be " << callback_.input(i).dimString() << " instead of expected " << output(i).dimString() << "."
       );
       callback_.input(i).setAll(0);
     }
   }
  }

  // Call the initialization method of the base class
  FXInternal::init();
}


   

  void NLPSolverInternal::reportConstraints(std::ostream &stream) { 
  
    stream << "Reporting NLP constraints" << endl;
    CasADi::reportConstraints(stream,output(NLP_X_OPT),input(NLP_LBX),input(NLP_UBX), "decision bounds");
    
    G_.input(0).set(output(NLP_X_OPT));
    if (G_.getNumInputs()==2)  G_.input(1).set(output(NLP_P));
    G_.evaluate();
    CasADi::reportConstraints(stream,G_.output(),input(NLP_LBG),input(NLP_UBG), "constraints");
  }


} // namespace CasADi
