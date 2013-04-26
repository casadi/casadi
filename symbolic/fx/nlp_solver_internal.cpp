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
  addOption("generate_hessian", OT_BOOLEAN,     false,         "Generate an exact Hessian of the Lagrangian if not supplied");
  addOption("generate_jacobian", OT_BOOLEAN,     true,         "Generate an exact Jacobian of the constraints if not supplied");
  addOption("iteration_callback", OT_FX,     FX(),            "A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.");
  addOption("iteration_callback_step", OT_INTEGER,     1,       "Only call the callback function every few iterations.");
  addOption("iteration_callback_ignore_errors", OT_BOOLEAN,     false,      "If set to true, errors thrown by iteration_callback will be ignored.");
  addOption("ignore_check_vec", OT_BOOLEAN,     false,            "If set to true, the input shape of F will not be checked.");
  addOption("warn_initial_bounds", OT_BOOLEAN,     false,       "Warn if the initial guess does not satisfy LBX and UBX");
  addOption("parametric", OT_BOOLEAN, false, "Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.");
  addOption("gauss_newton",      OT_BOOLEAN,  false,           "Use Gauss Newton Hessian approximation");

  nx_ = 0;
  ng_ = 0;
  
  inputScheme_ = SCHEME_NLPInput;
  outputScheme_ = SCHEME_NLPOutput;

}

NLPSolverInternal::~NLPSolverInternal(){
}

void NLPSolverInternal::init(){
  // Read options
  verbose_ = getOption("verbose");
  gauss_newton_ = getOption("gauss_newton");
  
  // Initialize the functions
  casadi_assert_message(!F_.isNull(),"No objective function");
  if(!F_.isInit()){
    F_.init();
    log("Objective function initialized");
  }
  if(!G_.isNull() && !G_.isInit()){
    G_.setOption("verbose",getOption("verbose"));
    G_.init();
    log("Constraint function initialized");
  }

  // Get dimensions
  nx_ = F_.input(0).numel();
  ng_ = G_.isNull() ? 0 : G_.output(0).numel();

  // Remove G if it has zero dimension
  if(ng_ == 0 && !G_.isNull())
    G_ = FX();

  parametric_ = getOption("parametric");
  
  if (parametric_) {
    casadi_assert_message(F_.getNumInputs()==2, "Wrong number of input arguments to F for parametric NLP. Must be 2, but got " << F_.getNumInputs());
  } else {
    casadi_assert_message(F_.getNumInputs()==1, "Wrong number of input arguments to F for non-parametric NLP. Must be 1, but got " << F_.getNumInputs() << " instead. Do you perhaps intend to use fixed parameters? Then use the 'parametric' option.");
  }

  // Basic sanity checks
  casadi_assert_message(F_.getNumInputs()==1 || F_.getNumInputs()==2, "Wrong number of input arguments to F. Must be 1 or 2");
  
  if (F_.getNumInputs()==2) parametric_=true;
  casadi_assert_message(getOption("ignore_check_vec") || gauss_newton_ || F_.input().size2()==1,
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
  casadi_assert_message(gauss_newton_  || F_.output().scalar(), "Output argument of F not scalar.");
  casadi_assert_message(F_.output().dense(), "Output argument of F not dense.");
  casadi_assert_message(F_.input().dense(), "Input argument of F must be dense. You supplied " << F_.input().dimString());
  
  if(!G_.isNull()) {
    if (parametric_) {
      casadi_assert_message(G_.getNumInputs()==2, "Wrong number of input arguments to G for parametric NLP. Must be 2, but got " << G_.getNumInputs());
    } else {
      casadi_assert_message(G_.getNumInputs()==1, "Wrong number of input arguments to G for non-parametric NLP. Must be 1, but got " << G_.getNumInputs() << " instead. Do you perhaps intend to use fixed parameters? Then use the 'parametric' option.");
    }
    casadi_assert_message(G_.getNumOutputs()>=1, "Wrong number of output arguments to G");
    casadi_assert_message(G_.input().numel()==nx_, "Inconsistent dimensions");
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
      G_.init();
    }
  }
  
  // Find out if we are to expand the constraint function in terms of scalar operations
  bool generate_hessian = getOption("generate_hessian");
  if(generate_hessian && H_.isNull()){
    casadi_assert_message(!gauss_newton_,"Automatic generation of Gauss-Newton Hessian not yet supported");
    log("generating hessian");
    
    if(G_.isNull()){ // unconstrained
      // Calculate Hessian and wrap to get required syntax
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
      H_.setOption("name","nlp_hessian");
      log("Unconstrained Hessian function generated");
      
    } else { // Constrained
      
      // SXFunction if both functions are SXFunction
      if(is_a<SXFunction>(F_) && is_a<SXFunction>(G_)){
        SXFunction F = shared_cast<SXFunction>(F_);
        SXFunction G = shared_cast<SXFunction>(G_);
        
        // Expression for f and g
        SXMatrix g = G.outputExpr(0);
        SXMatrix f = substitute(F.outputExpr(),F.inputExpr(),G.inputExpr()).front();
        
        // Lagrange multipliers
        SXMatrix lam = ssym("lambda",g.size1());

        // Objective function scaling
        SXMatrix sigma = ssym("sigma");        
        
        // Lagrangian function
        vector<SXMatrix> lfcn_in(parametric_? 4: 3);
        lfcn_in[0] = G.inputExpr(0);
        lfcn_in[1] = lam;
        lfcn_in[2] = sigma;
        if (parametric_) lfcn_in[3] = G.inputExpr(1);
        SXFunction lfcn(lfcn_in, sigma*f + inner_prod(lam,g));
        lfcn.setOption("verbose",verbose());
        lfcn.setOption("name","nlp_lagrangian");
        lfcn.init();
        if(verbose()) 
          cout << "SX Lagrangian function generated: algorithm size " << lfcn.getAlgorithmSize() << ", work size = " << lfcn.getWorkSize() << endl;
        
        // Hessian of the Lagrangian
        H_ = lfcn.hessian();
        H_.setOption("name","nlp_hessian");
        log("SX Hessian function generated");
        
      } else { // MXFunction otherwise

        // Try to cast into MXFunction
        MXFunction F = shared_cast<MXFunction>(F_);
        MXFunction G = shared_cast<MXFunction>(G_);

        // Expressions in F and G
        vector<MX> FG_in;
        MX f, g;
        
        // Convert to MX if cast failed and make sure that they use the same expressions if cast was successful
        if(!G.isNull()){
          FG_in = G.inputExpr();
          g = G.outputExpr(0);
          if(!F.isNull()){ // Both are MXFunction, make sure they use the same variables
            f = substitute(F.outputExpr(),F.inputExpr(),FG_in).front();
          } else { // G_ but not F_ MXFunction
            f = F_.call(FG_in).front();
          }
        } else {
          if(!F.isNull()){ // F_ but not G_ MXFunction
            FG_in = F.inputExpr();
            f = F.outputExpr(0);
            g = G_.call(FG_in).front();
          } else { // None of them MXFunction
            FG_in = G_.symbolicInput();
            g = G_.call(FG_in).front();
            f = F_.call(FG_in).front();
          }
        }
         
        // Lagrange multipliers
        MX lam = msym("lam",g.size1());
      
        // Objective function scaling
        MX sigma = msym("sigma");

        // Lagrangian function
        vector<MX> lfcn_in(parametric_? 4: 3);
        lfcn_in[0] = FG_in.at(0);
        lfcn_in[1] = lam;
        lfcn_in[2] = sigma;
        if (parametric_) lfcn_in[3] = FG_in.at(1);
        MXFunction lfcn(lfcn_in,sigma*f + inner_prod(lam,g));
        lfcn.setOption("verbose",verbose());
        lfcn.setOption("name","nlp_lagrangian");
        lfcn.init();
        if(verbose())
          cout << "MX Lagrangian function generated: algorithm size " << lfcn.getAlgorithmSize() << ", work size = " << lfcn.getWorkSize() << endl;
          
        // Hessian of the Lagrangian
        H_ = lfcn.hessian();
        H_.setOption("name","nlp_hessian");
        log("MX Lagrangian Hessian function generated");
          
      } // SXFunction/MXFunction
    } // constrained/unconstrained
  } // generate_hessian && H_.isNull()
  if(!H_.isNull() && !H_.isInit()) {
    H_.init();
    log("Hessian function initialized");
  }

  // Create a Jacobian if it does not already exists
  bool generate_jacobian = getOption("generate_jacobian");
  if(generate_jacobian && !G_.isNull() && J_.isNull()){
    log("Generating Jacobian");
    J_ = G_.jacobian();
    J_.setOption("name","nlp_jacobian");
    log("Jacobian function generated");
  }
    
  if(!J_.isNull() && !J_.isInit()){
    J_.init();
    log("Jacobian function initialized");
  }

  
  if(!H_.isNull()) {
    if (parametric_) {
      casadi_assert_message(H_.getNumInputs()>=2, "Wrong number of input arguments to H for parametric NLP. Must be at least 2, but got " << G_.getNumInputs());
    } else {
      casadi_assert_message(H_.getNumInputs()>=1, "Wrong number of input arguments to H for non-parametric NLP. Must be at least 1, but got " << G_.getNumInputs() << " instead. Do you perhaps intend to use fixed parameters? Then use the 'parametric' option.");
    }
    casadi_assert_message(H_.getNumOutputs()>=1, "Wrong number of output arguments to H");
    casadi_assert_message(H_.input(0).numel()==nx_,"Inconsistent dimensions");
    casadi_assert_message(H_.output().size1()==nx_,"Inconsistent dimensions");
    casadi_assert_message(H_.output().size2()==nx_,"Inconsistent dimensions");
  }

  if(!J_.isNull()){
    if (parametric_) {
      casadi_assert_message(J_.getNumInputs()==2, "Wrong number of input arguments to J for parametric NLP. Must be at least 2, but got " << G_.getNumInputs());
    } else {
      casadi_assert_message(J_.getNumInputs()==1, "Wrong number of input arguments to J for non-parametric NLP. Must be at least 1, but got " << G_.getNumInputs() << " instead. Do you perhaps intend to use fixed parameters? Then use the 'parametric' option.");
    }
    casadi_assert_message(J_.getNumOutputs()>=1, "Wrong number of output arguments to J");
    casadi_assert_message(J_.input().numel()==nx_,"Inconsistent dimensions");
    casadi_assert_message(J_.output().size2()==nx_,"Inconsistent dimensions");
  }

  if (parametric_) {
    np_ = F_->input(1).size();
    
    if (!G_.isNull()) casadi_assert_message(np_ == G_->input(G_->getNumInputs()-1).size(),"Parametric NLP has inconsistent parameter dimensions. F has got " << np_ << " parameters, while G has got " << G_->input(G_->getNumInputs()-1).size());
    if (!H_.isNull()) casadi_assert_message(np_ == H_->input(H_->getNumInputs()-1).size(),"Parametric NLP has inconsistent parameter dimensions. F has got " << np_ << " parameters, while H has got " << H_->input(H_->getNumInputs()-1).size());
    if (!J_.isNull()) casadi_assert_message(np_ == J_->input(J_->getNumInputs()-1).size(),"Parametric NLP has inconsistent parameter dimensions. F has got " << np_ << " parameters, while J has got " << J_->input(J_->getNumInputs()-1).size());
  } else {
    np_ = 0;
  }
  
  // Allocate space for inputs
  input_.resize(NLP_SOLVER_NUM_IN);
  input(NLP_SOLVER_X0)      =  DMatrix::zeros(nx_);
  input(NLP_SOLVER_LBX)         = -DMatrix::inf(nx_);
  input(NLP_SOLVER_UBX)         =  DMatrix::inf(nx_);
  input(NLP_SOLVER_LBG)         = -DMatrix::inf(ng_);
  input(NLP_SOLVER_UBG)         =  DMatrix::inf(ng_);
  input(NLP_SOLVER_LAM_G0) =  DMatrix::zeros(ng_);
  input(NLP_SOLVER_P)           =  DMatrix::zeros(np_);
  
  // Allocate space for outputs
  output_.resize(NLP_SOLVER_NUM_OUT);
  output(NLP_SOLVER_X)      = DMatrix::zeros(nx_);
  output(NLP_SOLVER_F)       = DMatrix::zeros(1);
  output(NLP_SOLVER_LAM_X)   = DMatrix::zeros(nx_);
  output(NLP_SOLVER_LAM_G)   = DMatrix::zeros(ng_);
  output(NLP_SOLVER_LAM_P)   = DMatrix::zeros(np_);
  output(NLP_SOLVER_G)          = DMatrix::zeros(ng_);
  
  if (hasSetOption("iteration_callback")) {
   callback_ = getOption("iteration_callback");
   if (!callback_.isNull()) {
     if (!callback_.isInit()) callback_.init();
     casadi_assert_message(callback_.getNumOutputs()==1, "Callback function should have one output, a scalar that indicates wether to break. 0 = continue");
     casadi_assert_message(callback_.output(0).size()==1, "Callback function should have one output, a scalar that indicates wether to break. 0 = continue");
     casadi_assert_message(callback_.getNumInputs()==NLP_SOLVER_NUM_OUT, "Callback function should have the output scheme of NLPSolver as input scheme. i.e. " <<NLP_SOLVER_NUM_OUT << " inputs instead of the " << callback_.getNumInputs() << " you provided." );
     for (int i=0;i<NLP_SOLVER_NUM_OUT;i++) {
       if (!callback_.input(i).empty()) {
         casadi_assert_message(callback_.input(i).sparsity()==output(i).sparsity(),
           "Callback function should have the output scheme of NLPSolver as input scheme. " << 
           describeInput(inputScheme_,i) << " was found to be " << callback_.input(i).dimString() << " instead of expected " << output(i).dimString() << "."
         );
       }
       callback_.input(i).setAll(0);
     }
   }
  }
  
  callback_step_ = getOption("iteration_callback_step");

  // Call the initialization method of the base class
  FXInternal::init();
}

void NLPSolverInternal::checkInitialBounds() { 
  if(bool(getOption("warn_initial_bounds"))){
    bool violated = false;
    for (int k=0;k<input(NLP_SOLVER_X0).size();++k) {
      if (input(NLP_SOLVER_X0).at(k)>input(NLP_SOLVER_UBX).at(k)) {
        violated = true;
      }
      if (input(NLP_SOLVER_X0).at(k)<input(NLP_SOLVER_LBX).at(k)) {
        violated = true;
      }
    }
    if (violated) casadi_warning("NLPSolver: The initial guess does not satisfy LBX and UBX. Option 'warn_initial_bounds' controls this warning.");
  }
}
   

  void NLPSolverInternal::reportConstraints(std::ostream &stream) { 
  
    stream << "Reporting NLP constraints" << endl;
    CasADi::reportConstraints(stream,output(NLP_SOLVER_X),input(NLP_SOLVER_LBX),input(NLP_SOLVER_UBX), "decision bounds");
    double tol = 1e-8;
    if (hasOption("constr_viol_tol")) tol = getOption("constr_viol_tol");
    CasADi::reportConstraints(stream,output(NLP_SOLVER_G),input(NLP_SOLVER_LBG),input(NLP_SOLVER_UBG), "constraints",getOption("constr_viol_tol"));
  }


} // namespace CasADi
