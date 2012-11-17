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

#include "kinsol_internal.hpp"
#include "symbolic/fx/sx_function_internal.hpp"
#include "symbolic/stl_vector_tools.hpp"
#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/linear_solver_internal.hpp"

using namespace std;
namespace CasADi{

KinsolInternal::KinsolInternal(const FX& f, int nrhs) : ImplicitFunctionInternal(f,nrhs){
  addOption("linear_solver",            OT_STRING, "dense");
  addOption("upper_bandwidth",          OT_INTEGER);
  addOption("lower_bandwidth",          OT_INTEGER);
  addOption("max_krylov",               OT_INTEGER, 0);
  addOption("exact_jacobian",           OT_BOOLEAN, true);
  addOption("iterative_solver",         OT_STRING,"gmres");
  addOption("f_scale",                  OT_REALVECTOR);
  addOption("u_scale",                  OT_REALVECTOR);
  addOption("pretype",                  OT_STRING, "none","","none|left|right|both");
  addOption("use_preconditioner",       OT_BOOLEAN, false); // precondition an iterative solver
  addOption("constraints",              OT_INTEGERVECTOR);
  addOption("strategy",                 OT_STRING, "none", "Globalization strateg","none|linesearch");
  addOption("linear_solver_creator",    OT_LINEARSOLVER, GenericType(), "User-defined linear solver class");
  addOption("linear_solver_options",    OT_DICTIONARY, GenericType(), "Options to be passed to the linear solver");
  addOption("disable_internal_warnings",   OT_BOOLEAN,false, "Disable KINSOL internal warning messages");
  addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "", "eval_f|eval_djac", true);
    
  mem_ = 0;
  u_ = 0;
  u_scale_ = 0;
  f_scale_ = 0;
  disable_internal_warnings_ = false;
}

KinsolInternal* KinsolInternal::clone() const{
  // Return a deep copy
  FX f = shared_cast<FX>(f_.clone());
  KinsolInternal* node = new KinsolInternal(f,nrhs_);
  node->setOption(dictionary());
  node->setJacobian(shared_cast<FX>(J_.clone()));
  node->setLinearSolver(shared_cast<LinearSolver>(linsol_.clone()));
  
  if(isInit())
    node->init();
  return node;
}

KinsolInternal::~KinsolInternal(){
  if(u_) N_VDestroy_Serial(u_);
  if(u_scale_) N_VDestroy_Serial(u_scale_);
  if(f_scale_) N_VDestroy_Serial(f_scale_);
  if(mem_) KINFree(&mem_);
}

void KinsolInternal::init(){
  ImplicitFunctionInternal::init();
  
  // Read options
  if(getOption("strategy")=="linesearch"){
    strategy_ = KIN_LINESEARCH;
  } else {
    casadi_assert(getOption("strategy")=="none");
    strategy_ = KIN_NONE;
  }

  // Get the linear solver creator function
  if(linsol_.isNull() && hasSetOption("linear_solver_creator")){
    linearSolverCreator linear_solver_creator = getOption("linear_solver_creator");
  
    // Allocate an NLP solver
    linsol_ = linear_solver_creator(CRSSparsity());
  
    // Pass options
    if(hasSetOption("linear_solver_options")){
      const Dictionary& linear_solver_options = getOption("linear_solver_options");
      linsol_.setOption(linear_solver_options);
    }
  }

  // Return flag
  int flag;
  
  // Use exact Jacobian?
  bool exact_jacobian = getOption("exact_jacobian");

  // Generate Jacobian if not provided
  if(J_.isNull()) J_ = f_.jacobian(0,0);
  J_.init();
  
  // Initialize the linear solver, if provided
  if(!linsol_.isNull()){
    linsol_.setSparsity(J_.output().sparsity());
    linsol_.init();
  }
  
  // Allocate N_Vectors
  if(u_) N_VDestroy_Serial(u_);
  if(u_scale_) N_VDestroy_Serial(u_scale_);
  if(f_scale_) N_VDestroy_Serial(f_scale_);
  u_ = N_VNew_Serial(N_);
  u_scale_ = N_VNew_Serial(N_);
  f_scale_ = N_VNew_Serial(N_);
  
  // Set scaling factors on variables
  if(hasSetOption("u_scale")){
    const vector<double>& u_scale = getOption("u_scale");
    casadi_assert(u_scale.size()==NV_LENGTH_S(u_scale_));
    copy(u_scale.begin(),u_scale.end(),NV_DATA_S(u_scale_));
  } else {
    N_VConst(1.0,u_scale_);
  }
  
  // Set scaling factors on equations
  if(hasSetOption("f_scale")){
    const vector<double>& f_scale = getOption("f_scale");
    casadi_assert(f_scale.size()==NV_LENGTH_S(f_scale_));
    copy(f_scale.begin(),f_scale.end(),NV_DATA_S(f_scale_));
  } else {
    N_VConst(1.0,f_scale_);
  }
  
  // Create KINSOL memory block
  if(mem_) KINFree(&mem_);
  mem_ = KINCreate();
  
  // KINSOL bugfix
  KINMem kin_mem = KINMem(mem_);
  kin_mem->kin_inexact_ls = FALSE;
    
  // Set optional inputs
  flag = KINSetUserData(mem_, this);
  casadi_assert_message(flag==KIN_SUCCESS, "KINSetUserData");

  // Disable internal warning messages?
  disable_internal_warnings_ = getOption("disable_internal_warnings");
  
  // Set error handler function
  flag = KINSetErrHandlerFn(mem_, ehfun_wrapper, this);
  casadi_assert_message(flag==KIN_SUCCESS, "KINSetErrHandlerFn");
  
  // Initialize KINSOL
  flag = KINInit(mem_,func_wrapper, u_);
  casadi_assert(flag==KIN_SUCCESS);

  // Set constraints
  if(hasSetOption("constraints")){
    // Get the user-set constraints
    const vector<int>& u_c = getOption("constraints");
    casadi_assert(u_c.size()==N_);
    
    // Copy to a temporary N_Vector
    N_Vector constraints = N_VNew_Serial(N_);
    copy(u_c.begin(),u_c.end(),NV_DATA_S(constraints));
    
    // Pass to KINSOL
    flag = KINSetConstraints(mem_, constraints);
    casadi_assert(flag==KIN_SUCCESS);
    
    // Free the temporary vector
    N_VDestroy_Serial(constraints);
  }

  // attach a linear solver
  if(getOption("linear_solver")=="dense"){
    // Dense jacobian
    flag = KINDense(mem_, N_);
    casadi_assert_message(flag==KIN_SUCCESS, "KINDense");
    
    if(exact_jacobian){
      flag = KINDlsSetDenseJacFn(mem_, djac_wrapper);
      casadi_assert_message(flag==KIN_SUCCESS, "KINDlsSetDenseJacFn");
    }
    
  } else if(getOption("linear_solver")=="banded") {
    // Banded jacobian
    flag = KINBand(mem_, N_, getOption("upper_bandwidth").toInt(), getOption("lower_bandwidth").toInt());
    casadi_assert_message(flag==KIN_SUCCESS, "KINBand");
    
    if(exact_jacobian){
      flag = KINDlsSetBandJacFn(mem_, bjac_wrapper);
      casadi_assert_message(flag==KIN_SUCCESS, "KINDlsBandJacFn");
    }
    
  } else if(getOption("linear_solver")=="iterative") {
    // Sparse (iterative) solver  
    // Max dimension of the Krylov space
    int maxl = getOption("max_krylov").toInt();

    // Attach the sparse solver  
    if(getOption("iterative_solver")=="gmres"){
      flag = KINSpgmr(mem_, maxl);
      casadi_assert_message(flag==KIN_SUCCESS, "KINSpgmr");
    } else if(getOption("iterative_solver")=="bcgstab") {
      flag = KINSpbcg(mem_, maxl);
      casadi_assert_message(flag==KIN_SUCCESS, "KINSpbcg");
    } else if(getOption("iterative_solver")=="tfqmr") {
      flag = KINSptfqmr(mem_, maxl);
      casadi_assert_message(flag==KIN_SUCCESS, "KINSptfqmr");
    } else {
      throw CasadiException("KINSOL: Unknown sparse solver");
    }
    
    // Attach functions for jacobian information
    if(exact_jacobian){
      flag = KINSpilsSetJacTimesVecFn(mem_, jtimes_wrapper);
      casadi_assert_message(flag==KIN_SUCCESS, "KINSpilsSetJacTimesVecFn");
    }
    
    // Add a preconditioner
    if(bool(getOption("use_preconditioner"))){
      // Make sure that a Jacobian has been provided
      casadi_assert_message(!J_.isNull(),"No Jacobian has been provided");

      // Make sure that a linear solver has been providided
      casadi_assert_message(!linsol_.isNull(), "No linear solver has been provided.");

      // Pass to IDA
      flag = KINSpilsSetPreconditioner(mem_, psetup_wrapper, psolve_wrapper);
      casadi_assert(flag==KIN_SUCCESS);
    }
    
  } else if(getOption("linear_solver")=="user_defined") {
    // Make sure that a Jacobian has been provided
    casadi_assert(!J_.isNull());

    // Make sure that a linear solver has been providided
    casadi_assert(!linsol_.isNull());

    // Set fields in the IDA memory
    KINMem kin_mem = KINMem(mem_);
    kin_mem->kin_lmem   = this;
    kin_mem->kin_lsetup = lsetup_wrapper;
    kin_mem->kin_lsolve = lsolve_wrapper;
    kin_mem->kin_setupNonNull = TRUE;
    
  } else {
    throw CasadiException("Unknown linear solver ");
  }
  
  // Set stop criterion
  if(hasSetOption("abstol")){
    flag = KINSetFuncNormTol(mem_,getOption("abstol"));
    casadi_assert(flag==KIN_SUCCESS);
  }
}

void KinsolInternal::evaluate(int nfdir, int nadir){
  casadi_assert(nfdir<=nfdir_fcn_);
  casadi_assert(nadir<=nadir_fcn_);
  
  // Reset the counters
  t_func_ = 0;
  t_jac_ = 0;

  if(verbose()){
    cout << "KinsolInternal::evaluate: Initial guess = " << output(0).data() << endl;
  }
  
  // Get the initial guess
  output().get(NV_DATA_S(u_));
  
  // Solve the nonlinear system of equations
  int flag = KINSol(mem_, u_, strategy_, u_scale_, f_scale_);
  if (flag<KIN_SUCCESS) kinsol_error("evaluate",flag);
  
  // Warn if not successful return
  if(verbose()){
    if(flag!=KIN_SUCCESS) kinsol_error("evaluate",flag,false);
  }
  
  // Save the solution
  output(0).set(NV_DATA_S(u_));

  // Save auxillary outputs
  for(int i=1; i<getNumOutputs(); ++i){
    output(i).set(f_.output(i));
  }
  
  // Print solution
  if(verbose()){
    cout << "KinsolInternal::evaluate: solution = " << output(0).data() << endl;
  }
  
  // End of function if no sensitivities
  if(nfdir==0 && nadir==0)
    return;
  
  // Make sure that a linear solver has been provided
  casadi_assert_message(!linsol_.isNull(),"Sensitivities of an implicit function requires a provided linear solver");
  casadi_assert_message(!J_.isNull(),"Sensitivities of an implicit function requires an exact Jacobian");

  // Factorize the linear system TODO: make optional
  psetup(u_, u_scale_, 0, f_scale_, 0, 0);

  // Pass inputs to function
  f_.setInput(output(0),0);
  for(int i=0; i<getNumInputs(); ++i)
    f_.input(i+1).set(input(i));

  // Pass input seeds to function
  for(int dir=0; dir<nfdir; ++dir){
    f_.fwdSeed(0,dir).setZero();
    for(int i=0; i<getNumInputs(); ++i){
      f_.fwdSeed(i+1,dir).set(fwdSeed(i,dir));
    }
  }
  
  // Solve for the adjoint seeds
  for(int dir=0; dir<nadir; ++dir){
    // Negate adjoint seed and pass to function
    Matrix<double>& faseed = f_.adjSeed(0,dir);
    faseed.set(adjSeed(0,dir));
    for(vector<double>::iterator it=faseed.begin(); it!=faseed.end(); ++it){
      *it = -*it;
    }
    
    // Solve the transposed linear system
    linsol_.solve(&faseed.front(),1,true);

    // Set auxillary adjoint seeds
    for(int oind=1; oind<getNumOutputs(); ++oind){
      f_.adjSeed(oind,dir).set(adjSeed(oind,dir));
    }
  }
  
  // Evaluate
  f_.evaluate(nfdir,nadir);
  
  // Get the forward sensitivities
  for(int dir=0; dir<nfdir; ++dir){
    // Negate intermediate result and copy to output
    Matrix<double>& fsens = fwdSens(0,dir);
    fsens.set(f_.fwdSens(0,dir));
    for(vector<double>::iterator it=fsens.begin(); it!=fsens.end(); ++it){
      *it = -*it;
    }
    
    // Solve the linear system
    linsol_.solve(&fsens.front());
  }
  
  // Get auxillary forward sensitivities
  if(getNumOutputs()>1){
    // Pass the seeds to the implicitly defined variables
    for(int dir=0; dir<nfdir; ++dir){
      f_.fwdSeed(0,dir).set(fwdSens(0,dir));
    }
    
    // Evaluate
    f_.evaluate(nfdir);
  
    // Get the sensitivities
    for(int dir=0; dir<nfdir; ++dir){
      for(int oind=1; oind<getNumOutputs(); ++oind){
        fwdSens(oind,dir).set(f_.fwdSens(oind,dir));
      }
    }
  }
  
  // Get the adjoint sensitivities
  for(int dir=0; dir<nadir; ++dir){
    for(int i=0; i<getNumInputs(); ++i){
      f_.adjSens(i+1,dir).get(adjSens(i,dir));
    }
  }
}

KinsolSolver KinsolInternal::jac(int iind, int oind){
  return jac(vector<int>(1,iind),oind);
}

KinsolSolver KinsolInternal::jac(const std::vector<int> iind, int oind){
  // Single output
  casadi_assert(oind==0);
  
  // Get the function
  SXFunction f = shared_cast<SXFunction>(f_);
  casadi_assert(!f.isNull());
  
  // Get the jacobians
  Matrix<SX> Jz = f.jac(0,0);

  // Number of equations
  int nz = f.input(0).numel();

  // All variables
  vector<Matrix<SX> > f_in(f.getNumInputs());
  f_in[0] = f.inputExpr(0);
  for(int i=1; i<f.getNumInputs(); ++i)
    f_in[i] = f.inputExpr(i);

  // Augmented nonlinear equation
  Matrix<SX> F_aug = f.outputExpr(0);

  // Number of right hand sides
  int nrhs = 1;
  
  // Augment variables and equations
  for(vector<int>::const_iterator it=iind.begin(); it!=iind.end(); ++it){

    // Get the jacobian
    Matrix<SX> Jx = f.jac(*it+1,0);
    
    // Number of variables
    int nx = f.input(*it+1).numel();

    // Sensitivities
    Matrix<SX> dz_dx = ssym("dz_dx", nz, nx);
    
    // Derivative equation
    Matrix<SX> f_der = mul(Jz,dz_dx) + Jx;

    // Append variables
    f_in[0].append(vec(dz_dx));
      
    // Augment nonlinear equation
    F_aug.append(vec(f_der));
    
    // Number of right hand sides
    nrhs += nx;
  }

  // Augmented nonlinear equation
  SXFunction f_aug(f_in, F_aug);

  // Create KINSOL instance
  KinsolSolver augsol(f_aug, nrhs);
  
  // Return the created solver
  return augsol;
}

void KinsolInternal::func(N_Vector u, N_Vector fval){
  // Get time
  time1_ = clock();

  // Pass input
  f_.setInput(NV_DATA_S(u),0);
  for(int i=0; i<getNumInputs(); ++i)
    f_.setInput(input(i),i+1);
  
  // Evaluate
  f_.evaluate();
    
  if(monitored("eval_f")){
      cout << "f = " << f_.output() << endl;
  }
    
  // Get results
  f_.getOutput(NV_DATA_S(fval));

  
  // Get a referebce to the nonzeros of the function
  const vector<double>& fdata = f_.output().data();
  
  // Make sure that all entries of the linear system are valid
  for(int k=0; k<fdata.size(); ++k){
	try{
	  casadi_assert_message(!isnan(fdata[k]),"Nonzero " << k << " is not-a-number");
	  casadi_assert_message(!isinf(fdata[k]),"Nonzero " << k << " is infinite");
	} catch(exception& ex){
	  stringstream ss;
	  ss << ex.what() << endl;
	  if(verbose()){
		ss << "u = " << f_.input() << endl;
		
		// Print the expression for f[Jrow] if f is an SXFunction instance
		SXFunction f_sx = shared_cast<SXFunction>(f_);
		if(!f_sx.isNull()){
		  f_sx.print(ss);
		  ss << f_sx->work_ << endl;
		  ss << "Equation " << k << " = " << f_sx.outputExpr(0).at(k) << endl;
		}
	  }

	  throw CasadiException(ss.str());
	}
  }
  
  // Log time
  time2_ = clock();
  t_func_ += double(time2_-time1_)/CLOCKS_PER_SEC;
}

int KinsolInternal::func_wrapper(N_Vector u, N_Vector fval, void *user_data){
try{
    casadi_assert(user_data);
    KinsolInternal *this_ = (KinsolInternal*)user_data;
    this_->func(u,fval);
    return 0;
  } catch(exception& e){
    cerr << "func failed: " << e.what() << endl;
    return 1;
  }
}

int KinsolInternal::djac_wrapper(long N, N_Vector u, N_Vector fu, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2){
  try{
    casadi_assert(user_data);
    KinsolInternal *this_ = (KinsolInternal*)user_data;
    this_->djac(N, u, fu, J, tmp1, tmp2);
    return 0;
  } catch(exception& e){
    cerr << "djac failed: " << e.what() << endl;;
    return 1;
  }
}

void KinsolInternal::djac(long N, N_Vector u, N_Vector fu, DlsMat J, N_Vector tmp1, N_Vector tmp2){
  // Get time
  time1_ = clock();

  // Pass inputs to the jacobian function
  J_.setInput(NV_DATA_S(u),0);
  for(int i=0; i<getNumInputs(); ++i)
    J_.setInput(input(i),i+1);

  // Evaluate
  J_.evaluate();

  if(monitored("eval_djac")){
      cout << "djac = " << J_.output() << endl;
  }
  
  // Get sparsity and non-zero elements
  const vector<int>& rowind = J_.output().rowind();
  const vector<int>& col = J_.output().col();
  const vector<double>& val = J_.output().data();

  // Loop over rows
  for(int i=0; i<rowind.size()-1; ++i){
    // Loop over non-zero entries
    for(int el=rowind[i]; el<rowind[i+1]; ++el){
      // Get column
      int j = col[el];
      
      // Set the element
      DENSE_ELEM(J,i,j) = val[el];
    }
  }
  
  if(monitored("eval_djac")){
      cout << "djac = ";
      PrintMat(J);
  }
  
  // Log time duration
  time2_ = clock();
  t_jac_ += double(time2_-time1_)/CLOCKS_PER_SEC;
}

int KinsolInternal::bjac_wrapper(long N, long mupper, long mlower, N_Vector u, N_Vector fu, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2){
  try{
    casadi_assert(user_data);
    KinsolInternal *this_ = (KinsolInternal*)user_data;
    this_->bjac(N, mupper, mlower, u, fu, J, tmp1, tmp2);
    return 0;
  } catch(exception& e){
    cerr << "bjac failed: " << e.what() << endl;;
    return 1;
  }
}

void KinsolInternal::bjac(long N, long mupper, long mlower, N_Vector u, N_Vector fu, DlsMat J, N_Vector tmp1, N_Vector tmp2){
  // Get time
  time1_ = clock();

  // Pass inputs to the jacobian function
  J_.setInput(NV_DATA_S(u),0);
  for(int i=0; i<getNumInputs(); ++i)
    J_.setInput(input(i),i+1);

  // Evaluate
  J_.evaluate();
  
  // Get sparsity and non-zero elements
  const vector<int>& rowind = J_.output().rowind();
  const vector<int>& col = J_.output().col();
  const vector<double>& val = J_.output().data();

  // Loop over rows
  for(int i=0; i<rowind.size()-1; ++i){
    // Loop over non-zero entries
    for(int el=rowind[i]; el<rowind[i+1]; ++el){
      // Get column
      int j = col[el];
      
      // Set the element
      if(i-j>=-mupper && i-j<=mlower)
        BAND_ELEM(J,i,j) = val[el];
    }
  }
  
  // Log time duration
  time2_ = clock();
  t_jac_ += double(time2_-time1_)/CLOCKS_PER_SEC;
}

int KinsolInternal::jtimes_wrapper(N_Vector v, N_Vector Jv, N_Vector u, int* new_u, void *user_data){
  try{
    casadi_assert(user_data);
    KinsolInternal *this_ = (KinsolInternal*)user_data;
    this_->jtimes(v,Jv,u,new_u);
    return 0;
  } catch(exception& e){
    cerr << "jtimes failed: " << e.what() << endl;;
    return 1;
  }
}

void KinsolInternal::jtimes(N_Vector v, N_Vector Jv, N_Vector u, int* new_u){
  // Get time
  time1_ = clock();

  // Pass inputs
  f_.setInput(NV_DATA_S(u),0);
  for(int i=0; i<getNumInputs(); ++i)
    f_.setInput(input(i),i+1);

  // Pass input seeds
  f_.setFwdSeed(NV_DATA_S(v),0);
  for(int i=0; i<getNumInputs(); ++i)
    f_.fwdSeed(i+1).setZero();
  
  // Evaluate
  f_.evaluate(1,0);

  // Get the output seeds
  f_.getFwdSens(NV_DATA_S(Jv));
  
  // Log time duration
  time2_ = clock();
  t_jac_ += double(time2_-time1_)/CLOCKS_PER_SEC;
}

int KinsolInternal::psetup_wrapper(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, void* user_data, N_Vector tmp1, N_Vector tmp2){
  try{
    casadi_assert(user_data);
    KinsolInternal *this_ = (KinsolInternal*)user_data;
    this_->psetup(u, uscale, fval, fscale, tmp1, tmp2);
    return 0;
  } catch(exception& e){
    cerr << "psetup failed: " << e.what() << endl;;
    return 1;
  }
}

void KinsolInternal::psetup(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector tmp1, N_Vector tmp2){
  // Get time
  time1_ = clock();

  // Pass inputs
  J_.setInput(NV_DATA_S(u),0);
  for(int i=0; i<getNumInputs(); ++i)
    J_.setInput(input(i),i+1);

  // Evaluate jacobian
  J_.evaluate();

  // Get a referebce to the nonzeros of Jacobian
  const vector<double>& Jdata = J_.output().data();
  
  // Make sure that all entries of the linear system are valid
  for(int k=0; k<Jdata.size(); ++k){
	try{
	  casadi_assert_message(!isnan(Jdata[k]),"Nonzero " << k << " is not-a-number");
	  casadi_assert_message(!isinf(Jdata[k]),"Nonzero " << k << " is infinite");
	} catch(exception& ex){
	  stringstream ss;
	  ss << ex.what() << endl;
	  
	  if(verbose()){
	  
		// Print inputs
		ss << "Input vector is " << J_.input().data() << endl;
		
		// Get the row
		int Jrow = J_.output().sparsity().getRow().at(k);

		// Get the column
		int Jcol = J_.output().sparsity().col(k);
		
		// Which equation
		ss << "This corresponds to the derivative of equation " << Jrow << " with respect to the variable " << Jcol << "." << endl;
		
		// Print the expression for f[Jrow] if f is an SXFunction instance
		SXFunction f_sx = shared_cast<SXFunction>(f_);
		if(!f_sx.isNull()){
		  ss << "Variable " << Jcol << " = " << f_sx.inputExpr(0).at(Jcol) << endl;
		  ss << "Equation " << Jrow << " = " << f_sx.outputExpr(0).at(Jrow) << endl;
		}
		
		// Print the expression for J[k] if J is an SXFunction instance
		SXFunction J_sx = shared_cast<SXFunction>(J_);
		if(!J_sx.isNull()){
		  ss << "J[" << Jrow << "," << Jcol << "] = " << J_sx.outputExpr(0).at(k) << endl;
		}
	  }

	  throw CasadiException(ss.str());
	}
  }
  
  // Log time duration
  time2_ = clock();
  t_lsetup_jac_ += double(time2_-time1_)/CLOCKS_PER_SEC;

  // Pass non-zero elements, scaled by -gamma, to the linear solver
  linsol_.setInput(J_.output(),0);

  // Prepare the solution of the linear system (e.g. factorize) -- only if the linear solver inherits from LinearSolver
  linsol_.prepare();

  // Log time duration
  time1_ = clock();
  t_lsetup_fac_ += double(time1_-time2_)/CLOCKS_PER_SEC;
}

int KinsolInternal::psolve_wrapper(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector v, void* user_data, N_Vector tmp){
  try{
    casadi_assert(user_data);
    KinsolInternal *this_ = (KinsolInternal*)user_data;
    this_->psolve(u, uscale, fval, fscale, v, tmp);
    return 0;
  } catch(exception& e){
    cerr << "psolve failed: " << e.what() << endl;;
    return 1;
  }
}

void KinsolInternal::psolve(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector v, N_Vector tmp){
  // Get time
  time1_ = clock();

  // Solve the factorized system 
  linsol_.solve(NV_DATA_S(v));
  
  // Log time duration
  time2_ = clock();
  t_lsolve_ += double(time2_-time1_)/CLOCKS_PER_SEC;
}

FX KinsolInternal::getJacobian(){
  return J_;
}

LinearSolver KinsolInternal::getLinearSolver(){
  return linsol_;
}

void KinsolInternal::setJacobian(const FX& jac){
  J_ = jac;
}

void KinsolInternal::setLinearSolver(const LinearSolver& linsol){
  linsol_ = linsol;
}

int KinsolInternal::lsetup_wrapper(KINMem kin_mem){
 try{
    KinsolInternal *this_ = (KinsolInternal*)(kin_mem->kin_lmem);
    casadi_assert(this_);
    this_->lsetup(kin_mem);
    return 0;
  } catch(exception& e){
    cerr << "lsetup failed: " << e.what() << endl;;
    return -1;
  }
}

void KinsolInternal::lsetup(KINMem kin_mem){
  N_Vector u =  kin_mem->kin_uu;
  N_Vector uscale = kin_mem->kin_uscale;
  N_Vector fval = kin_mem->kin_fval;
  N_Vector fscale = kin_mem->kin_fscale;
  N_Vector tmp1 = kin_mem->kin_vtemp1;
  N_Vector tmp2 = kin_mem->kin_vtemp2;
  psetup(u, uscale, fval, fscale, tmp1, tmp2);
}

int KinsolInternal::lsolve_wrapper(KINMem kin_mem, N_Vector x, N_Vector b, double *res_norm){
 try{
    KinsolInternal *this_ = (KinsolInternal*)(kin_mem->kin_lmem);
    casadi_assert(this_);
    this_->lsolve(kin_mem,x,b,res_norm);
    return 0;
  } catch(exception& e){
    cerr << "lsolve failed: " << e.what() << endl;;
    return -1;
  }
}


void KinsolInternal::lsolve(KINMem kin_mem, N_Vector x, N_Vector b, double *res_norm){
  // Get vectors
  N_Vector u =  kin_mem->kin_uu;
  N_Vector uscale = kin_mem->kin_uscale;
  N_Vector fval = kin_mem->kin_fval;
  N_Vector fscale = kin_mem->kin_fscale;
  N_Vector tmp1 = kin_mem->kin_vtemp1;
  N_Vector tmp2 = kin_mem->kin_vtemp2;
  
  // Solve the linear system
  N_VScale(1.0, b, x);
  psolve(u,uscale,fval,fscale,x,tmp1);

  // Calculate residual
  jtimes(x, tmp2, u, 0);

  // Calculate the error in residual norm
  N_VLinearSum(1.0, b, -1.0, tmp2, tmp1);
  *res_norm = sqrt(N_VDotProd(tmp1, tmp1));
}

void KinsolInternal::ehfun_wrapper(int error_code, const char *module, const char *function, char *msg, void *eh_data){
try{
    casadi_assert(eh_data);
    KinsolInternal *this_ = (KinsolInternal*)eh_data;
    this_->ehfun(error_code,module,function,msg);        
  } catch(exception& e){
    cerr << "ehfun failed: " << e.what() << endl;
  }
}
  
void KinsolInternal::ehfun(int error_code, const char *module, const char *function, char *msg){
  if(!disable_internal_warnings_){
    cerr << msg << endl;
  }
}

map<int,Message> KinsolInternal::flagmap = KinsolInternal::calc_flagmap();

void KinsolInternal::kinsol_error(const string& module, int flag, bool fatal){
  // Find the error
  map<int,Message>::const_iterator it = flagmap.find(flag);
  
  stringstream ss;
  if(it == flagmap.end()){
    ss << "Unknown " << (fatal? "error" : "warning") <<" (" << flag << ") from module \"" << module << "\".";
  } else {
    ss << "Module \"" << module << "\" returned flag \"" << it->second.first << "\"." << endl;
    ss << "The description of this flag is: " << endl;
    ss << "\"" << it->second.second << "\"" << endl;
  }
  ss << "Consult KINSOL documentation for more information.";
  if (fatal) {
    casadi_error(ss.str())
  } else {
    casadi_warning(ss.str());
  }
}


map<int,Message> KinsolInternal::calc_flagmap(){
  
  map<int, Message > f;
  f[KIN_SUCCESS] = Message("KIN_SUCCES","KINSol succeeded; the scaled norm of F(u) is less than fnormtol");
  f[KIN_INITIAL_GUESS_OK] = Message("KIN_INITIAL_GUESS_OK","The guess u = u0 satisfied the system F(u) = 0 within the tolerances specified.");
  f[KIN_STEP_LT_STPTOL] = Message("KIN_STEP_LT_STPTOL","KINSol stopped based on scaled step length. This means that the current iterate may be an approximate solution of the given nonlinear system, but it is also quite possible that the algorithm is 'stalled' (making insufficient progress) near an invalid solution, or that the scalar scsteptol is too large.");
  f[KIN_MEM_NULL] = Message("KIN_MEM_NULL","The kinsol memory block pointer was NULL.");
  f[KIN_ILL_INPUT] = Message("KIN_ILL_INPUT","An input parameter was invalid.");
  f[KIN_NO_MALLOC] = Message("KIN_NO_MALLOC","The kinsol memory was not allocated by a call to KINCreate.");
  f[KIN_LINESEARCH_NONCONV] = Message("KIN_LINESEARCH_NONCONV","The line search algorithm was unable to find an iterate sufficiently distinct from the current iterate, or could not find an iterate satisfying the sufficient decrease condition. Failure to satisfy the sufficient decrease condition could mean the current iterate is 'close' to an approximate solution of the given nonlinear system, the difference approximation of the matrix-vector product J(u)v is inaccurate, or the real scalar scsteptol is too large.");
  f[KIN_MAXITER_REACHED] = Message("KIN_MAXITER_REACHED","The maximum number of nonlinear iterations has been reached.");
  f[KIN_MXNEWT_5X_EXCEEDED] = Message("KIN_MXNEWT_5X_EXCEEDED","Five consecutive steps have been taken that satisfy the inequality  || D_u p ||_L2 > 0.99 mxnewtstep, where p denotes the current step and mxnewtstep is a scalar upper bound on the scaled step length. Such a failure may mean that || D_F F(u)||_L2 asymptotes from above to a positive value, or the real scalar mxnewtstep is too small.");
  f[KIN_LINESEARCH_BCFAIL] = Message("KIN_LINESEARCH_BCFAIL","The line search algorithm was unable to satisfy the “beta-condition” for MXNBCF +1 nonlinear iterations (not necessarily consecutive), which may indicate the algorithm is making poor progress.");
  f[KIN_LINSOLV_NO_RECOVERY] = Message("KIN_LINSOLV_NO_RECOVERY","The user-supplied routine psolve encountered a recoverable error, but the preconditioner is already current.");
  f[KIN_LINIT_FAIL] = Message("KIN_LINIT_FAIL","The linear solver initialization routine (linit) encountered an error.");
  f[KIN_LSETUP_FAIL] = Message("KIN_LSETUP_FAIL","The user-supplied routine pset (used to set up the preconditioner data) encountered an unrecoverable error.");
  f[KIN_LSOLVE_FAIL] = Message("KIN_LSOLVE_FAIL","Either the user-supplied routine psolve (used to to solve the preconditioned linear system) encountered an unrecoverable error, or the linear solver routine (lsolve) encountered an error condition.");
  f[KIN_SYSFUNC_FAIL] = Message("KIN_SYSFUNC_FAIL","The system function failed in an unrecoverable manner.");
  f[KIN_FIRST_SYSFUNC_ERR] = Message("KIN_FIRST_SYSFUNC_ERR","The system function failed recoverably at the first call.");
  f[KIN_REPTD_SYSFUNC_ERR] = Message("KIN_REPTD_SYSFUNC_ERR","The system function had repeated recoverable errors. No recovery is possible.");
  return f;
}

} // namespace CasADi

