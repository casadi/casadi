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

#include "fx_internal.hpp"
#include "../mx/evaluation.hpp"
#include <typeinfo> 
#include "../stl_vector_tools.hpp"
#include "jacobian.hpp"
#include "mx_function.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../matrix/sparsity_tools.hpp"

using namespace std;

namespace CasADi{
  
FXInternal::FXInternal(){
  setOption("name","unnamed_function"); // name of the function
  addOption("sparse",                   OT_BOOLEAN,             true,           "function is sparse");
  addOption("number_of_fwd_dir",        OT_INTEGER,             1,              "number of forward derivatives to be calculated simultanously");
  addOption("number_of_adj_dir",        OT_INTEGER,             1,              "number of adjoint derivatives to be calculated simultanously");
  addOption("verbose",                  OT_BOOLEAN,             false,          "verbose evaluation -- for debugging");
  addOption("store_jacobians",          OT_BOOLEAN,             false,          "keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times");
  addOption("numeric_jacobian",         OT_BOOLEAN,             false,          "Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method");
  addOption("numeric_hessian",          OT_BOOLEAN,             false,          "Calculate Hessians numerically (using directional derivatives) rather than with the built-in method");
  addOption("ad_mode",                  OT_STRING,              "automatic",    "How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)","forward|reverse|automatic");
  addOption("jacobian_generator",       OT_JACOBIANGENERATOR,   GenericType(),  "Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines");
  addOption("sparsity_generator",       OT_SPARSITYGENERATOR,   GenericType(),  "Function that provides sparsity for a given input output block, overrides internal routines");
  addOption("jac_for_sens",             OT_BOOLEAN,             false,          "Create the a Jacobian function and use this to calculate forward sensitivities");
  addOption("user_data",                OT_VOIDPTR,             GenericType(),  "A user-defined field that can be used to identify the function or pass additional information");
  addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "Monitors to be activated","inputs|outputs");
  
  verbose_ = false;
  numeric_jacobian_ = false;
  jacgen_ = 0;
  spgen_ = 0;
  user_data_ = 0;
  jac_for_sens_ = false;
  monitor_inputs_ = false;
  monitor_outputs_ = false;
}

FXInternal::~FXInternal(){
}

void FXInternal::init(){
  verbose_ = getOption("verbose");
  store_jacobians_ = getOption("store_jacobians");
  numeric_jacobian_ = getOption("numeric_jacobian");
  jac_for_sens_ = getOption("jac_for_sens");
  
  // Allocate data for sensitivities (only the method in this class)
  FXInternal::updateNumSens(false);
  
  // Generate storage for generated Jacobians
  if(store_jacobians_){
    jacs_.resize(getNumInputs());
    for(int i=0; i<jacs_.size(); ++i){
      jacs_[i].resize(getNumOutputs());
    }
  }

  // Resize the matrix that holds the sparsity of the Jacobian blocks
  jac_sparsity_compact_.resize(getNumInputs(),vector<CRSSparsity>(getNumOutputs()));
  jac_sparsity_.resize(getNumInputs(),vector<CRSSparsity>(getNumOutputs()));

  // Get the Jacobian generator function, if any
  if(hasSetOption("jacobian_generator")){
    jacgen_ = getOption("jacobian_generator");
  }
  
  // Get the sparsity detector function, if any
  if(hasSetOption("sparsity_generator")){
    spgen_ = getOption("sparsity_generator");
  }

  if(hasSetOption("user_data")){
    user_data_ = getOption("user_data").toVoidPointer();
  }
  
  // Pass monitors
  if(hasSetOption("monitor")){
    const std::vector<std::string>& monitors = getOption("monitor");
    for (std::vector<std::string>::const_iterator it=monitors.begin();it!=monitors.end();it++) {
      monitors_.insert(*it);
    }
  }
  
  monitor_inputs_ = monitored("inputs");
  monitor_outputs_ = monitored("outputs");
  

  // Mark the function as initialized
  is_init_ = true;
}

void FXInternal::updateNumSens(bool recursive){
  nfdir_ = getOption("number_of_fwd_dir");
  nadir_ = getOption("number_of_adj_dir");
  for(vector<FunctionIO>::iterator it=input_.begin(); it!=input_.end(); ++it){
    it->dataF.resize(nfdir_,it->data);
    it->dataA.resize(nadir_,it->data);
  }

  for(vector<FunctionIO>::iterator it=output_.begin(); it!=output_.end(); ++it){
    it->dataF.resize(nfdir_,it->data);
    it->dataA.resize(nadir_,it->data);
  }
}

void FXInternal::print(ostream &stream) const{
  if (getNumInputs()==1) {
    stream << " Input: " << input().dimString() << std::endl;
  } else{
    stream << " Inputs (" << getNumInputs() << "):" << std::endl;
    for (int i=0;i<getNumInputs();i++) {
       stream << "  " << i+1 << ". " << input(i).dimString() << std::endl;
    }
  }
  if (getNumOutputs()==1) {
    stream << " Output: " << output().dimString() << std::endl;
  } else {
    stream << " Outputs (" << getNumOutputs() << "):" << std::endl;
    for (int i=0;i<getNumOutputs();i++) {
       stream << "  " << i+1 << ". " << output(i).dimString() << std::endl;
    }
  }
}

void FXInternal::repr(ostream &stream) const{
  stream << "function(\"" << getOption("name") << "\")";
}

FX FXInternal::hessian(int iind, int oind){
  casadi_error("FXInternal::hessian: hessian not defined for class " << typeid(*this).name());
}

FunctionIO& FXInternal::inputStruct(int i){
  if (i<0 || i>=input_.size()) {
    std::stringstream ss;
    ss <<  "In function " << getOption("name") << ": input " << i << " not in interval [0," << input_.size() << "]";
    if (!isInit()) ss << endl << "Did you forget to initialize?";
    casadi_error(ss.str());
  }

  return input_.at(i);
}

const FunctionIO& FXInternal::inputStruct(int i) const{
  return const_cast<FXInternal*>(this)->inputStruct(i);
}
  
FunctionIO& FXInternal::outputStruct(int i){
  if (i<0 || i>=output_.size()) {
    std::stringstream ss;
    ss << "In function " << getOption("name") << ": output " << i << " not in interval [0," << output_.size() << "]";
    if (!isInit()) ss << endl << "Did you forget to initialize?";
    casadi_error(ss.str());
  }

  return output_.at(i);
}

const FunctionIO& FXInternal::outputStruct(int i) const{
  return const_cast<FXInternal*>(this)->outputStruct(i);
}

void FXInternal::log(const string& msg) const{
  if(verbose()){
    cout << "CasADi log message: " << msg << endl;
  }
}

void FXInternal::log(const string& fcn, const string& msg) const{
  if(verbose()){
    cout << "CasADi log message: In \"" << fcn << "\" --- " << msg << endl;
  }
}

bool FXInternal::verbose() const{
  return verbose_;
}

bool FXInternal::monitored(const string& mod) const{
  return monitors_.count(mod)>0;
}

Matrix<double>& FXInternal::input(int iind){
  return inputStruct(iind).data;
}
    
const Matrix<double>& FXInternal::input(int iind) const{
  return inputStruct(iind).data;
}

Matrix<double>& FXInternal::output(int oind){
  return outputStruct(oind).data;
}
    
const Matrix<double>& FXInternal::output(int oind) const{
  return outputStruct(oind).data;
}

Matrix<double>& FXInternal::fwdSeed(int iind, int dir){
  try{
    return inputStruct(iind).dataF.at(dir);
  } catch(out_of_range&){
    stringstream ss;
    if(inputStruct(iind).dataF.empty()){
      ss << "No forward directions ";
    } else {
      ss << "Forward direction " << dir << " is out of range [0," << inputStruct(iind).dataF.size() << ") ";
    }
    ss << "for function " << getOption("name");
    casadi_error(ss.str());
  }
}
    
const Matrix<double>& FXInternal::fwdSeed(int iind, int dir) const{
  return const_cast<FXInternal*>(this)->fwdSeed(iind,dir);
}

Matrix<double>& FXInternal::fwdSens(int oind, int dir){
  try{
    return outputStruct(oind).dataF.at(dir);
  } catch(out_of_range&){
    stringstream ss;
    if(outputStruct(oind).dataF.empty()){
      ss << "No forward directions ";
    } else {
      ss << "Forward direction " << dir << " is out of range [0," << outputStruct(oind).dataF.size() << ") ";
    }
    ss << "for function " << getOption("name");
    casadi_error(ss.str());
  }
}
    
const Matrix<double>& FXInternal::fwdSens(int oind, int dir) const{
  return const_cast<FXInternal*>(this)->fwdSens(oind,dir);
}

Matrix<double>& FXInternal::adjSeed(int oind, int dir){
  try{
    return outputStruct(oind).dataA.at(dir);
  } catch(out_of_range&){
    stringstream ss;
    if(outputStruct(oind).dataA.empty()){
      ss << "No adjoint directions ";
    } else {
      ss << "Adjoint direction " << dir << " is out of range [0," << outputStruct(oind).dataA.size() << ") ";
    }
    ss << "for function " << getOption("name");
    casadi_error(ss.str());
  }
}
    
const Matrix<double>& FXInternal::adjSeed(int oind, int dir) const{
  return const_cast<FXInternal*>(this)->adjSeed(oind,dir);
}

Matrix<double>& FXInternal::adjSens(int iind, int dir){
  try{
    return inputStruct(iind).dataA.at(dir);
  } catch(out_of_range&){
    stringstream ss;
    if(inputStruct(iind).dataA.empty()){
      ss << "No adjoint directions ";
    } else {
      ss << "Adjoint direction " << dir << " is out of range [0," << inputStruct(iind).dataA.size() << ") ";
    }
    ss << "for function " << getOption("name");
    casadi_error(ss.str());
  }
}
    
const Matrix<double>& FXInternal::adjSens(int iind, int dir) const{
  return const_cast<FXInternal*>(this)->adjSens(iind,dir);
}

void FXInternal::setNumInputs(int num_in){
  return input_.resize(num_in);
}

void FXInternal::setNumOutputs(int num_out){
  return output_.resize(num_out);  
}

int FXInternal::getNumInputs() const{
  return input_.size();
}

int FXInternal::getNumOutputs() const{
  return output_.size();
}

const Dictionary & FXInternal::getStats() const {
  return stats_;
}

GenericType FXInternal::getStat(const string & name) const {
  // Locate the statistic
  Dictionary::const_iterator it = stats_.find(name);

  // Check if found
  if(it == stats_.end()){
    casadi_error("Statistic: " << name << " has not been set." << endl <<  "Note: statistcs are only set after an evaluate call");
  }

  return GenericType(it->second);
}

std::vector<MX> FXInternal::symbolicInput() const{
  vector<MX> ret(getNumInputs());
  assertInit();
  for(int i=0; i<ret.size(); ++i){
    stringstream name;
    name << "x_" << i;
    ret[i] = MX(name.str(),input(i).sparsity());
  }
  return ret;
}

std::vector<SXMatrix> FXInternal::symbolicInputSX() const{
  vector<SXMatrix> ret(getNumInputs());
  assertInit();
  for(int i=0; i<ret.size(); ++i){
    stringstream name;
    name << "x_" << i;
    ret[i] = ssym(name.str(),input(i).sparsity());
  }
  return ret;
}

FX FXInternal::jacobian_switch(const std::vector<std::pair<int,int> >& jblocks){
  if(numeric_jacobian_){
    return numeric_jacobian(jblocks);
  } else {
    if(jacgen_==0){
      // Use internal routine to calculate Jacobian
      return jacobian(jblocks);
    } else {
      // Create a temporary FX instance
      FX tmp;
      tmp.assignNode(this);

      // Use user-provided routine to calculate Jacobian
      return jacgen_(tmp,jblocks,user_data_);
    }
  }
}

FX FXInternal::jacobian(const vector<pair<int,int> >& jblocks){
  return numeric_jacobian(jblocks);
}

FX FXInternal::numeric_jacobian(const vector<pair<int,int> >& jblocks){
  // Symbolic input
  vector<MX> j_in = symbolicInput();
  
  // Nondifferentiated function
  FX fcn;
  fcn.assignNode(this);
  vector<MX> fcn_eval = fcn.call(j_in);
  
  // Less overhead if only jacobian is requested
  if(jblocks.size()==1 && jblocks.front().second>=0){
    return Jacobian(fcn,jblocks.front().second,jblocks.front().first);
  }
  
  // Outputs
  vector<MX> j_out;
  j_out.reserve(jblocks.size());
  for(vector<pair<int,int> >::const_iterator it=jblocks.begin(); it!=jblocks.end(); ++it){
    // If variable index is -1, we want nondifferentiated function output
    if(it->second==-1){
      // Nondifferentiated function
      j_out.push_back(fcn_eval[it->first]);
      
    } else {
      // Create jacobian for block
      Jacobian J(fcn,it->second,it->first);
      
      if(!J.isNull()){
        J.init();
      
        // Evaluate symbolically
        j_out.push_back(J.call(j_in).at(0));
      } else {
        j_out.push_back(MX::sparse(output(it->first).numel(),input(it->second).numel()));
      }
    }
  }
  
  // Create function
  return MXFunction(j_in,j_out);
}

CRSSparsity FXInternal::getJacSparsity(int iind, int oind){
  // Check if we are able to propagate dependencies throught he function
  if(spCanEvaluate(true) || spCanEvaluate(false)){
    
    // Number of nonzero inputs
    int nz_in = input(iind).size();
    
    // Number of nonzero outputs
    int nz_out = output(oind).size();

    // Number of forward sweeps we must make
    int nsweep_fwd = nz_in/bvec_size;
    if(nz_in%bvec_size>0) nsweep_fwd++;
    
    // Number of adjoint sweeps we must make
    int nsweep_adj = nz_out/bvec_size;
    if(nz_out%bvec_size>0) nsweep_adj++;
    
    // Use forward mode?
    bool use_fwd = spCanEvaluate(true) && nsweep_fwd <= nsweep_adj;
    
    // Reset the virtual machine
    spInit(use_fwd);

    // Clear the forward seeds/adjoint sensitivities
    for(int ind=0; ind<getNumInputs(); ++ind){
      vector<double> &v = input(ind).data();
      if(!v.empty()) fill_n(get_bvec_t(v),v.size(),0);
    }

    // Clear the adjoint seeds/forward sensitivities
    for(int ind=0; ind<getNumOutputs(); ++ind){
      vector<double> &v = output(ind).data();
      if(!v.empty()) fill_n(get_bvec_t(v),v.size(),0);
    }
    
    // Get seeds and sensitivities
    bvec_t* seed_v = use_fwd ? get_bvec_t(input(iind).data()) : get_bvec_t(output(oind).data());
    bvec_t* sens_v = use_fwd ? get_bvec_t(output(oind).data()) : get_bvec_t(input(iind).data());
    
    // Number of sweeps needed
    int nsweep = use_fwd ? nsweep_fwd : nsweep_adj;
    
    // The number of zeros in the seed and sensitivity directions
    int nz_seed = use_fwd ? nz_in  : nz_out;
    int nz_sens = use_fwd ? nz_out : nz_in;

    // Input/output index
    int ind_seed = use_fwd ? iind : oind;
    int ind_sens = use_fwd ? oind : iind;

    // Print
    if(verbose()){
      std::cout << "XFunctionInternal::getJacSparsity: using " << (use_fwd ? "forward" : "adjoint") << " mode: ";
      std::cout << nsweep << " sweeps needed for " << nz_seed << " directions" << std::endl;
    }
    
    // Progress
    int progress = -10;

    // Temporary vectors
    std::vector<int> jrow, jcol;
    
    // Loop over the variables, ndir variables at a time
    for(int s=0; s<nsweep; ++s){
      // Print progress
      if(verbose()){
	int progress_new = (s*100)/nsweep;
	// Print when entering a new decade
	if(progress_new / 10 > progress / 10){
	  progress = progress_new;
	  std::cout << progress << " %"  << std::endl;
	}
      }
      
      // Nonzero offset
      int offset = s*bvec_size;

      // Number of local seed directions
      int ndir_local = std::min(bvec_size,nz_seed-offset);

      for(int i=0; i<ndir_local; ++i){
	seed_v[offset+i] |= bvec_t(1)<<i;
      }
      
      // Propagate the dependencies
      spEvaluate(use_fwd);
	    
      // Loop over the nonzeros of the output
      for(int el=0; el<nz_sens; ++el){

	// Get the sparsity sensitivity
	bvec_t spsens = sens_v[el];

	// Clear the sensitivities for the next sweep
	if(!use_fwd){
	  sens_v[el] = 0;
	}
	
	// If there is a dependency in any of the directions
	if(0!=spsens){
	  
	  // Loop over seed directions
	  for(int i=0; i<ndir_local; ++i){
	    
	    // If dependents on the variable
	    if((bvec_t(1) << i) & spsens){
	      // Add to pattern
	      jrow.push_back(el);
	      jcol.push_back(i+offset);
	    }
	  }
	}
      }
      
      // Remove the seeds
      if(use_fwd){
	for(int i=0; i<bvec_size && offset+i<nz_in; ++i){
	  seed_v[offset+i] = 0;
	}
      }
    }
    
    // Set inputs and outputs to zero
    for(int ind=0; ind<getNumInputs(); ++ind) input(ind).setZero();
    for(int ind=0; ind<getNumOutputs(); ++ind) output(ind).setZero();

    // Construct sparsity pattern
    CRSSparsity ret = sp_triplet(nz_out, nz_in,use_fwd ? jrow : jcol, use_fwd ? jcol : jrow);
    
    // Return sparsity pattern
    if(verbose()) std::cout << "XFunctionInternal::getJacSparsity end " << std::endl;
    return ret;

  } else {
    // Dense sparsity by default
    return CRSSparsity(output(oind).size(),input(iind).size(),true);
  }
}

void FXInternal::setJacSparsity(const CRSSparsity& sp, int iind, int oind, bool compact){
  if(compact){
    jac_sparsity_compact_[iind][oind] = sp;
  } else {
    jac_sparsity_[iind][oind] = sp;
  }
}

CRSSparsity& FXInternal::jacSparsity(int iind, int oind, bool compact){
  casadi_assert_message(isInit(),"Function not initialized.");

  // Get a reference to the block
  CRSSparsity& jsp = compact ? jac_sparsity_compact_[iind][oind] : jac_sparsity_[iind][oind];

  // Generate, if null
  if(jsp.isNull()){
    if(compact){
      if(spgen_==0){
        // Use internal routine to determine sparsity
        jsp = getJacSparsity(iind,oind);
      } else {
        // Create a temporary FX instance
        FX tmp;
        tmp.assignNode(this);

        // Use user-provided routine to determine sparsity
        jsp = spgen_(tmp,iind,oind,user_data_);
      }
    } else {
      
      // Get the compact sparsity pattern
      CRSSparsity sp = jacSparsity(iind,oind,true);

      // Enlarge if sparse output 
      if(output(oind).numel()!=sp.size1()){
        casadi_assert(sp.size1()==output(oind).size());
        
        // New row for each old row 
        vector<int> row_map = output(oind).sparsity().getElements();
    
        // Insert rows 
        sp.enlargeRows(output(oind).numel(),row_map);
      }
  
      // Enlarge if sparse input 
      if(input(iind).numel()!=sp.size2()){
        casadi_assert(sp.size2()==input(iind).size());
        
        // New column for each old column
        vector<int> col_map = input(iind).sparsity().getElements();
        
        // Insert columns
        sp.enlargeColumns(input(iind).numel(),col_map);
      }
      
      // Save
      jsp = sp;
    }
  }
  
  // If still null, not dependent
  if(jsp.isNull()){
    jsp = CRSSparsity(output(oind).size(),input(iind).size());
  }
  
  // Return a reference to the block
  return jsp;
}

void FXInternal::getPartition(int iind, int oind, CRSSparsity& D1, CRSSparsity& D2, bool compact, bool symmetric){
  
  // Sparsity pattern with transpose
  CRSSparsity &A = jacSparsity(iind,oind,compact);
  vector<int> mapping;
  CRSSparsity AT = symmetric ? A : A.transpose(mapping);
  mapping.clear();
  
  // Which AD mode?
  bool test_ad_fwd=true, test_ad_adj=true;
  if(getOption("ad_mode") == "forward"){
    test_ad_adj = false;
  } else if(getOption("ad_mode") == "reverse"){
    test_ad_fwd = false;
  } else if(getOption("ad_mode") != "automatic"){
    casadi_error("FXInternal::jac: Unknown ad_mode \"" << getOption("ad_mode") << "\". Possible values are \"forward\", \"reverse\" and \"automatic\".");
  }
  
  // Get seed matrices by graph coloring
  if(symmetric){
    // Star coloring if symmetric
    D1 = A.starColoring();
    
  } else {
    
    // Test unidirectional coloring using forward mode
    if(test_ad_fwd){
      D1 = AT.unidirectionalColoring(A);
    }
      
    // Test unidirectional coloring using reverse mode
    if(test_ad_adj){
      D2 = A.unidirectionalColoring(AT);
    }

    // Use whatever required less colors if we tried both (with preference to forward mode)
    if(test_ad_fwd && test_ad_adj){
      if((D1.size1() <= D2.size1())){
        D2=CRSSparsity();
      } else {
        D1=CRSSparsity();
      }
    }
  }
}

void FXInternal::getFullJacobian(){
  // Quick return if it has been calculated
  if(!full_jacobian_.isNull()) return;

  // We want all the jacobian blocks
  vector<pair<int,int> > jblocks;
  jblocks.reserve((1+getNumInputs())+getNumOutputs());
  for(int oind=0; oind<getNumOutputs(); ++oind){
    for(int iind=-1; iind<getNumInputs(); ++iind){
      if(iind>=0) jacSparsity(iind,oind,false);
      jblocks.push_back(pair<int,int>(oind,iind));
    }
  }
  full_jacobian_ = jacobian_switch(jblocks);
  full_jacobian_.init();
}

void FXInternal::evaluate_switch(int nfdir, int nadir){
  if (monitor_inputs_) {
    std::cout << "FXInternal::evaluate_switch:Inputs:" << std::endl;
    for (int i=0;i<getNumInputs();++i) {
      std::cout << "#" << i << ": " << input(i) << std::endl;
    }
  }
  if(!jac_for_sens_ || (nfdir==0 && nadir==0)){     // Default, directional derivatives
    evaluate(nfdir,nadir);
  } else { // Calculate complete Jacobian and multiply
    // Generate the Jacobian if it does not exist
    if(full_jacobian_.isNull()){
      getFullJacobian();
    }
        
    // Make sure inputs and outputs dense
    for(int iind=0; iind<getNumInputs(); ++iind){
      casadi_assert_message(isDense(input(iind)),"sparse input currently not supported");
    }
    for(int oind=0; oind<getNumOutputs(); ++oind){
      casadi_assert_message(isDense(output(oind)),"sparse output currently not supported");
    }
    
    // Pass inputs to the Jacobian function
    for(int iind=0; iind<getNumInputs(); ++iind){
      full_jacobian_.setInput(input(iind),iind);
    }
    
    // Evaluate Jacobian function
    full_jacobian_.evaluate();

    // Reset forward sensitivities
    for(int dir=0; dir<nfdir; ++dir){
      for(int oind=0; oind<getNumOutputs(); ++oind){
        fwdSens(oind,dir).setAll(0.0);
      }
    }
    
    // Reset adjoint sensitivities
    for(int dir=0; dir<nadir; ++dir){
      for(int iind=0; iind<getNumInputs(); ++iind){
        adjSens(iind,dir).setAll(0.0);
      }
    }
    
    // Jacobian output index
    int oind_jac = 0;

    // Loop over outputs
    for(int oind=0; oind<getNumOutputs(); ++oind){
      // Get undifferentiated result
      full_jacobian_.getOutput(output(oind),oind_jac++);
      
      // Get get contribution from Jacobian block
      for(int iind=0; iind<getNumInputs(); ++iind){
        // Get Jacobian block
        const Matrix<double>& Jblock = full_jacobian_.output(oind_jac++);
        
        // Forward sensitivities
        for(int dir=0; dir<nfdir; ++dir){
          addMultiple(Jblock,fwdSeed(iind,dir).data(),fwdSens(oind,dir).data(),false);
        }
        
        // Adjoint sensitivities
        for(int dir=0; dir<nadir; ++dir){
          addMultiple(Jblock,adjSeed(oind,dir).data(),adjSens(iind,dir).data(),true);
        }
      }
    }
  }
}

void FXInternal::evalSX(const std::vector<SXMatrix>& input, std::vector<SXMatrix>& output, 
			const std::vector<std::vector<SXMatrix> >& fwdSeed, std::vector<std::vector<SXMatrix> >& fwdSens, 
			const std::vector<std::vector<SXMatrix> >& adjSeed, std::vector<std::vector<SXMatrix> >& adjSens,
			bool output_given, bool eliminate_constants){
  casadi_error("FXInternal::evalSX not defined for class " << typeid(*this).name());
}

void FXInternal::evalMX(const std::vector<MX>& input, std::vector<MX>& output, 
			const std::vector<std::vector<MX> >& fwdSeed, std::vector<std::vector<MX> >& fwdSens, 
			const std::vector<std::vector<MX> >& adjSeed, std::vector<std::vector<MX> >& adjSens,
			bool output_given, bool eliminate_constants){
  casadi_error("FXInternal::evalMX not defined for class " << typeid(*this).name());
}

void FXInternal::spEvaluate(bool fwd){
  casadi_error("FXInternal::spEvaluate not defined for class " << typeid(*this).name());
}

} // namespace CasADi

