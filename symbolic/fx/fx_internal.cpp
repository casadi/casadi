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
#include "../mx/evaluation_mx.hpp"
#include <typeinfo> 
#include "../stl_vector_tools.hpp"
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
  addOption("max_number_of_fwd_dir",    OT_INTEGER,             optimized_num_dir,  "Allow \"number_of_fwd_dir\" to grow until it reaches this number");
  addOption("max_number_of_adj_dir",    OT_INTEGER,             optimized_num_dir,  "Allow \"number_of_adj_dir\" to grow until it reaches this number");
  addOption("verbose",                  OT_BOOLEAN,             false,          "verbose evaluation -- for debugging");
  addOption("store_jacobians",          OT_BOOLEAN,             false,          "keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times");
  addOption("numeric_jacobian",         OT_BOOLEAN,             false,          "Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method");
  addOption("numeric_hessian",          OT_BOOLEAN,             false,          "Calculate Hessians numerically (using directional derivatives) rather than with the built-in method");
  addOption("ad_mode",                  OT_STRING,              "automatic",    "How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)","forward|reverse|automatic");
  addOption("jacobian_generator",       OT_JACOBIANGENERATOR,   GenericType(),  "Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines");
  addOption("sparsity_generator",       OT_SPARSITYGENERATOR,   GenericType(),  "Function that provides sparsity for a given input output block, overrides internal routines");
  addOption("user_data",                OT_VOIDPTR,             GenericType(),  "A user-defined field that can be used to identify the function or pass additional information");
  addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "Monitors to be activated","inputs|outputs");
  addOption("regularity_check",         OT_BOOLEAN,             true,          "Throw exceptions when NaN or Inf appears during evaluation");
  addOption("gather_stats",             OT_BOOLEAN,             false,         "Flag to indicate wether statistics must be gathered");
  
  verbose_ = false;
  jacgen_ = 0;
  spgen_ = 0;
  user_data_ = 0;
  monitor_inputs_ = false;
  monitor_outputs_ = false;
  
  inputScheme  = SCHEME_unknown;
  outputScheme = SCHEME_unknown;
}



FXInternal::~FXInternal(){
}

void FXInternal::init(){
  verbose_ = getOption("verbose");
  regularity_check_ = getOption("regularity_check");
  bool store_jacobians = getOption("store_jacobians");
  casadi_assert_warning(!store_jacobians,"Option \"store_jacobians\" has been deprecated. Jacobians are now always cached.");
  
  // Allocate data for sensitivities (only the method in this class)
  FXInternal::updateNumSens(false);
  
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
  
  gather_stats_ = getOption("gather_stats");

  // Mark the function as initialized
  is_init_ = true;
}

void FXInternal::updateNumSens(bool recursive){
  // Get the new number
  nfdir_ = getOption("number_of_fwd_dir");
  nadir_ = getOption("number_of_adj_dir");
  
  // Warn if the number exceeds the maximum
  casadi_assert_warning(nfdir_ <= int(getOption("max_number_of_fwd_dir")), "The number of forward directions exceeds the maximum number. Decrease \"number_of_fwd_dir\" or increase \"max_number_of_fwd_dir\"");
  casadi_assert_warning(nadir_ <= int(getOption("max_number_of_adj_dir")), "The number of adjoint directions exceeds the maximum number. Decrease \"number_of_adj_dir\" or increase \"max_number_of_adj_dir\"");
  
  // Allocate memory for the seeds and sensitivities
  for(vector<FunctionIO>::iterator it=input_.begin(); it!=input_.end(); ++it){
    it->dataF.resize(nfdir_,it->data);
    it->dataA.resize(nadir_,it->data);
  }
  for(vector<FunctionIO>::iterator it=output_.begin(); it!=output_.end(); ++it){
    it->dataF.resize(nfdir_,it->data);
    it->dataA.resize(nadir_,it->data);
  }
}

void FXInternal::requestNumSens(int nfwd, int nadj){
  // Request the number of directions to the number that we would ideally have
  int nfwd_new = std::max(nfwd,std::max(nfdir_,int(getOption("number_of_fwd_dir"))));
  int nadj_new = std::max(nadj,std::max(nadir_,int(getOption("number_of_adj_dir"))));
  
  // Cap at the maximum
  nfwd_new = std::min(nfwd_new,int(getOption("max_number_of_fwd_dir")));
  nadj_new = std::min(nadj_new,int(getOption("max_number_of_adj_dir")));
  
  // Update the number of directions, if needed
  if(nfwd_new>nfdir_ || nadj_new>nadir_){
    setOption("number_of_fwd_dir",nfwd_new);
    setOption("number_of_adj_dir",nadj_new);
    updateNumSens(true);
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

FX FXInternal::gradient(int iind, int oind){
  // Assert scalar
  casadi_assert_message(output(oind).scalar(),"Only gradients of scalar functions allowed. Use jacobian instead.");
  
  // Generate gradient function
  FX ret = getGradient(iind,oind);
  
  // Give it a suitable name
  stringstream ss;
  ss << "gradient_" << getOption("name") << "_" << iind << "_" << oind;
  ret.setOption("name",ss.str());
  
  return ret;
}
  
FX FXInternal::hessian(int iind, int oind){
  log("FXInternal::hessian");
  
  // Assert scalar
  casadi_assert_message(output(oind).scalar(),"Only hessians of scalar functions allowed.");
  
  // Generate gradient function
  FX ret = getHessian(iind,oind);
  
  // Give it a suitable name
  stringstream ss;
  ss << "hessian_" << getOption("name") << "_" << iind << "_" << oind;
  ret.setOption("name",ss.str());
  
  return ret;
}
  
FX FXInternal::getGradient(int iind, int oind){
  casadi_error("FXInternal::getGradient: getGradient not defined for class " << typeid(*this).name());
}
  
FX FXInternal::getHessian(int iind, int oind){
  log("FXInternal::getHessian");

  // Create gradient function
  log("FXInternal::getHessian generating gradient");
  FX g = gradient(iind,oind);
  g.setOption("numeric_jacobian",getOption("numeric_hessian"));
  g.setOption("verbose",getOption("verbose"));
  g.init();
  
  // Return the Jacobian of the gradient, exploiting symmetry (the gradient has output index 0)
  log("FXInternal::getHessian generating Jacobian of gradient");
  return g.jacobian(iind,0,false,true);
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

void FXInternal::setNumInputs(int num_in){
  input_.resize(num_in);
}

void FXInternal::setNumOutputs(int num_out){
  output_.resize(num_out);  
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
    
    // Override default behavior?
    if(getOption("ad_mode") == "forward"){
      use_fwd = true;
    } else if(getOption("ad_mode") == "reverse"){
      use_fwd = false;
    }
    
    // Reset the virtual machine
    spInit(use_fwd);

    // Clear the forward seeds/adjoint sensitivities
    for(int ind=0; ind<getNumInputs(); ++ind){
      vector<double> &v = inputNoCheck(ind).data();
      if(!v.empty()) fill_n(get_bvec_t(v),v.size(),bvec_t(0));
    }

    // Clear the adjoint seeds/forward sensitivities
    for(int ind=0; ind<getNumOutputs(); ++ind){
      vector<double> &v = outputNoCheck(ind).data();
      if(!v.empty()) fill_n(get_bvec_t(v),v.size(),bvec_t(0));
    }
    
    // Get seeds and sensitivities
    bvec_t* input_v = get_bvec_t(inputNoCheck(iind).data());
    bvec_t* output_v = get_bvec_t(outputNoCheck(oind).data());
    bvec_t* seed_v = use_fwd ? input_v : output_v;
    bvec_t* sens_v = use_fwd ? output_v : input_v;
    
    // Number of sweeps needed
    int nsweep = use_fwd ? nsweep_fwd : nsweep_adj;
    
    // The number of zeros in the seed and sensitivity directions
    int nz_seed = use_fwd ? nz_in  : nz_out;
    int nz_sens = use_fwd ? nz_out : nz_in;

    // Print
    if(verbose()){
      std::cout << "FXInternal::getJacSparsity: using " << (use_fwd ? "forward" : "adjoint") << " mode: ";
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
      for(int i=0; i<bvec_size && offset+i<nz_seed; ++i){
        seed_v[offset+i] = 0;
      }
    }
    
    // Set inputs and outputs to zero
    for(int ind=0; ind<getNumInputs(); ++ind) input(ind).setZero();
    for(int ind=0; ind<getNumOutputs(); ++ind) output(ind).setZero();

    // Construct sparsity pattern
    CRSSparsity ret = sp_triplet(nz_out, nz_in,use_fwd ? jrow : jcol, use_fwd ? jcol : jrow);
    
    // Return sparsity pattern
    if(verbose()){
      std::cout << "Formed Jacobian sparsity pattern (dimension " << ret.shape() << ", " << 100*double(ret.size())/ret.numel() << " \% nonzeros)." << endl;
      std::cout << "FXInternal::getJacSparsity end " << std::endl;
    }
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
  log("FXInternal::getPartition begin");
  
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
    log("FXInternal::getPartition starColoring");
    D1 = A.starColoring();
    if(verbose()){
      cout << "Star coloring completed: " << D1.size1() << " directional derivatives needed (" << A.size2() << " without coloring)." << endl;
    }
    
  } else {
    
    // Test unidirectional coloring using forward mode
    if(test_ad_fwd){
      log("FXInternal::getPartition unidirectional coloring (forward mode)");
      D1 = AT.unidirectionalColoring(A);
      if(verbose()){
        cout << "Forward mode coloring completed: " << D1.size1() << " directional derivatives needed (" << A.size2() << " without coloring)." << endl;
      }
    }
      
    // Test unidirectional coloring using reverse mode
    if(test_ad_adj){
      log("FXInternal::getPartition unidirectional coloring (adjoint mode)");
      D2 = A.unidirectionalColoring(AT);
      if(verbose()){
        cout << "Adjoint mode coloring completed: " << D2.size1() << " directional derivatives needed (" << A.size1() << " without coloring)." << endl;
      }
    }
    
    // Adjoint mode penalty factor (adjoint mode is usually more expensive to calculate)
    int adj_penalty = 2;

    // Use whatever required less colors if we tried both (with preference to forward mode)
    if(test_ad_fwd && test_ad_adj){
      if((D1.size1() <= adj_penalty*D2.size1())){
        D2=CRSSparsity();
        log("Forward mode chosen");
      } else {
        D1=CRSSparsity();
        log("Adjoint mode chosen");
      }
    }
    log("FXInternal::getPartition end");
  }
}

void FXInternal::evalSX(const std::vector<SXMatrix>& arg, std::vector<SXMatrix>& res, 
      const std::vector<std::vector<SXMatrix> >& fseed, std::vector<std::vector<SXMatrix> >& fsens, 
      const std::vector<std::vector<SXMatrix> >& aseed, std::vector<std::vector<SXMatrix> >& asens,
      bool output_given){
  casadi_error("FXInternal::evalSX not defined for class " << typeid(*this).name());
}

void FXInternal::evalMX(const std::vector<MX>& arg, std::vector<MX>& res, 
      const std::vector<std::vector<MX> >& fseed, std::vector<std::vector<MX> >& fsens, 
      const std::vector<std::vector<MX> >& aseed, std::vector<std::vector<MX> >& asens,
      bool output_given){
  casadi_error("FXInternal::evalMX not defined for class " << typeid(*this).name());
}

void FXInternal::spEvaluate(bool fwd){
  // By default, everything is assumed to depend on everything
  
  // Variable which depends on all everything
  bvec_t all_depend(0);
  if(fwd){
    // Get dependency on all inputs
    for(int iind=0; iind<getNumInputs(); ++iind){
      const DMatrix& m = inputNoCheck(iind);
      const bvec_t* v = reinterpret_cast<const bvec_t*>(m.ptr());
      for(int i=0; i<m.size(); ++i){
        all_depend |= v[i];
      }
    }
    
    // Propagate to all outputs
    for(int oind=0; oind<getNumOutputs(); ++oind){
      DMatrix& m = outputNoCheck(oind);
      bvec_t* v = reinterpret_cast<bvec_t*>(m.ptr());
      for(int i=0; i<m.size(); ++i){
        v[i] = all_depend;
      }
    }
    
  } else {
    
    // Get dependency on all outputs
    for(int oind=0; oind<getNumOutputs(); ++oind){
      const DMatrix& m = outputNoCheck(oind);
      const bvec_t* v = reinterpret_cast<const bvec_t*>(m.ptr());
      for(int i=0; i<m.size(); ++i){
        all_depend |= v[i];
      }
    }
    
    // Propagate to all inputs
    for(int iind=0; iind<getNumInputs(); ++iind){
      DMatrix& m = inputNoCheck(iind);
      bvec_t* v = reinterpret_cast<bvec_t*>(m.ptr());
      for(int i=0; i<m.size(); ++i){
        v[i] |= all_depend;
      }
    }
  }
}

FX FXInternal::jacobian(int iind, int oind, bool compact, bool symmetric){
  // Return value
  FX ret;
  
  // Generate Jacobian
  if(jacgen_!=0){
    // Use user-provided routine to calculate Jacobian
    FX fcn = shared_from_this<FX>();
    ret = jacgen_(fcn,iind,oind,user_data_);
  } else if(bool(getOption("numeric_jacobian"))){
    ret = getNumericJacobian(iind,oind,compact,symmetric);
  } else {
    // Use internal routine to calculate Jacobian
    ret = getJacobian(iind,oind,compact, symmetric);
  }
  
  // Give it a suitable name
  stringstream ss;
  ss << "jacobian_" << getOption("name") << "_" << iind << "_" << oind;
  ret.setOption("name",ss.str());
  ret.setOption("verbose",getOption("verbose"));
  return ret;
}

FX FXInternal::getJacobian(int iind, int oind, bool compact, bool symmetric){
  return getNumericJacobian(iind,oind,compact,symmetric);
}

FX FXInternal::derivative(int nfwd, int nadj){
  // Quick return if 0x0
  if(nfwd==0 && nadj==0) return shared_from_this<FX>();

  // Check if there are enough forward directions allocated
  if(nfwd>=derivative_fcn_.size()){
    derivative_fcn_.resize(nfwd+1);
  }

  // Check if there are enough adjoint directions allocated
  if(nadj>=derivative_fcn_[nfwd].size()){
    derivative_fcn_[nfwd].resize(nadj+1);
  }

  // Return value
  FX& ret = derivative_fcn_[nfwd][nadj];

  // Check if already cached
  if(ret.isNull()){
    // Generate a new function
    ret = getDerivative(nfwd,nadj);
    
    // Give it a suitable name
    stringstream ss;
    ss << "derivative_" << getOption("name") << "_" << nfwd << "_" << nadj;
    ret.setOption("name",ss.str());
    
    // Initialize it
    ret.init();
  }

  // Return cached or generated function
  return ret;
}

FX FXInternal::getDerivative(int nfwd, int nadj){
  casadi_error("FXInternal::getDerivative not defined for class " << typeid(*this).name());
}

int FXInternal::getNumScalarInputs() const{
  int ret=0;
  for(int iind=0; iind<getNumInputs(); ++iind){
    ret += input(iind).size();
  }
  return ret;
}

int FXInternal::getNumScalarOutputs() const{
  int ret=0;
  for(int oind=0; oind<getNumOutputs(); ++oind){
    ret += output(oind).size();
  }
  return ret;
}

int FXInternal::inputSchemeEntry(const std::string &name) const {
  return schemeEntry(inputScheme,name);
}

int FXInternal::outputSchemeEntry(const std::string &name) const {
  return schemeEntry(outputScheme,name);
}

int FXInternal::schemeEntry(InputOutputScheme scheme, const std::string &name) const {
  if (scheme==SCHEME_unknown) casadi_error("Unable to look up '" <<  name<< "' in input scheme, as the input scheme of this function is unknown. You can only index with integers.");
  if (name=="") casadi_error("FXInternal::inputSchemeEntry: you supplied an empty string as the name of a entry in " << getSchemeName(scheme) << ". Available names are: " << getSchemeEntryNames(scheme) << ".");
  int n = getSchemeEntryEnum(scheme,name);
  if (n==-1) casadi_error("FXInternal::inputSchemeEntry: could not find entry '" << name << "' in " << getSchemeName(scheme) << ". Available names are: " << getSchemeEntryNames(scheme) << ".");
  return n;
}

void FXInternal::call(const MXVector& arg, MXVector& res,  const MXVectorVector& fseed, MXVectorVector& fsens, 
                      const MXVectorVector& aseed, MXVectorVector& asens, 
                      bool output_given, bool always_inline, bool never_inline){
  casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");
  
  // Lo logic for inlining yet
  bool inline_function = always_inline;
  
  if(inline_function){
    // Evaluate the function symbolically
    evalMX(arg,res,fseed,fsens,aseed,asens,output_given);
    
  } else {
    // Create a call-node
    assertInit();
    
    // Argument checking
    casadi_assert_message(arg.size()<=getNumInputs(), "FX::call: number of passed-in dependencies (" << arg.size() << ") should not exceed the number of inputs of the function (" << getNumInputs() << ").");

    // Assumes initialised
    for(int i=0; i<arg.size(); ++i){
      if(arg[i].isNull() || arg[i].empty() || input(i).isNull() || input(i).empty()) continue;
      casadi_assert_message(arg[i].size1()==input(i).size1() && arg[i].size2()==input(i).size2(),
                            "Evaluation::shapes of passed-in dependencies should match shapes of inputs of function." << 
                            std::endl << "Input argument " << i << " has shape (" << input(i).size1() << 
                            "," << input(i).size2() << ") while a shape (" << arg[i].size1() << "," << arg[i].size2() << 
                            ") was supplied.");
    }
    EvaluationMX::create(shared_from_this<FX>(),arg,res,fseed,fsens,aseed,asens,output_given);
  }
}

void FXInternal::call(const std::vector<SXMatrix>& arg, std::vector<SXMatrix>& res, 
                      const std::vector<std::vector<SXMatrix> >& fseed, std::vector<std::vector<SXMatrix> >& fsens, 
                      const std::vector<std::vector<SXMatrix> >& aseed, std::vector<std::vector<SXMatrix> >& asens,
                      bool output_given, bool always_inline, bool never_inline){
  casadi_assert_message(!(always_inline && never_inline), "Inconsistent options");
  casadi_assert_message(!never_inline, "SX expressions do not support call-nodes");
  evalSX(arg,res,fseed,fsens,aseed,asens,output_given);
}

FX FXInternal::getNumericJacobian(int iind, int oind, bool compact, bool symmetric){
  vector<MX> arg = symbolicInput();
  vector<MX> res = shared_from_this<FX>().call(arg);
  FX f = MXFunction(arg,res);
  f.setOption("numeric_jacobian", false);
  f.init();
  return f->getNumericJacobian(iind,oind,compact,symmetric);
}

} // namespace CasADi

