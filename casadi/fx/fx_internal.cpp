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
#include "../sx/evaluation_sx.hpp"
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
  monitor_inputs_ = false;
  monitor_outputs_ = false;
  
  inputScheme  = SCHEME_unknown;
  outputScheme = SCHEME_unknown;
}



FXInternal::~FXInternal(){
}

void FXInternal::init(){
  verbose_ = getOption("verbose");
  bool store_jacobians = getOption("store_jacobians");
  casadi_assert_warning(!store_jacobians,"Option \"store_jacobians\" has been deprecated. Jacobians are now always cached.");
  
  numeric_jacobian_ = getOption("numeric_jacobian");
  bool jac_for_sens = getOption("jac_for_sens");
  casadi_assert_warning(jac_for_sens==false,"The option \"jac_for_sens\" has been deprecated. Ignored.");
  
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
    Jacobian ret(fcn,jblocks.front().second,jblocks.front().first);
    ret.setOption("verbose",getOption("verbose"));
    return ret;
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
      vector<double> &v = input<false>(ind).data();
      if(!v.empty()) fill_n(get_bvec_t(v),v.size(),bvec_t(0));
    }

    // Clear the adjoint seeds/forward sensitivities
    for(int ind=0; ind<getNumOutputs(); ++ind){
      vector<double> &v = output<false>(ind).data();
      if(!v.empty()) fill_n(get_bvec_t(v),v.size(),bvec_t(0));
    }
    
    // Get seeds and sensitivities
    bvec_t* input_v = get_bvec_t(input<false>(iind).data());
    bvec_t* output_v = get_bvec_t(output<false>(oind).data());
    bvec_t* seed_v = use_fwd ? input_v : output_v;
    bvec_t* sens_v = use_fwd ? output_v : input_v;
    
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
      std::cout << "Formed Jacobian sparsity pattern (" << 100*double(ret.size())/ret.numel() << " \% nonzeros)." << endl;
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
  log("MXFunctionInternal::getPartition begin");
  
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

void FXInternal::evalSX(const std::vector<SXMatrix>& arg, std::vector<SXMatrix>& res, 
      const std::vector<std::vector<SXMatrix> >& fseed, std::vector<std::vector<SXMatrix> >& fsens, 
      const std::vector<std::vector<SXMatrix> >& aseed, std::vector<std::vector<SXMatrix> >& asens,
      bool output_given, int offset_begin, int offset_end){
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
      const DMatrix& m = input<false>(iind);
      const bvec_t* v = reinterpret_cast<const bvec_t*>(m.ptr());
      for(int i=0; i<m.size(); ++i){
        all_depend |= v[i];
      }
    }
    
    // Propagate to all outputs
    for(int oind=0; oind<getNumOutputs(); ++oind){
      DMatrix& m = output<false>(oind);
      bvec_t* v = reinterpret_cast<bvec_t*>(m.ptr());
      for(int i=0; i<m.size(); ++i){
        v[i] = all_depend;
      }
    }
    
  } else {
    
    // Get dependency on all outputs
    for(int oind=0; oind<getNumOutputs(); ++oind){
      const DMatrix& m = output<false>(oind);
      const bvec_t* v = reinterpret_cast<const bvec_t*>(m.ptr());
      for(int i=0; i<m.size(); ++i){
        all_depend |= v[i];
      }
    }
    
    // Propagate to all inputs
    for(int iind=0; iind<getNumInputs(); ++iind){
      DMatrix& m = input<false>(iind);
      bvec_t* v = reinterpret_cast<bvec_t*>(m.ptr());
      for(int i=0; i<m.size(); ++i){
        v[i] |= all_depend;
      }
    }
  }
}

FX FXInternal::jacobian_new(int iind, int oind){
  casadi_assert(iind>=0 && iind<getNumInputs());
  casadi_assert(oind>=0 && oind<getNumOutputs());
  
  // Make sure that enough cache entries have been allocated
  if(jacobian_fcn_.size()!=getNumInputs() || jacobian_fcn_[iind].size()!=getNumOutputs()){
    jacobian_fcn_.resize(getNumInputs());
    for(vector<vector<WeakRef> >::iterator it=jacobian_fcn_.begin(); it!=jacobian_fcn_.end(); ++it){
      it->resize(getNumOutputs());
    }
  }
  
  // Weak reference
  WeakRef& ref = jacobian_fcn_[iind][oind];

  // Return value
  FX ret;

  // Check if already cached
  if(ref.isNull()){
    // Generate a new function
    ret = getJacobian(iind,oind);
    
    // Give it a suitable name
    stringstream ss;
    ss << "jacobian[" << iind << "," << oind << "](" << getOption("name") << ")";
    ret.setOption("name",ss.str());
    
    // Initialize it
    ret.init();

    // Cache function for later reference
    ref = ret;
    
  } else {
    // Retrieve cached function
    ret = ref;
  }

  // Return cached or generated function
  return ret;
}

FX FXInternal::getJacobian(int iind, int oind){
  vector<pair<int,int> > jblocks;
  jblocks.push_back(pair<int,int>(oind,iind));
  return jacobian_switch(jblocks);
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
    ss << "derivative[" << nfwd << "," << nadj << "](" << getOption("name") << ")";
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
  if (scheme=SCHEME_unknown) casadi_error("Unable to look up '" <<  name<< "' in input scheme, as the input scheme of this function is unknown. You can only index with integers.");
  if (name=="") casadi_error("FXInternal::inputSchemeEntry: you supplied an empty string as the name of a entry in " << getSchemeName(scheme) << ". Available names are: " << getSchemeEntryNames(scheme) << ".");
  int n = getSchemeEntryEnum(scheme,name);
  if (n==-1) casadi_error("FXInternal::inputSchemeEntry: could not find entry '" << name << "' in " << getSchemeName(scheme) << ". Available names are: " << getSchemeEntryNames(scheme) << ".");
  return n;
}

} // namespace CasADi

