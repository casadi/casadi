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


#include "code_generator.hpp"
#include "fx_internal.hpp"
#include "../matrix/sparsity_tools.hpp"
#include <iomanip>
#include "symbolic/runtime/runtime_embedded.hpp"

using namespace std;
namespace CasADi{
  
  void CodeGenerator::flush(std::ostream& s) const{
    s << includes_.str();
    s << endl;

    // Space saving macro
    s << "#define d double" << endl << endl;
    
    s << auxiliaries_.str();

    // Print integer constants
    stringstream name;
    for(int i=0; i<integer_constants_.size(); ++i){
      name.str(string());
      name << "s" << i;
      printVector(s,name.str(),integer_constants_[i]);
    }

    // Print double constants
    for(int i=0; i<double_constants_.size(); ++i){
      name.str(string());
      name << "c" << i;
      printVector(s,name.str(),double_constants_[i]);
    }

    s << dependencies_.str();
    s << function_.str();
    s << finalization_.str();
  }
  
  std::string CodeGenerator::numToString(int n){
    stringstream ss;
    ss << n;
    return ss.str();
  }
  
  int CodeGenerator::addDependency(const FX& f){
    // Get the current number of functions before looking for it
    size_t num_f_before = added_dependencies_.size();
    
    // Get index of the pattern
    const void* h = static_cast<const void*>(f.get());
    int& ind = added_dependencies_[h];
    
    // Generate it if it does not exist
    if(added_dependencies_.size() > num_f_before){
      // Add at the end
      ind = num_f_before;
      
      // Give it a name
      stringstream name;
      name << "f" << ind;
      
      // Print to file
      f->generateFunction(dependencies_, name.str(), "const d*","d*","d",*this);
    }
    
    return ind;
  }

  void CodeGenerator::printVector(std::ostream &s, const std::string& name, const vector<int>& v){
    s << "int " << name << "[] = {";
    for(int i=0; i<v.size(); ++i){
      if(i!=0) s << ",";
      s << v[i];
    }
    s << "};" << endl;
  }

  void CodeGenerator::printVector(std::ostream &s, const std::string& name, const vector<double>& v){
    s << "d " << name << "[] = {";
    for(int i=0; i<v.size(); ++i){
      if(i!=0) s << ",";
      printConstant(s,v[i]);
    }
    s << "};" << endl;
  }

  void CodeGenerator::addInclude(const std::string& new_include, bool relative_path){
    // Register the new element
    bool added = added_includes_.insert(new_include).second;

    // Quick return if it already exists
    if(!added) return;

    // Print to the header section
    if(relative_path){
      includes_ << "#include \"" << new_include << "\"" << endl;
    } else {
      includes_ << "#include <" << new_include << ">" << endl;
    }
  }
  
  int CodeGenerator::addSparsity(const CRSSparsity& sp){
    // Get the current number of patterns before looking for it
    size_t num_patterns_before = added_sparsities_.size();

    // Get index of the pattern
    const void* h = static_cast<const void*>(sp.get());
    int& ind = added_sparsities_[h];

    // Generate it if it does not exist
    if(added_sparsities_.size() > num_patterns_before){

      // Compact version of the sparsity pattern
      std::vector<int> sp_compact = sp_compress(sp);

      // Codegen vector
      ind = getConstant(sp_compact,true);
    }
    
    return ind;
  }

  int CodeGenerator::getSparsity(const CRSSparsity& sp) const{
    const void* h = static_cast<const void*>(sp.get());
    PointerMap::const_iterator it=added_sparsities_.find(h);
    casadi_assert(it!=added_sparsities_.end());
    return it->second;  
  }

  size_t CodeGenerator::hash(const std::vector<double>& v){
    // Calculate a hash value for the vector
    std::size_t seed=0;
    if(!v.empty()){
      casadi_assert(sizeof(double) % sizeof(size_t)==0);
      const int int_len = v.size()*(sizeof(double)/sizeof(size_t));
      const size_t* int_v = reinterpret_cast<const size_t*>(&v.front());
      for(size_t i=0; i<int_len; ++i){
        hash_combine(seed,int_v[i]);
      }
    }
    return seed;
  }  

  size_t CodeGenerator::hash(const std::vector<int>& v){
    size_t seed=0;
    hash_combine(seed,v);
    return seed;
  }

  int CodeGenerator::getConstant(const std::vector<double>& v, bool allow_adding){
    // Hash the vector
    size_t h = hash(v);
    
    // Try to locate it in already added constants
    pair<multimap<size_t,size_t>::iterator,multimap<size_t,size_t>::iterator> eq = added_double_constants_.equal_range(h);
    for(multimap<size_t,size_t>::iterator i=eq.first; i!=eq.second; ++i){
      if(equal(v,double_constants_[i->second])) return i->second;
    }
    
    if(allow_adding){
      // Add to constants
      int ind = double_constants_.size();
      double_constants_.push_back(v);
      added_double_constants_.insert(pair<size_t,size_t>(h,ind));
      return ind;
    } else {
      casadi_error("Constant not found");
      return -1;
    }
  }

  int CodeGenerator::getConstant(const std::vector<int>& v, bool allow_adding){
    // Hash the vector
    size_t h = hash(v);
    
    // Try to locate it in already added constants
    pair<multimap<size_t,size_t>::iterator,multimap<size_t,size_t>::iterator> eq = added_integer_constants_.equal_range(h);
    for(multimap<size_t,size_t>::iterator i=eq.first; i!=eq.second; ++i){
      if(equal(v,integer_constants_[i->second])) return i->second;
    }
    
    if(allow_adding){
      // Add to constants
      int ind = integer_constants_.size();
      integer_constants_.push_back(v);
      added_integer_constants_.insert(pair<size_t,size_t>(h,ind));
      return ind;
    } else {
      casadi_error("Constant not found");
      return -1;
    }
  }

  int CodeGenerator::getDependency(const FX& f) const{
    const void* h = static_cast<const void*>(f.get());
    PointerMap::const_iterator it=added_dependencies_.find(h);
    casadi_assert(it!=added_dependencies_.end());
    return it->second;
  }

  void CodeGenerator::addAuxiliary(Auxiliary f){
    // Register the new auxiliary
    bool added = added_auxiliaries_.insert(f).second;
    
    // Quick return if it already exists
    if(!added) return;

    // Add the appropriate function
    switch(f){
    case AUX_COPY: 
      auxiliaries_ << codegen_str_copy << endl; 
      break;
    case AUX_SWAP: 
      auxiliaries_ << codegen_str_swap << endl; 
      break;
    case AUX_SCAL: 
      auxiliaries_ << codegen_str_scal << endl; 
      break;
    case AUX_AXPY: 
      auxiliaries_ << codegen_str_axpy << endl; 
      break;
    case AUX_DOT: 
      auxiliaries_ << codegen_str_dot << endl; 
      break;
    case AUX_ASUM: 
      auxiliaries_ << codegen_str_asum << endl; 
      break;
    case AUX_IAMAX: 
      auxiliaries_ << codegen_str_iamax << endl; 
      break;
    case AUX_NRM2: 
      auxiliaries_ << codegen_str_nrm2 << endl; 
      break;
    case AUX_FILL: 
      auxiliaries_ << codegen_str_fill << endl; 
      break;
    case AUX_MM_NT_SPARSE: 
      auxiliaries_ << codegen_str_mm_nt_sparse << endl; 
      break;
    case AUX_SQ: 
      auxSq(); 
      break;
    case AUX_SIGN:
      auxSign(); 
      break;
    case AUX_COPY_SPARSE:
      addAuxiliary(AUX_COPY);
      auxiliaries_ << codegen_str_copy_sparse << endl; 
      break;
    case AUX_TRANS: 
      auxiliaries_ << codegen_str_trans << endl; 
      break;
    }
  }

  void CodeGenerator::auxSq(){
    auxiliaries_ << "inline d sq(d x){ return x*x;}" << endl << endl;
  }

  void CodeGenerator::auxSign(){
    auxiliaries_ << "inline d sign(d x){ return x<0 ? -1 : x>0 ? 1 : x;}" << endl << endl;
  }

  void CodeGenerator::printConstant(std::ostream& s, double v){
    int v_int(v);
    if(v_int==v){
      // Print integer
      s << v_int << ".0";
    } else {
      // Print real
      s << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1) << v;
    }
  }

  void CodeGenerator::copyVector(std::ostream &s, const std::string& arg, std::size_t n, const std::string& res, const std::string& it, bool only_if_exists) const{
    // Quick return if nothing to do
    if(n==0) return;

    // Indent
    s << "  ";

    // Print condition
    if(only_if_exists){
      s << "if(" << res << "!=0) ";
    }

    if(n==1){
      // Simplify if scalar assignment
      s << "*" << res << "=*" << arg << ";" << endl;
    } else {
      // For loop
      s << "for(" << it << "=0; " << it << "<" << n << "; ++" << it << ") " << res << "[" << it << "]=" << arg << "[" << it << "];" << endl;
    }    
  }

  std::string CodeGenerator::casadi_dot(int n, const std::string& x, int inc_x, const std::string& y, int inc_y){
    addAuxiliary(AUX_DOT);
    stringstream ss;
    ss << "casadi_dot(" << n << "," << x << "," << inc_x << "," << y << "," << inc_y << ")";
    return ss.str();
  }

} // namespace CasADi

