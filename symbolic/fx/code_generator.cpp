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
    case AUX_COPY: auxCopy(); break;
    case AUX_SWAP: auxSwap(); break;
    case AUX_SCAL: auxScal(); break;
    case AUX_AXPY: auxAxpy(); break;
    case AUX_DOT: auxDot(); break;
    case AUX_ASUM: auxAsum(); break;
    case AUX_IAMAX: auxIamax(); break;
    case AUX_NRM2: auxNrm2(); break;
    case AUX_FILL: auxFill(); break;
    case AUX_MM_NT_SPARSE: auxMmNtSparse(); break;
    case AUX_SIGN: auxSign(); break;
    case AUX_COPY_SPARSE: auxCopySparse(); break;
    case AUX_TRANS: auxTrans(); break;
    }
  }

  void CodeGenerator::auxCopy(){
    stringstream& s = auxiliaries_;

    s << "void casadi_copy(int n, const d* x, int inc_x, d* y, int inc_y){" << endl;
    s << "  int i;" << endl;
    s << "  for(i=0; i<n; ++i){" << endl;
    s << "    *y = *x;" << endl;
    s << "    x += inc_x;" << endl;
    s << "    y += inc_y;" << endl;
    s << "  }" << endl;
    s << "}" << endl;
    s << endl;
  }

  void CodeGenerator::auxSwap(){
    stringstream& s = auxiliaries_;

    s << "void casadi_swap(int n, d* x, int inc_x, d* y, int inc_y){" << endl;
    s << "  d t;" << endl;
    s << "  int i;" << endl;
    s << "  for(i=0; i<n; ++i){" << endl;
    s << "    t = *x;" << endl;
    s << "    *x = *y;" << endl;
    s << "    *y = t;" << endl;
    s << "    x += inc_x;" << endl;
    s << "    y += inc_y;" << endl;
    s << "  }" << endl;
    s << "}" << endl;
    s << endl;
  }

  void CodeGenerator::auxCopySparse(){
    addAuxiliary(AUX_COPY);
    stringstream& s = auxiliaries_;

    s << "void casadi_copy_sparse(const d* x, const int* sp_x, d* y, const int* sp_y){" << endl;
    s << "  int nrow_x = sp_x[0];" << endl;
    s << "  int ncol_x = sp_x[1];" << endl;
    s << "  const int* rowind_x = sp_x+2;" << endl;
    s << "  const int* col_x = sp_x + 2 + nrow_x+1;" << endl;
    s << "  int nnz_x = rowind_x[nrow_x];" << endl;
    s << "  int nrow_y = sp_y[0];" << endl;
    s << "  int ncol_y = sp_y[1];" << endl;
    s << "  const int* rowind_y = sp_y+2;" << endl;
    s << "  const int* col_y = sp_y + 2 + nrow_y+1;" << endl;
    s << "  int nnz_y = rowind_y[nrow_y];" << endl;
    s << "  if(sp_x==sp_y){" << endl;
    s << "    casadi_copy(nnz_x,x,1,y,1);" << endl;
    s << "  } else {" << endl;
    s << "    int i;" << endl;
    s << "    for(i=0; i<nrow_x; ++i){" << endl;
    s << "      int el_x = rowind_x[i];" << endl;
    s << "      int el_x_end = rowind_x[i+1];" << endl;
    s << "      int j_x = el_x<el_x_end ? col_x[el_x] : ncol_x;" << endl;
    s << "      int el_y;" << endl;
    s << "      for(el_y=rowind_y[i]; el_y!=rowind_y[i+1]; ++el_y){" << endl;
    s << "        int j=col_y[el_y];" << endl;
    s << "        while(j_x<j){" << endl;
    s << "          el_x++;" << endl;
    s << "          j_x = el_x<el_x_end ? col_x[el_x] : ncol_x;" << endl;
    s << "        }" << endl;
    s << "        if(j_x==j){" << endl;
    s << "          y[el_y] = x[el_x++];" << endl;
    s << "          j_x = el_x<el_x_end ? col_x[el_x] : ncol_x;" << endl;
    s << "        } else {" << endl;
    s << "          y[el_y] = 0;" << endl;
    s << "        }" << endl;
    s << "      }" << endl;
    s << "    }" << endl;
    s << "  }" << endl;
    s << "}" << endl;
    s << endl;
  }

  void CodeGenerator::auxSign(){
    stringstream& s = auxiliaries_;

    s << "inline d sign(d x){ return x<0 ? -1 : x>0 ? 1 : x;}" << endl;
    s << endl;
  }

  void CodeGenerator::auxScal(){
    stringstream& s = auxiliaries_;
    
    s << "void casadi_scal(int n, d alpha, d* x, int inc_x){" << endl;
    s << "  int i;" << endl;
    s << "  for(i=0; i<n; ++i){" << endl;
    s << "    *x *= alpha;" << endl;
    s << "    x += inc_x;" << endl;
    s << "  }" << endl;
    s << "}" << endl;
    s << endl;
  }

  void CodeGenerator::auxAxpy(){
    stringstream& s = auxiliaries_;

    s << "void casadi_axpy(int n, d alpha, const d* x, int inc_x, d* y, int inc_y){" << endl;
    s << "  int i;" << endl;
    s << "  for(i=0; i<n; ++i){" << endl;
    s << "    *y += alpha**x;" << endl;
    s << "    x += inc_x;" << endl;
    s << "    y += inc_y;" << endl;
    s << "  }" << endl;
    s << "}" << endl;
    s << endl;
  }

  void CodeGenerator::auxDot(){
    stringstream& s = auxiliaries_;

    s << "d casadi_dot(int n, const d* x, int inc_x, d* y, int inc_y){" << endl;
    s << "  d r = 0;" << endl;
    s << "  int i;" << endl;
    s << "  for(i=0; i<n; ++i){" << endl;
    s << "    r += *x**y;" << endl;
    s << "    x += inc_x;" << endl;
    s << "    y += inc_y;" << endl;
    s << "  }" << endl;
    s << "  return r;" << endl;
    s << "}" << endl;
    s << endl;
  }

  void CodeGenerator::auxAsum(){
    stringstream& s = auxiliaries_;

    s << "d casadi_asum(int n, const d* x, int inc_x){" << endl;
    s << "  d r = 0;" << endl;
    s << "  int i;" << endl;
    s << "  for(i=0; i<n; ++i){" << endl;
    s << "    r += fabs(*x);" << endl;
    s << "    x += inc_x;" << endl;
    s << "  }" << endl;
    s << "  return r;" << endl;
    s << "}" << endl;
    s << endl;
  }

  void CodeGenerator::auxIamax(){
    stringstream& s = auxiliaries_;
    
    s << "int casadi_iamax(int n, const d* x, int inc_x){" << endl;
    s << "  d t;" << endl;
    s << "  d largest_value = -1.0;" << endl;
    s << "  int largest_index = -1;" << endl;
    s << "  int i;" << endl;
    s << "  for(i=0; i<n; ++i){" << endl;
    s << "    t = fabs(*x);" << endl;
    s << "    x += inc_x;" << endl;
    s << "    if(t>largest_value){" << endl;
    s << "      largest_value = t;" << endl;
    s << "      largest_index = i;" << endl;
    s << "    }" << endl;
    s << "  }" << endl;
    s << "  return largest_index;" << endl;
    s << "}" << endl;
    s << endl;
  }

  void CodeGenerator::auxFill(){
    stringstream& s = auxiliaries_;

    s << "void casadi_fill(int n, d alpha, d* x, int inc_x){" << endl;
    s << "  int i;" << endl;
    s << "  for(i=0; i<n; ++i){" << endl;
    s << "    *x = alpha;" << endl;
    s << "    x += inc_x;" << endl;
    s << "  }" << endl;
    s << "}" << endl;
    s << endl;
  }

  void CodeGenerator::auxMmNtSparse(){
    stringstream& s = auxiliaries_;

    s << "void casadi_mm_nt_sparse(const d* x, const int* sp_x, const d* trans_y, const int* sp_trans_y, d* z, const int* sp_z){" << std::endl;
    
    s << "  int nrow_x = sp_x[0];" << endl;
    s << "  int ncol_x = sp_x[1];" << endl;
    s << "  const int* rowind_x = sp_x+2;" << endl;
    s << "  const int* col_x = sp_x + 2 + nrow_x+1;" << endl;
    //  s << "  int nnz_x = rowind_x[nrow_x];" << endl;
    
    s << "  int ncol_y = sp_trans_y[0];" << endl;
    s << "  int nrow_y = sp_trans_y[1];" << endl;
    s << "  const int* colind_y = sp_trans_y+2;" << endl;
    s << "  const int* row_y = sp_trans_y + 2 + ncol_y+1;" << endl;
    //  s << "  int nnz_y = colind_y[ncol_y];" << endl;
    
    s << "  int nrow_z = sp_z[0];" << endl;
    s << "  int ncol_z = sp_z[1];" << endl;
    s << "  const int* rowind_z = sp_z+2;" << endl;
    s << "  const int* col_z = sp_z + 2 + nrow_z+1;" << endl;
    //  s << "  int nnz_z = rowind_z[nrow_z];" << endl;
    
    s << "  int i;" << endl;
    s << "  for(i=0; i<nrow_z; ++i){" << endl;
    s << "    int el;" << endl;
    s << "    for(el=rowind_z[i]; el<rowind_z[i+1]; ++el){" << endl;
    s << "      int j = col_z[el];" << endl;
    s << "      int el1 = rowind_x[i];" << endl;
    s << "      int el2 = colind_y[j];" << endl;
    s << "      while(el1 < rowind_x[i+1] && el2 < colind_y[j+1]){ " << endl;
    s << "        int j1 = col_x[el1];" << endl;
    s << "        int i2 = row_y[el2];" << endl;
    s << "        if(j1==i2){" << endl;
    s << "          z[el] += x[el1++] * trans_y[el2++];" << endl;
    s << "        } else if(j1<i2) {" << endl;
    s << "          el1++;" << endl;
    s << "        } else {" << endl;
    s << "          el2++;" << endl;
    s << "        }" << endl;
    s << "      }" << endl;
    s << "    }" << endl;
    s << "  }" << endl;
    s << "}" << endl;
    s << endl;
  }    


  void CodeGenerator::auxNrm2(){
    stringstream& s = auxiliaries_;

    s << "d casadi_nrm2(int n, const d* x, int inc_x){" << endl;
    s << "  d r = 0;" << endl;
    s << "  int i;" << endl;
    s << "  for(i=0; i<n; ++i){" << endl;
    s << "    r += *x**x;" << endl;
    s << "    x += inc_x;" << endl;
    s << "  }" << endl;
    s << "  return sqrt(r);" << endl;
    s << "}" << endl;
    s << endl;
  }

  void CodeGenerator::auxTrans(){
    stringstream& s = auxiliaries_;
    s << "void casadi_trans(const d* x, const int* sp_x, d* y, const int* sp_y, int *tmp){" << endl;
    s << "  int nrow_x = sp_x[0];" << endl;
    s << "  int nnz_x = sp_x[2 + nrow_x];" << endl;
    s << "  const int* col_x = sp_x + 2 + nrow_x+1;" << endl;
    s << "  int nrow_y = sp_y[0];" << endl;
    s << "  const int* rowind_y = sp_y+2;" << endl;
    s << "  int k;" << endl;
    s << "  for(k=0; k<nrow_y; ++k) tmp[k] = rowind_y[k];" << endl;   
    s << "  for(k=0; k<nnz_x; ++k){" << endl;
    s << "    y[tmp[col_x[k]]++] = x[k];" << endl;
    s << "  }" << endl;
    s << "}" << endl;
    s << endl;
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

} // namespace CasADi

