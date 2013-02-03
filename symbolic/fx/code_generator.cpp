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

using namespace std;
namespace CasADi{
  
  void CodeGenerator::flush(std::ostream& s){
    s << includes_.str();
    s << endl;

    // Space saving macro
    s << "#define d double" << endl << endl;
    
    s << auxiliaries_.str();
    s << sparsities_.str();
    s << dependents_.str();
    s << function_.str();
    s << finalization_.str();
  }
  
void CodeGenerator::generateCopySparse(){
  stringstream& s = auxiliaries_;

  // COPY: y <- x
  s << "inline void casadi_copy(int n, const d* x, int inc_x, d* y, int inc_y){" << endl;
  s << "  int i;" << endl;
  s << "  for(i=0; i<n; ++i){" << endl;
  s << "    *y = *x;" << endl;
  s << "    x += inc_x;" << endl;
  s << "    y += inc_y;" << endl;
  s << "  }" << endl;
  s << "}" << endl;
  s << endl;

  s << "inline void casadi_copy_sparse(const d* x, const int* sp_x, d* y, const int* sp_y){" << endl;
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

std::string CodeGenerator::numToString(int n){
  stringstream ss;
  ss << n;
  return ss.str();
}

int CodeGenerator::findSparsity(const CRSSparsity& sp, const std::map<const void*,int>& sparsity_index){
  const void* h = static_cast<const void*>(sp.get());
  std::map<const void*,int>::const_iterator it=sparsity_index.find(h);
  casadi_assert(it!=sparsity_index.end());
  return it->second;  
}

int CodeGenerator::addDependent(const FX& f){
  // Get the current number of functions before looking for it
  size_t num_f_before = added_dependents_.size();

  // Get index of the pattern
  const void* h = static_cast<const void*>(f.get());
  int& ind = added_dependents_[h];

  // Generate it if it does not exist
  if(added_dependents_.size() > num_f_before){
    // Add at the end
    ind = num_f_before;
          
    // Give it a name
    stringstream name;
    name << "f" << ind;
    
    // Print to file
    f->generateFunction(dependents_, name.str(), "const d*","d*","d",*this);
  }

  return ind;
}

int CodeGenerator::findDependent(const FX& f, const std::map<const void*,int>& dependent_index){
  const void* h = static_cast<const void*>(f.get());
  std::map<const void*,int>::const_iterator it=dependent_index.find(h);
  casadi_assert(it!=dependent_index.end());
  return it->second;
}

void CodeGenerator::printVector(std::ostream &cfile, const std::string& name, const vector<int>& v){
  cfile << "int " << name << "[] = {";
  for(int i=0; i<v.size(); ++i){
    if(i!=0) cfile << ",";
    cfile << v[i];
  }
  cfile << "};" << endl;
}

int CodeGenerator::printSparsity(std::ostream &stream, const CRSSparsity& sp, std::map<const void*,int>& sparsity_index){
  // Get the current number of patterns before looking for it
  size_t num_patterns_before = sparsity_index.size();

  // Get index of the pattern
  const void* h = static_cast<const void*>(sp.get());
  int& ind = sparsity_index[h];

  // Generate it if it does not exist
  if(sparsity_index.size() > num_patterns_before){
    // Add at the end
    ind = num_patterns_before;
    
    // Compact version of the sparsity pattern
    std::vector<int> sp_compact = sp_compress(sp);
      
    // Give it a name
    stringstream name;
    name << "s" << ind;
    
    // Print to file
    printVector(stream,name.str(),sp_compact);
    
    // Separate with an empty line
    stream << endl;
  }

  return ind;
}

  void CodeGenerator::addInclude(const std::string& new_include, bool relative_path){
    // Register the new element in the map
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
    printSparsity(sparsities_,sp,added_sparsities_);
  }

  int CodeGenerator::getSparsity(const CRSSparsity& sp) const{
    findSparsity(sp,added_sparsities_);
  }


} // namespace CasADi

