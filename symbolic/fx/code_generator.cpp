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
  
  CodeGenerator::CodeGenerator(std::ostream& s) : s_(s){
  }
  
  void CodeGenerator::flush(){
    s_ << includes_.str();
    s_ << endl;

    // Space saving macro
    s_ << "#define d double" << endl << endl;
    
    s_ << auxiliaries_.str();
    s_ << sparsities_.str();
    s_ << dependents_.str();
    s_ << function_.str();
    s_ << finalization_.str();
  }
  
void CodeGenerator::generateCopySparse(std::ostream &stream){
  // COPY: y <- x
  stream << "inline void casadi_copy(int n, const d* x, int inc_x, d* y, int inc_y){" << endl;
  stream << "  int i;" << endl;
  stream << "  for(i=0; i<n; ++i){" << endl;
  stream << "    *y = *x;" << endl;
  stream << "    x += inc_x;" << endl;
  stream << "    y += inc_y;" << endl;
  stream << "  }" << endl;
  stream << "}" << endl;
  stream << endl;

  stream << "inline void casadi_copy_sparse(const d* x, const int* sp_x, d* y, const int* sp_y){" << endl;
  stream << "  int nrow_x = sp_x[0];" << endl;
  stream << "  int ncol_x = sp_x[1];" << endl;
  stream << "  const int* rowind_x = sp_x+2;" << endl;
  stream << "  const int* col_x = sp_x + 2 + nrow_x+1;" << endl;
  stream << "  int nnz_x = rowind_x[nrow_x];" << endl;
  stream << "  int nrow_y = sp_y[0];" << endl;
  stream << "  int ncol_y = sp_y[1];" << endl;
  stream << "  const int* rowind_y = sp_y+2;" << endl;
  stream << "  const int* col_y = sp_y + 2 + nrow_y+1;" << endl;
  stream << "  int nnz_y = rowind_y[nrow_y];" << endl;
  stream << "  if(sp_x==sp_y){" << endl;
  stream << "    casadi_copy(nnz_x,x,1,y,1);" << endl;
  stream << "  } else {" << endl;
  stream << "    int i;" << endl;
  stream << "    for(i=0; i<nrow_x; ++i){" << endl;
  stream << "      int el_x = rowind_x[i];" << endl;
  stream << "      int el_x_end = rowind_x[i+1];" << endl;
  stream << "      int j_x = el_x<el_x_end ? col_x[el_x] : ncol_x;" << endl;
  stream << "      int el_y;" << endl;
  stream << "      for(el_y=rowind_y[i]; el_y!=rowind_y[i+1]; ++el_y){" << endl;
  stream << "        int j=col_y[el_y];" << endl;
  stream << "        while(j_x<j){" << endl;
  stream << "          el_x++;" << endl;
  stream << "          j_x = el_x<el_x_end ? col_x[el_x] : ncol_x;" << endl;
  stream << "        }" << endl;
  stream << "        if(j_x==j){" << endl;
  stream << "          y[el_y] = x[el_x++];" << endl;
  stream << "          j_x = el_x<el_x_end ? col_x[el_x] : ncol_x;" << endl;
  stream << "        } else {" << endl;
  stream << "          y[el_y] = 0;" << endl;
  stream << "        }" << endl;
  stream << "      }" << endl;
  stream << "    }" << endl;
  stream << "  }" << endl;
  stream << "}" << endl;
  stream << endl;
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

int CodeGenerator::printDependent(std::ostream &stream, const FX& f, const std::map<const void*,int>& sparsity_index, std::map<const void*,int>& dependent_index){
  // Get the current number of functions before looking for it
  size_t num_f_before = dependent_index.size();

  // Get index of the pattern
  const void* h = static_cast<const void*>(f.get());
  int& ind = dependent_index[h];

  // Generate it if it does not exist
  if(dependent_index.size() > num_f_before){
    // Add at the end
    ind = num_f_before;
          
    // Give it a name
    stringstream name;
    name << "f" << ind;
    
    // Print to file
    f->generateFunction(stream, name.str(), "const d*","d*","d",sparsity_index,dependent_index);
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



} // namespace CasADi

