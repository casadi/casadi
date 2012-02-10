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

#include "x_function_internal.hpp"
#include "../matrix/sparsity_tools.hpp"

namespace CasADi{

using namespace std;


XFunctionInternal::XFunctionInternal(){
  addOption("topological_sorting",OT_STRING,"breadth-first","Topological sorting algorithm","depth-first|breadth-first");
}

XFunctionInternal::~XFunctionInternal(){
}

CRSSparsity XFunctionInternal::spDetect(int iind, int oind){
  
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
  bool use_fwd = nsweep_fwd <= nsweep_adj;
  
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
    cout << "XFunctionInternal::spDetect: using " << (use_fwd ? "forward" : "adjoint") << " mode: ";
    cout << nsweep << " sweeps needed for " << nz_seed << " directions" << endl;
  }
  
  // Nonzero offset
  int offset = 0;

  // Progress
  int progress = -10;

  // Temporary vectors
  vector<int> jrow, jcol;
  
  // Loop over the variables, ndir variables at a time
  for(int s=0; s<nsweep; ++s){
    // Print progress
    if(verbose()){
      int progress_new = (s*100)/nsweep;
      // Print when entering a new decade
      if(progress_new / 10 > progress / 10){
        progress = progress_new;
        cout << progress << " %"  << endl;
      }
    }
    
    // Give seeds to a set of directions
    for(int i=0; i<bvec_size && offset+i<nz_seed; ++i){
      spGet(use_fwd,ind_seed,offset+i) |= bvec_t(1)<<i;
    }
    
    // Propagate the dependencies
    spProp(use_fwd);
      
    // Number of local seed directions
    int ndir_local = std::min(bvec_size,nz_seed-offset);
    
    // Loop over the nonzeros of the output
    for(int el=0; el<nz_sens; ++el){

      // Get the sparsity sensitivity
      bvec_t spsens = spGet(!use_fwd,ind_sens,el);

      // Clear the seeds for the next sweep
      if(!use_fwd){
        spGet(true,iind,el) = 0; 
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
        spGet(true,iind,offset+i) = 0;
      }
    }
    
    // Update offset
    offset += bvec_size;
  }

  // Construct sparsity pattern
  CRSSparsity ret = sp_triplet(nz_out, nz_in,use_fwd ? jrow : jcol, use_fwd ? jcol : jrow);
  
  // Return sparsity pattern
  if(verbose()) cout << "XFunctionInternal::spDetect end " << endl;
  return ret;
}

void XFunctionInternal::eval(const std::vector<SXMatrix>& input, std::vector<SXMatrix>& output, 
                             const std::vector<std::vector<SXMatrix> >& fwdSeed, std::vector<std::vector<SXMatrix> >& fwdSens, 
                             const std::vector<std::vector<SXMatrix> >& adjSeed, std::vector<std::vector<SXMatrix> >& adjSens,
                             bool output_given, bool eliminate_constants){ 
  evalSX(input,output,fwdSeed,fwdSens,adjSeed,adjSens,output_given,eliminate_constants);
}

void XFunctionInternal::eval(const std::vector<MX>& input, std::vector<MX>& output, 
                             const std::vector<std::vector<MX> >& fwdSeed, std::vector<std::vector<MX> >& fwdSens, 
                             const std::vector<std::vector<MX> >& adjSeed, std::vector<std::vector<MX> >& adjSens,
                             bool output_given, bool eliminate_constants){ 
  evalMX(input,output,fwdSeed,fwdSens,adjSeed,adjSens,output_given,eliminate_constants);
}


} // namespace CasADi

