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
  
  // Nonzero offset
  int offset = 0;

  // Progress
  int progress = -10;

  // Temporary vectors
  vector<int> jrow, jcol;
  
  // We choose forward or adjoint based on whichever requires less sweeps
  if(!sp_adj_ok_ || nsweep_fwd <= nsweep_adj){ // forward mode
    if(verbose()) cout << "XFunctionInternal::spDetect: using forward mode: " << nsweep_fwd << " sweeps needed for " << nz_in << " directions" << endl;
    
    // Loop over the variables, ndir variables at a time
    for(int s=0; s<nsweep_fwd; ++s){
      // Print progress
      if(verbose()){
        int progress_new = (s*100)/nsweep_fwd;
        // Print when entering a new decade
        if(progress_new / 10 > progress / 10){
          progress = progress_new;
          cout << progress << " %"  << endl;
        }
      }
      
      // Give seeds to a set of directions
      for(int i=0; i<bvec_size && offset+i<nz_in; ++i){
        spGet(true,iind,offset+i) |= bvec_t(1)<<i;
      }

      // Propagate the dependencies
      spProp(true);
      
      // Number of local seed directions
      int ndir_local = std::min(bvec_size,nz_in-offset);
      
      // Loop over the nonzeros of the output
      for(int el=0; el<nz_out; ++el){

        // Get the sparsity sensitivity
        bvec_t spsens = spGet(false,oind,el);
        
        // If there is a dependency in any of the directions
        if(0 != spsens){
        
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
      for(int i=0; i<bvec_size && offset+i<nz_in; ++i){
        spGet(true,iind,offset+i) = 0;
      }

      // Update offset
      offset += bvec_size;
    }
    
  } else { // Adjoint mode
    if(verbose()) cout << "XFunctionInternal::spDetect: using adjoint mode: " << nsweep_adj << " sweeps needed for " << nz_out << " directions" << endl;
    
    // Loop over the variables, ndir variables at a time
    for(int s=0; s<nsweep_adj; ++s){
      
      // Print progress
      if(verbose()){
        int progress_new = (s*100)/nsweep_adj;
        // Print when entering a new decade
        if(progress_new / 10 > progress / 10){
          progress = progress_new;
          cout << progress << " %"  << endl;
        }
      }
      
      // Give seeds to a set of directions
      for(int i=0; i<bvec_size && offset+i<nz_out; ++i){
        spGet(false,oind,offset+i) |= bvec_t(1)<<i;
      }
      
      // Propagate the dependencies
      spProp(false);

      // Number of local seed directions
      int ndir_local = std::min(bvec_size,nz_out-offset);
      
      // Loop over the nonzeros of the input
      for(int el=0; el<nz_in; ++el){
        // Get the sparsity sensitivity
        bvec_t spsens = spGet(true,iind,el);
        spGet(true,iind,el) = 0; // Clear the seeds for the next sweep

        // If there is a dependency in any of the directions
        if(0 != spsens){
          
          // Loop over seed directions
          for(int i=0; i<ndir_local ; ++i){
              
            // If the output is influenced by the variable
            if((bvec_t(1) << i) & spsens){
            
              // Add to pattern
              jrow.push_back(i+offset);
              jcol.push_back(el);
            }
          }
        }
      }
          
      // Update offset
      offset += bvec_size;
    }
  }

  // Construct sparsity pattern
  CRSSparsity ret = sp_triplet(nz_out,nz_in,jrow,jcol);
  
  // Return sparsity pattern
  if(verbose()) cout << "XFunctionInternal::spDetect end " << endl;
  return ret;
}

} // namespace CasADi

