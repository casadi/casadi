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

#include "csparse_tools.hpp"

extern "C"{
#include "external_packages/CSparse/Include/cs.h"
}

using namespace std;
namespace CasADi{
  namespace Interfaces{
    
    BLT::BLT(CRSSparsity s, int seed){
      nb = s.dulmageMendelsohn(rowperm,colperm,rowblock,colblock,coarse_rowblock,coarse_colblock,seed);
      
      return;
      
      // Call CSparse to perform the BLT
      cs AT; // NOTE: transpose (since CasADi row major, CSparse column major)
      AT.nzmax = s.size();
      AT.m = s.size2();
      AT.n = s.size1();
      AT.p = &s.rowindRef().front();
      AT.i = &s.colRef().front();
      AT.x = 0;
      AT.nz = -1;
      cout << "calling cs_dmperm " << endl;
      csd *perm = cs_dmperm (&AT, seed);
       
      // Save to BLT structure // NOTE: swapping row<>col due to row/column major
      rowperm = vector<int>(perm->q, perm->q + s.size1());
      colperm = vector<int>(perm->p, perm->p + s.size2());
      nb = perm->nb;
      rowblock = vector<int>(perm->s, perm->s + nb + 1);
      colblock = vector<int>(perm->r, perm->r + nb + 1);
      coarse_rowblock = vector<int>(perm->cc, perm->cc+5);
      coarse_colblock = vector<int>(perm->rr, perm->rr+5);
      
      // Free allocated memory and return
      cs_dfree(perm);
    }

  } // namespace Interfaces
} // namespace CasADi
