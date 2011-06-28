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

#include "ooqp_internal.hpp"

#include "QpGenData.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "GondzioSolver.h"
#include "QpGenSparseMa27.h"
//#include "QpGenSparseMa57.h"

using namespace std;
namespace CasADi {
namespace Interfaces {

OOQPInternal* OOQPInternal::clone() const{
  // Return a deep copy
  OOQPInternal* node = new OOQPInternal(H,G,A);
  if(!node->is_init)
    node->init();
  return node;
}
  
OOQPInternal::OOQPInternal(const CRSSparsity & H, const CRSSparsity & G, const CRSSparsity & A) : QPSolverInternal(H,G,A){
  std::cout << "Warning: OOQP is highly experimental" << std::endl;
}

OOQPInternal::~OOQPInternal(){ 
}

  
void OOQPInternal::init(){
  
const int nx   = 2;
double    c[]  = { 1.5,  -2 };

double  xupp[] = { 20,   0 };
char   ixupp[] = {  1,   0 };

double  xlow[] = {  0,   0 };
char   ixlow[] = {  1,   1 };

const int nnzQ = 3;
int    irowQ[] = {  0,   1,   1 }; 
int    jcolQ[] = {  0,   0,   1 };
double    dQ[] = {  8,   2,  10 };


int my         = 0;
double * b     = 0;
int nnzA       = 0;
int * irowA    = 0;
int * jcolA    = 0;
double * dA    = 0;

const int mz   = 2;
double clow[]  = { 2,   0 };
char  iclow[]  = { 1,   0 };

double cupp[]  = { 0,   6 };
char  icupp[]  = { 0,   1 };

const int nnzC = 4;
int   irowC[]  = { 0,   0,   1,   1};
int   jcolC[]  = { 0,   1,   0,   1};
double   dC[]  = { 2,   1,  -1,   2};
   int usage_ok = 1, quiet = 0;
    
  QpGenSparseMa27 * qp 
    = new QpGenSparseMa27( nx, my, mz, nnzQ, nnzA, nnzC );
  
  QpGenData      * prob = (QpGenData * ) qp->copyDataFromSparseTriple(
        c,      irowQ,  nnzQ,   jcolQ,  dQ,
        xlow,   ixlow,  xupp,   ixupp,
        irowA,  nnzA,   jcolA,  dA,     b,
        irowC,  nnzC,   jcolC,  dC,
        clow,   iclow,  cupp,   icupp );

  QpGenVars      * vars 
    = (QpGenVars *) qp->makeVariables( prob );
  QpGenResiduals * resid 
    = (QpGenResiduals *) qp->makeResiduals( prob );
  
  GondzioSolver  * s     = new GondzioSolver( qp, prob );
  
  if( !quiet ) s->monitorSelf();
  int ierr = s->solve(prob,vars, resid);
  
  QPSolverInternal::init();

  
}

} // namespace Interfaces
} // namespace CasADi

