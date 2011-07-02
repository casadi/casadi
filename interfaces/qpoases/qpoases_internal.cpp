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

#include "qpoases_internal.hpp"

#include "casadi/mx/mx_tools.hpp"
#include "casadi/fx/mx_function.hpp"
#include "casadi/mx/densification.hpp"

using namespace std;
namespace CasADi {
namespace Interfaces {

QPOasesInternal* QPOasesInternal::clone() const{
  // Return a deep copy
  QPOasesInternal* node = new QPOasesInternal(H,G,A);
  if(!node->is_init)
    node->init();
  return node;
}
  
QPOasesInternal::QPOasesInternal(const CRSSparsity & H_, const CRSSparsity & G_, const CRSSparsity & A_) : QPSolverInternal(H_,G_,A_){
  addOption("nWSR",OT_INTEGER,-1,"The maximum number of working set recalculations to be performed during the initial homotopy. Default (-1) leads to 5(nx + nc)");
  addOption("hotstart",OT_BOOLEAN,false,"Set to true if you need to solve several QP's in a row that are similar");
  addOption("CPUtime",OT_REAL,0,"The maximum allowed CPU time in seconds for the whole initialisation (and the actually required one on output). Set to zero to disable.");
  std::cout << "Warning: QPOases interface is highly experimental" << std::endl;

}

QPOasesInternal::~QPOasesInternal(){ 
}

void QPOasesInternal::evaluate(int nfdir, int nadir) {
  if (nfdir!=0 || nadir!=0) throw CasadiException("QPOasesInternal::evaluate() not implemented for forward or backward mode");
  
  H_.input().set(input(QP_H));
  H_.evaluate();
  
  A_.input().set(input(QP_A));
  A_.evaluate();
  
  
  CPUtime = getOption("CPUtime").toDouble();
  double* CPUtime_pt=&CPUtime;
  
  if (CPUtime==0) CPUtime_pt = 0;

  int flag = QP.init(&H_.output().data()[0], &input(QP_G).data()[0], &A_.output().data()[0], &input(QP_LBX).data()[0], &input(QP_UBX).data()[0], &input(QP_LBA).data()[0], &input(QP_UBA).data()[0], nWSR,CPUtime_pt);
  
  if(flag!=SUCCESSFUL_RETURN) qpoases_error("Solve",flag);
  
    
  QP.getPrimalSolution( &output(QP_X_OPT).data()[0]);
  output(QP_COST).set(QP.getObjVal());

}

void QPOasesInternal::init(){
  QPSolverInternal::init();

  // Provide densifier functions for H and A
  
  MX Hx = MX("H",H);
  
  MX Hxd = MX::create(new Densification(Hx));
  
  H_ = MXFunction(Hx,Hxd);
  H_.init();
  
  MX Ax = MX("A",A);
  
  MX Axd = MX::create(new Densification(Ax));
  
  A_ = MXFunction(Ax,Axd);
  A_.init();
  
  nWSR = -1;
  if(hasSetOption("nWSR")) nWSR = getOption("nWSR").toInt();
  
  if (nWSR=-1) nWSR = 5 *(nx + nc);
  
  QP = SQProblem( nx,nc);

}

map<int,string> QPOasesInternal::calc_flagmap(){
  map<int,string> f;

  f[SUCCESSFUL_RETURN] = "SUCCESSFUL_RETURN";
  //f[NOT_FINISHED] = "NOT_FINISHED";
  f[RET_MAX_NWSR_REACHED] = "RET_MAX_NWSR_REACHED";
  f[RET_INIT_FAILED] = "RET INIT FAILED";
  f[RET_HOTSTART_FAILED] = "RET_HOTSTART_FAILED";
  return f;
}
  
map<int,string> QPOasesInternal::flagmap = QPOasesInternal::calc_flagmap();

void QPOasesInternal::qpoases_error(const string& module, int flag){
  // Find the error
  map<int,string>::const_iterator it = flagmap.find(flag);
  
  stringstream ss;
  if(it == flagmap.end()){
    ss << "Unknown error (" << flag << ") from module \"" << module << "\".";
  } else {
    ss << "Module \"" << module << "\" returned flag \"" << it->second << "\".";
  }
  ss << " Consult qpOASES documentation.";
  throw CasadiException(ss.str());
}

} // namespace Interfaces
} // namespace CasADi

