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

#include "simple_homotopy_nlp_internal.hpp"
#include "casadi/symbolic/std_vector_tools.hpp"
#include "casadi/symbolic/matrix/sparsity_tools.hpp"
#include "casadi/symbolic/matrix/matrix_tools.hpp"
#include "casadi/symbolic/mx/mx_tools.hpp"
#include "casadi/symbolic/function/mx_function.hpp"
#include "casadi/symbolic/casadi_calculus.hpp"
#include <ctime>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cfloat>

using namespace std;
namespace casadi{

  SimpleHomotopyNLPInternal::SimpleHomotopyNLPInternal(const Function& nlp)
      : HomotopyNLPInternal(nlp)
  {
    addOption("nlp_solver",         OT_NLPSOLVER,   GenericType(),
              "The NLP solver to be used by the Homotopy solver");
    addOption("nlp_solver_options", OT_DICTIONARY, GenericType(),
              "Options to be passed to the Homotopy solver");

    //addOption("max_step",            OT_REAL,      1.,
    //          "The maximal homotopy step that is allowed");
    addOption("num_steps",            OT_INTEGER,      10,
              "Take this many steps to go from tau=0 to tau=1.");

  }


  SimpleHomotopyNLPInternal::~SimpleHomotopyNLPInternal(){
  }

  void SimpleHomotopyNLPInternal::init(){
    casadi_warning("SimpleHomotopyNLPSolver is experimental");
    // Call the init method of the base class
    HomotopyNLPInternal::init();

    // Bundle p and tau together;
    MX P = MX::sym("p",hnlp_.input(HNL_P).size()+1);
    std::vector<int> split;
    split.push_back(0);
    split.push_back(np_);
    split.push_back(np_+1);
    std::vector<MX> p_tau = vertsplit(P,split);
    MX x = MX::sym("x",hnlp_.input(HNL_X).size());

    std::vector<MX> nlp_in = nlpIn("x",x,"p",P);

    MXFunction nlp(
                    nlpIn("x",x,"p",P),
                    hnlp_.call(hnlpIn("x",x,"p",p_tau[0],"tau",p_tau[1]))
                  );

    // Create an nlpsolver instance
    NLPSolverCreator nlpsolver_creator = getOption("nlp_solver");
    nlpsolver_ = nlpsolver_creator(nlp);

    if(hasSetOption("nlp_solver_options")){
      nlpsolver_.setOption(getOption("nlp_solver_options"));
    }

    // Initialize the NLP solver
    nlpsolver_.init();

    num_steps_ = getOption("num_steps");

  }

  void SimpleHomotopyNLPInternal::evaluate(){
    nlpsolver_.setInput(input(NLP_SOLVER_X0),NLP_SOLVER_X0);
    nlpsolver_.setInput(input(NLP_SOLVER_LAM_X0),NLP_SOLVER_LAM_X0);
    nlpsolver_.setInput(input(NLP_SOLVER_LBX),NLP_SOLVER_LBX);
    nlpsolver_.setInput(input(NLP_SOLVER_UBX),NLP_SOLVER_UBX);
    nlpsolver_.setInput(input(NLP_SOLVER_LBG),NLP_SOLVER_LBG);
    nlpsolver_.setInput(input(NLP_SOLVER_UBG),NLP_SOLVER_UBG);
    nlpsolver_.setInput(input(NLP_SOLVER_LAM_G0),NLP_SOLVER_LAM_G0);
    //nlpsolver_.setInput(input(NLP_SOLVER_LAM_P0),NLP_SOLVER_LAM_P0);

    std::copy(input(NLP_SOLVER_P).begin(),
              input(NLP_SOLVER_P).end(),
              nlpsolver_.input(NLP_SOLVER_P).begin());

    for (int i=0;i<num_steps_;++i) {
      nlpsolver_.input(NLP_SOLVER_P)[np_] = i*1.0/(num_steps_-1);
      nlpsolver_.evaluate();

      nlpsolver_.setInput(nlpsolver_.output(NLP_SOLVER_LAM_G),NLP_SOLVER_LAM_G0);
      nlpsolver_.setInput(nlpsolver_.output(NLP_SOLVER_LAM_X),NLP_SOLVER_LAM_X0);
     // nlpsolver_.setInput(nlpsolver_.output(NLP_SOLVER_LAM_P),NLP_SOLVER_LAM_P0);

      nlpsolver_.setInput(nlpsolver_.output(NLP_SOLVER_X),NLP_SOLVER_X0);

    }

    setOutput(nlpsolver_.output(NLP_SOLVER_X),NLP_SOLVER_X);
    setOutput(nlpsolver_.output(NLP_SOLVER_LAM_X),NLP_SOLVER_LAM_X);
    setOutput(nlpsolver_.output(NLP_SOLVER_LAM_G),NLP_SOLVER_LAM_G);


    std::copy(nlpsolver_.output(NLP_SOLVER_P).begin(),
              nlpsolver_.output(NLP_SOLVER_P).begin()+np_,
              output(NLP_SOLVER_LAM_P).begin());


  }


} // namespace casadi
