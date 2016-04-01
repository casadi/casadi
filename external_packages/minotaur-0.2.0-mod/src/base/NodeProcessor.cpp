// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file NodeProcessor.cpp
 * \brief Define two functions of NodeProcessor class. 
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#include "MinotaurConfig.h"
#include "NodeProcessor.h"

using namespace Minotaur;

void NodeProcessor::setBrancher(BrancherPtr brancher)
{
  brancher_ = brancher;
}


void NodeProcessor::processRootNode(NodePtr node, RelaxationPtr rel, 
    SolutionPoolPtr s_pool)
{
  process(node, rel, s_pool);
}


// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
