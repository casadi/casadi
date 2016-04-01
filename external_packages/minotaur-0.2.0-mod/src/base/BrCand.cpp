// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file BrCand.cpp
 * \author Ashutosh Mahajan, Argonne National Laboratory.
 * \brief Define the classes BranchCand for storing candidates for branching.
 */

#include <iostream>

#include "MinotaurConfig.h"
#include "BrCand.h"
#include "Variable.h"

using namespace Minotaur;


bool Minotaur::CompareScore(BrCandPtr c1, BrCandPtr c2)
{
  if (c1->getScore() > c2->getScore()) {
    return false;
  } else if (c1->getScore() < c2->getScore()) {
    return true;
  } else if (c1->getPCostIndex() < c2->getPCostIndex()) {
    return true;
  } else if (c1->getPCostIndex() > c2->getPCostIndex()) {
    return false;
  } else {
    //std::cout << "problem here: " << c1->getPCostIndex() << " " <<
    // c2->getPCostIndex() << std::endl;
    return false;
  }
}


void BrCand::setDist(double, double )
{
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
