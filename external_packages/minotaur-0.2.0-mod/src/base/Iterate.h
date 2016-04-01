//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

#ifndef MINOTAURITERATE_H
#define MINOTAURITERATE_H

// This is an iterate the countains the status, activities, and 
// variable values for the primal variables and multipliers.  What 
// about solvers or
// engines that do not compute multipliers, such as derivative-free
// optimization methods?

#include <vector>

#include "Types.h"

namespace Minotaur {
  class Iterate {
  public:
    Iterate(std::vector<double> &x, std::vector<double> &l);
    ~Iterate() { };

  private:
    std::vector<double> x_;
    std::vector<double> l_;
  };
}

#endif

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
