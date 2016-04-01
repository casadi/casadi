//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file YEqLFs.h
 * \brief Declare the class for storing auxiliary variables for
 * linear functions.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURYEQLFS_H
#define MINOTAURYEQLFS_H

#include "Types.h"

namespace Minotaur {
class LinearFunction;
typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;

class YEqLFs
{
public:
  YEqLFs(UInt n);
  VariablePtr findY(LinearFunctionPtr lf, double k);
  void insert(VariablePtr auxvar, LinearFunctionPtr lf, double k);

private:
  DoubleVector k_;
  std::vector<LinearFunctionPtr> lf_;
  UInt n_;
  DoubleVector rand_;
  DoubleVector hash_;
  VarVector y_;
  double evalHash_(LinearFunctionPtr lf);
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
