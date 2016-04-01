// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file BrVarCand.h
 * \author Ashutosh Mahajan, IIT Bombay
 * \brief Declare the classes BrVarCand  storing candidates for branching.
 */

#ifndef MINOTAURBRVARCAND_H
#define MINOTAURBRVARCAND_H

#include <string>
#include "Types.h"
#include "BrCand.h"

namespace Minotaur {

/** 
 * \brief Derived class of BrCand, it defines candidates for branching on
 * bounds of a variable.
 */
class BrVarCand : public BrCand {
public:
  /// Constructor.
  BrVarCand(VariablePtr var, int i, double d, double u);

  // Destroy.
  ~BrVarCand();

  // base class method.
  double getDDist();

  // base class method.
  std::string getName() const;

  // base class method.
  double getUDist();

  /// Get the variable that we are branching on.
  VariablePtr getVar();

  // base class method
  void setDist(double ddist, double udist);

private:
  /// distance from lower-side branch (obtained by changing ub)
  double dDist_;

  /// distance from upper-side branch (obtained by changing lb)
  double uDist_;

  /// The variable on which we should branch.
  VariablePtr var_;

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
