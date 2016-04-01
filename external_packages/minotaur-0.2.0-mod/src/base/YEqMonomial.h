//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file YEqMonomial.h
 * \brief Declare class for storing auxiliary variables for monomial
 * expressions.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURYEQMONOMIAL_H
#define MINOTAURYEQMONOMIAL_H

#include "Types.h"

namespace Minotaur {

class MultilinearTermsHandler;
class MonomialFunction;
class Variable;
typedef boost::shared_ptr<MonomialFunction> MonomialFunPtr;
typedef boost::shared_ptr<MultilinearTermsHandler> MultilinearTermsHandlerPtr;
typedef boost::shared_ptr<Variable> VariablePtr;


class YEqMonomial
{
public:
  YEqMonomial(UInt n);
  VariablePtr findY(MonomialFunPtr mf);
  void insert(VariablePtr auxvar, MonomialFunPtr mf);

private:
  std::vector<MonomialFunPtr> mf_;
  DoubleVector hash_;
  UInt n_;
  DoubleVector rand_;
  VarVector y_;
  double evalHash_(MonomialFunPtr mf_);
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
