//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//


#ifndef MINOTAURLPRELAXATION_H
#define MINOTAURLPRELAXATION_H

#include "Modification.h"
#include "Relaxation.h"

namespace Minotaur {

  // /**
  // LPRelaxation is a relaxation in which all constraints and objectives
  // have only linear functions.
  // */
  class LPRelaxation : public Relaxation {
    public:
      // /**
      // Construct an empty Nonlinear Relaxation.
      // */
      LPRelaxation();

      // /** 
      // Construct from a problem.
      // */
      LPRelaxation(ProblemPtr problem);

      // /** 
      // Destroy
      // */
      virtual ~LPRelaxation() { };

    protected:
      // /**
      // Pointer to the original problem
      // */
      ConstProblemPtr originalProblem_;
  };

  typedef boost::shared_ptr <LPRelaxation>       LPRelaxationPtr;
  typedef boost::shared_ptr <const LPRelaxation> ConstLPRelaxationPtr;
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
