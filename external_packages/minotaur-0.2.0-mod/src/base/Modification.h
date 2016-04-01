// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file Problem.h
 * \brief Declare the base class Modification.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURMODIFICATION_H
#define MINOTAURMODIFICATION_H

#include "Types.h"

namespace Minotaur {
  class   Modification;
  class   Problem;
  class   Relaxation;
  typedef boost::shared_ptr <Modification> ModificationPtr;
  typedef boost::shared_ptr <Problem> ProblemPtr;
  typedef boost::shared_ptr <Relaxation> RelaxationPtr;

  /**
   * Modification is a (pure) abstract class for changes that can be done to a 
   * problem   
   */
  class Modification {
    protected:
      /// Default constructor
      Modification() {}

    public:
      /// Default destroy
      virtual ~Modification() {};

      /**
       * \brief Covert a modification for a relaxation  to one for its original
       * problem. 
       * \param[in] rel Relaxation for which this mod is applicable.
       * \param[in] p Problem for which the new mod will be applicable.
       * \returns Modification  applicable to p
       */
      virtual ModificationPtr fromRel(RelaxationPtr rel, ProblemPtr p) const
        = 0;

      /**
       * \brief Covert a modification for a problem  to one for its 
       * relaxation. 
       * \param[in] p Problem for which this mod is applicable.
       * \param[in] rel Relaxation for which the new mod will be applicable.
       * \returns Modification  applicable to rel.
       */
      virtual ModificationPtr toRel(ProblemPtr p, RelaxationPtr rel) const = 0;

      /// Apply it to the problem.
      virtual void applyToProblem(ProblemPtr problem) = 0;

      /// Restore the modification for a problem.
      virtual void undoToProblem(ProblemPtr problem) = 0;

      /// Write it to 'out'.
      virtual void write(std::ostream &out) const = 0;
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
