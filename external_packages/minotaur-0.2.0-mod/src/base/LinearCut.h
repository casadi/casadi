//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
//

/**
 * \file LinearCut.h
 * \brief Define the class of valid inequalities of the form \f$ ax \leq b\f$
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURLINEARCUT_H
#define MINOTAURLINEARCUT_H

#include "Cut.h"


namespace Minotaur {

  class Function;
  typedef boost::shared_ptr<Function> FunctionPtr;
  /**
   * A linear cut is a valid inequality of the form \f$ l \leq ax \leq u \f$,
   * where \f$a \in \mathcal{R}^n, l \in \mathcal{R}, u \in \mathcal{R}\f$.
   */
  class LinearCut : public Cut {

    public:
      /// Default Constructor.
      LinearCut();

      /// Construct a cut of the form: \f$ax \leq b \f$
      LinearCut(LinearFunctionPtr lf, double lb, double ub);

      /// Destroy
      ~LinearCut();

      /// By how much does a given point x violate this cut.
      double getViolation(const double *){ return 0.; };

      /// Get ub of the inequality.
      double getUb(){ return ub_; };

      /// Get lb of the inequality.
      double getLb(){ return lb_; };

      /// Add this cut to problem.
      void applyToProblem(ProblemPtr);

      /// Remove this cut from the problem.
      void undoToProblem(ProblemPtr);

      /// Write this cut to the outstream.
      void write(std::ostream &out) const;


    protected:
      /// Pointer to the constraint that is added to the problem.
      ConstraintPtr cons_;

      /// Function or the lhs.
      FunctionPtr f_;

      /// lb.
      double lb_;

      /// Linear function or the lhs.
      LinearFunctionPtr lf_;

      /// ub.
      double ub_;
  };
  typedef boost::shared_ptr<LinearCut> LinearCutPtr;
  typedef std::vector< LinearCutPtr > LinearCutVector;
  typedef CutVector::iterator LinearCutIterator;
  typedef CutVector::const_iterator LinearCutConstIterator;
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
