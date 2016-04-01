// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file LinConMod.h
 * \brief Declare the class LinConMod. It is used to save modifications in the
 * linear functions of constraints.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURLINCONMOD_H
#define MINOTAURLINCONMOD_H

#include "Modification.h"

namespace Minotaur {

class Engine;
class LinConMod;
class LinearFunction;
typedef boost::shared_ptr<Engine> EnginePtr;
typedef boost::shared_ptr< LinConMod > LinConModPtr;
typedef boost::shared_ptr< const LinConMod > ConstLinConModPtr;
typedef std::vector < LinConModPtr > LinConModVector;
typedef LinConModVector::iterator LinConModIter;
typedef LinConModVector::const_iterator LinConModConstIter;
typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;

/** 
 * LinConMod keeps the changes in the linear function part, and the right
 * hand side (only the upper bound) of a constraint.
 */
class LinConMod : public Modification {
public:
  /// Constructor.
  LinConMod(ConstraintPtr con, LinearFunctionPtr new_lf, double new_ub,
            double new_lb);

  /// Destroy.
  ~LinConMod();

  /// Apply it to the problem.
  void applyToProblem(ProblemPtr problem);

  // base class method.
  ModificationPtr fromRel(RelaxationPtr, ProblemPtr ) const
    {return LinConModPtr();};

  // base class method.
  ModificationPtr toRel(ProblemPtr, RelaxationPtr) const
    {return LinConModPtr();};

  /// Restore the modification for a problem.
  void undoToProblem(ProblemPtr problem);

  /// Write it to 'out'.
  void write(std::ostream &out) const;

private:
  /// Pointer to the constraint whose linear part is changed.
  ConstraintPtr con_;

  /// New linear function.
  LinearFunctionPtr newLf_;

  /// New constraint lower bound
  double newLb_;

  /// New constraint upper bound.
  double newUb_;

  /// Old linear function.
  LinearFunctionPtr oldLf_;

  /// Old constraint lower bound
  double oldLb_;

  /// Old constraint upper bound.
  double oldUb_;

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
