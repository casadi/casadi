// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file Problem.h
 * \brief Declare the derived class SecantMod for modifications associated
 * with a secant.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURLINCONSTRMOD_H
#define MINOTAURLINCONSTRMOD_H

#include "LinearFunction.h"
#include "LinConMod.h"
#include "VarBoundMod.h"

namespace Minotaur {

  class SecantMod;
  typedef boost::shared_ptr<SecantMod> SecantModPtr;
  typedef boost::shared_ptr<const SecantMod> ConstSecantModPtr;  

  /// Modification of a single linear constraint and its upper bound. An old
  /// constraint:
  //
  /// nlf + qf + lf <= ub 
  /// 
  /// is changed to:
  ///
  /// nlf + qf + new_lf <= new_ub 
  ///   
  class SecantMod : public Modification {
  public:
    /// Construct.
    SecantMod(ConstraintPtr con, LinearFunctionPtr new_lf, double new_rhs, 
              VariablePtr x, BoundType lu, double new_b, VariablePtr y);

    /// Destroy.
    ~SecantMod();

    // Implement Modification::applyToProblem()
    void applyToProblem(ProblemPtr problem);

    // base class method.
    ModificationPtr toRel(ProblemPtr, RelaxationPtr) const
    {return SecantModPtr();};

    // base class method.
    ModificationPtr fromRel(RelaxationPtr, ProblemPtr) const
    {return SecantModPtr();};

    // Implement Modification::undoToProblem()
    void undoToProblem(ProblemPtr problem);

    /// Get the auxiliary variable or variable y where \f$y <= x^2\f$.
    VariablePtr getY();

    // Implement Modification::write()
    void write(std::ostream &) const {};

  private:
    /// Changes in linear function in con_
    LinConModPtr lmod_;

    /// Bound changes on x.
    VarBoundModPtr xmod_;

    /// Bound changes on y.
    VarBoundMod2Ptr ymod_;
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
