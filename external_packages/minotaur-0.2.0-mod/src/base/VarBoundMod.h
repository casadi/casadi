// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file VarBoundMod.h
 * \brief Declare the class VarBoundMod. It is used to save modifications in a
 * bound of a variable. Also declare VarBoundMod2 that is used to change both
 * bounds of a variable.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURVARBOUNDMOD_H
#define MINOTAURVARBOUNDMOD_H

#include "Modification.h"

namespace Minotaur {
  class Engine;
  class VarBoundMod;
  typedef boost::shared_ptr <Engine> EnginePtr;
  typedef boost::shared_ptr<VarBoundMod> VarBoundModPtr;
  typedef std::vector < VarBoundModPtr > VarBoundModVector;
  typedef VarBoundModVector::iterator VarBoundModIter;
  typedef VarBoundModVector::const_iterator VarBoundModConstIter;


  /// Modification of a single bound on a variable.
  class VarBoundMod : public Modification {
    public:
      /// Construct.
      VarBoundMod(VariablePtr var, BoundType lu, double new_val);

      /// Destroy.
      ~VarBoundMod();

      // base class method.
      ModificationPtr fromRel(RelaxationPtr, ProblemPtr ) const;

      // base class method.
      ModificationPtr toRel(ProblemPtr, RelaxationPtr) const;

      /// Get the variable whose bound is changed.
      VariablePtr getVar() const;

      /// Get the type of bound that is changed: lower or upper.
      BoundType getLU() const;

      /// Get new value of the bound.
      double getNewVal() const;

      // Implement Modification::applyToProblem().
      void applyToProblem(ProblemPtr problem);

      // Implement Modification::undoToProblem().
      void undoToProblem(ProblemPtr problem);

      // Implement Modification::write().
      void write(std::ostream &out) const;

    private:
      /// Lower or upper bound.
      BoundType lu_;

      /// The new value of the bound.
      double newVal_;

      /// The old value of the bound.
      double oldVal_;

      /// The variable whose bounds are modified.
      VariablePtr var_;
  };

  class VarBoundMod2;
  typedef boost::shared_ptr<VarBoundMod2> VarBoundMod2Ptr;
  typedef boost::shared_ptr<const VarBoundMod2> ConstVarBoundMod2Ptr;  
  typedef std::vector < VarBoundMod2Ptr > VarBoundMod2Vector;
  typedef VarBoundMod2Vector::iterator VarBoundMod2Iter;
  typedef VarBoundMod2Vector::const_iterator VarBoundMod2ConstIter;

  /// Modification of a both bounds on a variable.
  class VarBoundMod2 : public Modification {
    public:
      /// Construct.
      VarBoundMod2(VariablePtr var, double new_lb, double new_ub);

      /// Destroy.
      ~VarBoundMod2();

      // Implement Modification::applyToProblem().
      void applyToProblem(ProblemPtr problem);

      // base class method.
      ModificationPtr fromRel(RelaxationPtr, ProblemPtr ) const;

      /// Get the variable whose bound is changed.
      VariablePtr getVar() const;

      /// Get new value of the bound.
      double getNewLb() const;

      /// Get new value of the bound.
      double getNewUb() const;

      // base class method.
      ModificationPtr toRel(ProblemPtr, RelaxationPtr) const;

      // Implement Modification::undoToProblem().
      void undoToProblem(ProblemPtr problem);

      // Implement Modification::write().
      void write(std::ostream &) const {};

    private:
      /// The new lower bound.
      double newLb_;

      /// The new upper bound.
      double newUb_;

      /// Old lower bound.
      double oldLb_;

      /// Old upper bound.
      double oldUb_;

      /// The variable whose bounds are modified.
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
