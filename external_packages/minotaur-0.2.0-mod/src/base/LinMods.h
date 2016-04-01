// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file LinMods.h
 * \brief Declare the class LinMods. It is used to save several modifications 
 * in the linear functions of constraints and changes in bounds of variables.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURLINMODS_H
#define MINOTAURLINMODS_H

#include "LinConMod.h"
#include "VarBoundMod.h"

namespace Minotaur {

  class LinMods;
  typedef boost::shared_ptr< LinMods > LinModsPtr;
  typedef boost::shared_ptr< const LinMods > ConstLinModsPtr;

  /** 
   * LinMods keeps record of changes in linear functions of different
   * constraints and the changes in bounds of different variables.
   */
  class LinMods : public Modification {
    public:
      /// Constructor.
      LinMods();

      /// Destroy.
      ~LinMods();

      /// Apply it to the problem.
      void applyToProblem(ProblemPtr problem);

      // base class method.
      ModificationPtr fromRel(RelaxationPtr, ProblemPtr ) const
      {return LinModsPtr();};

      /// Insert a new VarBoundMod
      void insert(VarBoundModPtr bmod);

      /// Insert a new VarBoundMod2
      void insert(VarBoundMod2Ptr b2mod);

      /// Insert a new LinConMod
      void insert(LinConModPtr lmod);

      /**
       * Return true if there are no modifications stored in this object.
       * False otherwise.
       */
      bool isEmpty() const;

      // base class method.
      ModificationPtr toRel(ProblemPtr, RelaxationPtr) const
      {return LinModsPtr();};

      /// Restore the modification for a problem.
      void undoToProblem(ProblemPtr problem);

      /// Write it to 'out'.
      void write(std::ostream &) const;

    private:
      /// Vector of simple bound changes.
      std::vector <VarBoundModPtr> bmods_;

      /// Vector of bound changes where both bounds of a variable are changed.
      std::vector <VarBoundMod2Ptr> bmods2_;

      /// Vector of changes in linear functions of constraints.
      std::vector <LinConModPtr> lmods_;
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
