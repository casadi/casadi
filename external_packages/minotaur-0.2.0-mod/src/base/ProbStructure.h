//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//


/**
 * \file ProbStructure.h 
 * \brief Declare base class ProbStructure. 
 * \author Serdar Yildiz, Argonne National Laboratory 
*/

#ifndef MINOTAURPROBSTRUCTURE_H
#define MINOTAURPROBSTRUCTURE_H

#include <vector>

#include "Problem.h"
#include "Constraint.h"
#include "Variable.h"
#include "Types.h"
#include "Environment.h"

namespace Minotaur{

struct ProbStructStats
{
  UInt totalcons; // Total number of constraints checked.
  UInt totalGUBs;  // Total number of GUBs obtained.
}; 

class ProbStructure;
typedef boost::shared_ptr<const ProbStructure> ConstProbStructPtr;
typedef boost::shared_ptr<ProbStructure> ProbStructPtr;
typedef ProbStructStats* ProbStructStatsPtr;
typedef ProbStructStats const * ConstProbStructStatsPtr;

typedef std::pair<ConstVariablePtr, ConstConstraintVectorPtr> VarConsPair;
typedef std::map<ConstVariablePtr, ConstConstraintVectorPtr> VarCons;
typedef boost::shared_ptr< VarCons > VarConsPtr;
typedef VarCons::iterator VarConsIterator;
typedef VarCons::const_iterator VarConsConstIterator;


/**
 * Serdar Yildiz: This class identifies GUB constraints in the given problem.
 * One list contains all the GUB constraints just as they found in the
 * constraint list.
 * Second list is a vector of pairs that is in the format of
 * pair<variable,GUBS> pointers to all the GUBs for the corresponding variable.
 */
class ProbStructure{
public:
  /// Default constructor
  ProbStructure();

  /// Constructs from a given problem.
  ProbStructure(ProblemPtr p, EnvPtr env);

  /// Destructor.
  ~ProbStructure();

  /// Checks if a variable is a GUB constraint.
  bool evalConstraint(ConstConstraintPtr cons);
  
  /// Add a constraint to the lists.
  void addConstraint(ConstConstraintPtr cons);
  
  /// Generate the lists for GUBs and GUBs corresponding to variables.
  void generateLists();

  /// Get total number of GUBs.
  UInt getNumGUB() const {return list_->size();};
  
  /// Get total number of GUBs for a given variable.
  UInt getNumVarGUB(ConstVariablePtr var) const;
  
  /// Get total number of constraints checked.
  UInt getNumConsChecked() const {return stats_->totalcons;};
  
  /// Get the statistics about GUB identification.
  ConstProbStructStatsPtr getStats() const {return stats_;};
  
  /// Get a pointer to the vector that contains GUBs for a given variable.
  ConstConstraintVectorPtr getVarGUBs(ConstVariablePtr var) const;

  /// Get a pointer to the vector that contains GUBs.
  ConstConstraintVectorPtr getGUBs() const {return list_;};
  
private:
  // Environment.
  EnvPtr env_;
  // Problem that we identify GUB constraints.
  ProblemPtr p_;
  // Constraint list that contains the pointers to GUB constraints.
  ConstConstraintVectorPtr list_;
  // Vector that contains pair<variable,GUBlist>.
  VarConsPtr varlist_;
  // Statistics about GUBs.
  ProbStructStatsPtr stats_;
};

} // end of namespace

#endif // MINOTAURPROBSTRUCTURE_H


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


