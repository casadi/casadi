//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//


/**
 * \file PerspList.h 
 * \brief Declare base class PerspList. 
 * \author Serdar Yildiz, Argonne National Laboratory 
*/

#ifndef MINOTAURPERSPLIST_H
#define MINOTAURPERSPLIST_H

#include <map>
#include <fstream>
using std::ofstream;
#include <string>
using std::string;

#include "Relaxation.h"
#include "Constraint.h"
#include "Variable.h"
#include "Types.h"
#include "Environment.h"

namespace Minotaur {

  struct PerspListStats
  {
    UInt totalcons; // Total number of constraints checked.
    UInt totalpersp; // Total number of perspective constraints obtained.
  };

  class PerspList;
  typedef boost::shared_ptr<PerspList> PerspListPtr;
  typedef boost::shared_ptr<const PerspList> ConstPerspListPtr;
  typedef PerspListStats* PerspListStatsPtr;
  typedef PerspListStats const * ConstPerspListStatsPtr;

  typedef std::map<VariablePtr, std::pair<ConstConstraintPtr, ConstConstraintPtr> > VarUbLb;
  typedef boost::shared_ptr<VarUbLb> VarUbLbPtr;
  typedef std::pair<ConstConstraintPtr, ConstVariablePtr> ConsVar;
  //typedef std::map<ConstConstraintPtr, VarUbLbPtr> PerspCons;
  typedef std::map<ConsVar, VarUbLbPtr> PerspCons;
  typedef boost::shared_ptr<PerspCons> PerspConsPtr;
  typedef boost::shared_ptr<const PerspCons> ConstPerspConsPtr;

  /** 
   * This class identifies the constraints that can be used for 
   * perspective cut generation. 
   */
  class PerspList {
  public:
    /// Default constructor.
    PerspList();

    /// Constructs from a given relaxation.
    PerspList(RelaxationPtr rel, EnvPtr env);

    /// Destructor.
    ~PerspList();

    /// Checks if a constraint is a Perspective constraint.
    bool evalConstraint(ConstConstraintPtr cons, VarUbLbPtr boundcons,
                        VariablePtr& binvar);

    /// Checks if the constraint does not have any multi-variable terms
    /// that includes u and x together, i.e. constraint is separable 
    /// such that f(x) + cu <= 0.
    bool separable(ConstConstraintPtr cons, ConstVariablePtr binvar);

    /// Checks if all the variables are continuous or at most one binary.
    /// Otherwise, we cannot generate perspective cuts.
    bool checkVarTypes(const FunctionPtr f, ConstVariablePtr& binvar);

    /// Checks if the variables are bounded by only one binary variable.
    bool checkVarsBounds(const FunctionPtr f, ConstVariablePtr binvar,
                         VarUbLbPtr boundcons);

    /// Checks if a given variable is bounded by binary variable.
    bool checkVarBounds(ConstVariablePtr var, ConstVariablePtr binvar,
                        VarUbLbPtr varbounds);

    /// Add a constraint to the lists.
    void addConstraint(ConstConstraintPtr cons, VarUbLbPtr boundcons,
                       VariablePtr binvar);

    /// Generate list of perspective constraints.
    void generateList();

    /// Generate map of variables that are in the initial variable's constraint
    /// set.
    bool initialBinary(ConstVariablePtr var, VarSetPtr binaries);

    /// Get total number of perspective constraints.
    UInt getNumPersp() const {return list_->size();};

    /// Get total number of constraints checked.
    UInt getNumConsChecked() const {return stats_->totalcons;};

    /// Get the statistics aboud perspective identification.
    ConstPerspListStatsPtr getStats() const {return stats_;};

    /// Get a pointer to the vector that constains perspective constraints.
    ConstPerspConsPtr getPerspCons() const {return list_;};

    /// Print out the perspective structure.
    void printPersp(ConstConstraintPtr cons, VarUbLbPtr boundcons, 
                    ConstVariablePtr binvar);

  private:
    /// Environment.
    EnvPtr env_;
    /// Relaxation that we identify perspective constraints.
    RelaxationPtr rel_;
    /// Perspective constraint list.
    PerspConsPtr list_;
    /// Statistics about perspective constraints.
    PerspListStatsPtr stats_;
    /// Output file.
    ofstream output_;
    /// Output file name.
    string outfile_;
  }; // end of class PerspList.


} // end of namespace.



#endif // MINOTAURPERSPLIST_H




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
