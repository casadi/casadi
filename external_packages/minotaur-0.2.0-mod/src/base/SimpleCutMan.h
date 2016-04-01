//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
//

/**
 * \file SimpleCutMan.h
 * \brief Very simple implementation of a cut manager.
 * \author Ashutosh Mahajan, IIT Bombay
 */

#ifndef MINOTAURSIMPLECUTMAN_H
#define MINOTAURSIMPLECUTMAN_H

#include <list>
#include "CutManager.h"
#include "Types.h"


namespace Minotaur {

  /**
   * \brief Derived class for managing cuts. Adds all violated cuts from the
   * storage to the relaxation and never removes any. If a new cut is reported
   * but not violated by the current solution then it is added to storage.
   * This manager does not check for duplicacy or any other numerical problems
   * in cuts.
   */
  class SimpleCutMan : public CutManager {

  public:
    /// Empty constructor.
    SimpleCutMan();

    /**
     * \brief Default constructor.
     *
     * \param [in] env Minotaur Environment pointer.
     * \param [in] p Problem pointer to which cuts will be added or deleted.
     */
    SimpleCutMan(EnvPtr env, ProblemPtr p);

    /// Destroy.
    ~SimpleCutMan();

    // Base class method.
    void addCut(CutPtr c);

    // Base class method.
    ConstraintPtr addCut(ProblemPtr p, FunctionPtr f, double lb,
                         double ub, bool directToRel, bool neverDelete);

    // Base class method.
    void addCuts(CutVectorIter cbeg, CutVectorIter cend);

    // Base class method.
    UInt getNumCuts() const;

    // Base class method.
    UInt getNumEnabledCuts() const;

    // Base class method.
    UInt getNumDisabledCuts() const;

    // Base class method.
    UInt getNumNewCuts() const;

    // Base class method.
    void postSolveUpdate(ConstSolutionPtr sol, EngineStatus eng_status);

    // Base class method.
    void separate(ConstSolutionPtr sol, bool *separated, UInt *added);

    // Base class method.
    void write(std::ostream &out) const;

    // Base class method.
    void writeStats(std::ostream &out) const;

  private:
    /// Cuts already in the relaxation. 
    UInt enCuts_;

    /// Environment.
    EnvPtr env_;

    /// For logging.
    LoggerPtr logger_;

    /// For logging.
    const static std::string me_;

    /**
     * Cut storage for new cuts, i.e. those sent to pool after previous
     * separate or postSolveUpdate().
     */
    CutList newCuts_;

    /// The relaxation problem that cuts are added to and deleted from.
    ProblemPtr p_;

    /// Pool of cuts that were left unviolated. They may be added in the future.
    CutList pool_;

    /// A cut will be added only if the violation exceeds violAbs_.
    double violAbs_;

    /// A cut will be added only if the violation exceeds absolute
    /// value of the the activity times the violRel_.
    double violRel_;

    /// Append the newCuts_ to pool_ and clear newCuts_.
    void mvNewToPool_();
  };

//typedef boost::shared_ptr<CutManager> CutManagerPtr;
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
