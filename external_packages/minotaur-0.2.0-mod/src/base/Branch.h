// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file Branch.h
 * \author Ashutosh Mahajan, Argonne National Laboratory.
 * \brief Declare the class Branch for describing branches in branch-and-bound.
 */


#ifndef MINOTAURBRANCH_H
#define MINOTAURBRANCH_H

#include "Types.h"

namespace Minotaur {

class   BrCand;
class   Modification;
typedef boost::shared_ptr <BrCand> BrCandPtr;
typedef boost::shared_ptr <Modification> ModificationPtr;

/**
 * \brief Base class for storing branching modifications.
 *
 * A Branch just has a vector of modifications. This vector of
 * modifications can be applied to obtain a child node from the parent's
 * relaxation. For each child node, we must have an associated Branch object.
 * The object can also have other information (estimates on lower bounds of
 * child etc).
 */
class Branch {
public:    
  /// Constructor
  Branch();

  /// Destroy
  ~Branch();

  /**
   * \brief Add a problem modification to the current vector of modifications
   * associated with this branch.
   * \param[in] mod The modification that must be added to the child node.
   */
  void addPMod(ModificationPtr mod);

  /**
   * \brief Add a relaxation modification to the current vector of modifications
   * associated with this branch.
   * \param[in] mod The modification that must be added to the child node.
   */
  void addRMod(ModificationPtr mod);

  /** 
   * \brief Set the candidate that was used to generate this branch.
   * \param[in] cand The branching candidate that was used to create this
   * branch.
   */
  void setBrCand(BrCandPtr cand) {brCand_ = cand;};

  /// The first modification in the vector of modifications (problem).
  ModificationConstIterator pModsBegin() const
  { return pMods_.begin(); }

  /// The first modification in the vector of modifications (relaxation).
  ModificationConstIterator rModsBegin() const
  { return rMods_.begin(); }

  /// The last iterator of the vector of modifications (problem).
  ModificationConstIterator pModsEnd() const
  { return pMods_.end(); }

  /// The last iterator of the vector of modifications (relaxation).
  ModificationConstIterator rModsEnd() const
  { return rMods_.end(); }

  /**
   * \brief The reverse iterators are used for undoing the changes to the
   * problem. It is important that the changes are reverted in the reverse
   * order.
   */
  ModificationRConstIterator pModsRBegin() const 
  { return pMods_.rbegin(); }

  /**
   * \brief The reverse iterators are used for undoing the changes to the
   * relaxation. It is important that the changes are reverted in the reverse
   * order.
   */
  ModificationRConstIterator rModsRBegin() const 
  { return rMods_.rbegin(); }

  /**
   * \brief The last reverse iterator for modifications for the problem.
   * Corresponds to the first modification in the vector.
   */
  ModificationRConstIterator pModsREnd() const { return pMods_.rend(); }

  /**
   * \brief The last reverse iterator for modifications for the relaxation.
   * Corresponds to the first modification in the vector.
   */
  ModificationRConstIterator rModsREnd() const { return rMods_.rend(); }

  /**
   * \brief Return the activity or the value of the branching expression
   * before we branched. Used for updating pseudo-costs.
   */
  double getActivity() const;

  /**
   * \brief Set the activity or the value of the branching expression
   * before we branched.
   * \param[in] value The value of activity.
   */
  void setActivity(double value);

  /// Return the branching candidate that was used to create this branch.
  BrCandPtr getBrCand() {return brCand_;};

  /// Write the branch to 'out'
  void write(std::ostream &out) const;

protected:
  /// Name.
  const static std::string me_; 

  /**
   * \brief A vector of modifications of Problem that define this branch.
   *
   * A branch may have more than one modifications. For instance, fixing a
   * variable to zero may have implications on bounds of other variables
   * as well. 
   */
  ModVector pMods_;

  /**
   * \brief A vector of modifications of Relaxation that define this branch.
   *
   * A branch may have more than one modifications. For instance, fixing a
   * variable to zero may have implications on bounds of other variables
   * as well. 
   */
  ModVector rMods_;

  /**
   * \brief The value of the branching expression before we branched.
   *
   * If an integer variable has value 1.3 before branching, and we
   * branch on it, the value is set at 1.3
   */
  double activity_;

  /// Branching candidate that is used to create this branch. 
  BrCandPtr brCand_;
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
