//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file Brancher.h
 * \brief Declare the base class Brancher for finding and creating branches in
 * Branch-and-Bound.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURBRANCHER_H
#define MINOTAURBRANCHER_H

#include "Types.h"

namespace Minotaur {

  class   Node;
  class   Relaxation;
  class   Solution;
  class   SolutionPool;
  typedef boost::shared_ptr <Node> NodePtr;
  typedef boost::shared_ptr <Relaxation> RelaxationPtr;
  typedef boost::shared_ptr <const Solution> ConstSolutionPtr;
  typedef boost::shared_ptr <SolutionPool> SolutionPoolPtr;

  /**
   * \brief A brancher is used to find suitable branches for a given node. e.g.
   * LexicoBrancher. This class is abstract.
   */
  class Brancher {
    
    public:
      /// Default constructor.
      Brancher();

      /// Destroy.
      virtual ~Brancher();

      /**
       * \brief Find a branching candidate. 
       *
       * \return NULL if x does not have any
       * fractional values for integer constrained variables or if no branching
       * candidates are found (e.g. when we realize that problem is infeasible). 
       * \param[in] rel Relaxation at the current node.
       * \param[in] node The current node.
       * \param[in] sol The solution at the current node.
       * \param[in] s_pool Solution pool containing known feasible solutions.
       * \param[out] br_status Status returned by this brancher.
       * \param[out] mods Modification returned by this brancher. NULL if none
       * found.
       */
      virtual Branches findBranches(RelaxationPtr rel, NodePtr node, 
                                    ConstSolutionPtr sol,
                                    SolutionPoolPtr s_pool, 
                                    BrancherStatus & br_status,
                                    ModVector &mods) = 0;

      /// Return the name of this brancher.
      virtual std::string getName() const = 0;

      /**
       * \brief Update pseudo-costs after LP is solved.
       *
       * \param[in] node The node for which relaxation is solved. The pseudo cost of
       * branching candidate used in node->parent is updated.
       * \param[in] sol The solution of the relaxation at this node (not the
       * parent).
       */
      virtual void updateAfterLP(NodePtr node, ConstSolutionPtr sol);

      /// Write statistics to the given out stream.
      virtual void writeStats(std::ostream &) const {};

    protected:
      /// Log manager.
      LoggerPtr logger_;

  };
  typedef boost::shared_ptr<Brancher> BrancherPtr;

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
