// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file QPDRelaxer.h
 * \brief Declare the QPDRelaxer class. 
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURQPDRELAXER_H
#define MINOTAURQPDRELAXER_H

#include "NodeRelaxer.h"

namespace Minotaur {

  class Engine;
  class Logger;
  class Problem;
  typedef boost::shared_ptr<Engine> EnginePtr;
  typedef boost::shared_ptr<const Problem> ConstProblemPtr;


  /**
   * QPDRelaxer creates ``relaxation'' by 
   * creating a QP approximation. It does not really create a relaxation.
   * However, the Relaxation class suffices for its purposes. If we dived on a
   * node, then we don't create another QP, but just change the bounds.
   */
  class QPDRelaxer : public NodeRelaxer {
  public:
    /// Default constructor.
    QPDRelaxer(EnvPtr env, ProblemPtr p, EnginePtr qe, EnginePtr e);

    /// Destroy.
    ~QPDRelaxer();

    // Implement NodeRelaxer::CreateRootRelaxation().
    RelaxationPtr createRootRelaxation(NodePtr rootNode, bool &prune);

    // Implement NodeRelaxer::CreateNodeRelaxation().
    RelaxationPtr createNodeRelaxation(NodePtr node, bool dived, bool &prune); 

    // Implement NodeRelaxer::reset().
    void reset(NodePtr node, bool diving);

    // get the relaxation pointer, qp_.
    RelaxationPtr getRelaxation();

  private:
    /// Environment
    EnvPtr env_;

    /// Original Problem.
    ProblemPtr p_;

    /**
     * We only keep one QP formulation. We will clear all constraints (but
     * not variables) and rebuild them if needed.
     */
    RelaxationPtr qp_;

    /// Engine used to solve QP.
    EnginePtr qpe_;

    /// Engine used to solve NLP.
    EnginePtr e_;

    /// Logger.
    Logger *logger_;

  };

  typedef boost::shared_ptr <QPDRelaxer> QPDRelaxerPtr;
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
