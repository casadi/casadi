// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file NodeFullRelaxer.h
 * \brief Declare the NodeFullRelaxer class. It creates relaxation by 
 * fully recreating the relaxation at each node
 * \author Jeff Linderoth, UW-Madison
 */


#ifndef MINOTAURNODEFULLRELAXER_H
#define MINOTAURNODEFULLRELAXER_H

#include "Handler.h"
#include "NodeRelaxer.h"

namespace Minotaur {

  /**
   * The root relaxation is stored as rel_. In each node, we apply all
   * modifications stored in each ancestor of the node. When we are done
   * processing the node, we undo all these changes. 
   *
   * If we dive after processing a node, we do not need to undo all changes
   * and apply them again. We just apply the modifications of the parent.
   */
  class NodeFullRelaxer : public NodeRelaxer {
  public:

    /// 
    NodeFullRelaxer();

    /// Default constructor
    NodeFullRelaxer(EnvPtr env, HandlerVector handlers);

    /// If you know the engine, you can initialize it here
    NodeFullRelaxer(EnvPtr env, EnginePtr e, HandlerVector handlers);

    /// Destroy
    ~NodeFullRelaxer();

    // Implement NodeRelaxer::CreateRootRelaxation()
    RelaxationPtr createRootRelaxation(NodePtr rootNode, bool &prune);

    // Implement NodeRelaxer::CreateNodeRelaxation()
    RelaxationPtr createNodeRelaxation(NodePtr node, bool dived, bool &prune);

    // Implement NodeRelaxer::reset()
    void reset(NodePtr node, bool diving);

    /**
     * Set the engine that is used to solve the relaxations. We need to set
     * it in order to be able to load warm-starts at a node.
     */
    void setEngine(EnginePtr e);

    // get the relaxation pointer, rel_.
    RelaxationPtr getRelaxation();

    /// Set your own relaxation pointer.
    void setRelaxation(RelaxationPtr rel);

  private:
    /// Environment
    EnvPtr env_;

    /**
     * According to NodeRelaxer base class, we must keep the last relaxation
     * created
     */
    RelaxationPtr rel_;

    /// Pointer engine used to solve the relaxation.
    EnginePtr engine_;

    /// Vector of handlers that will make the relaxation.
    HandlerVector handlers_;

    /**
     * We don't update the Node bounds unless improve by at leat this (absolute)
     * amount     
     */
    double updateBoundsTol_;

    /**
    */
    bool isOriginalVariable_(ConstVariablePtr rv, ConstVariablePtr &ov);

    /**
     * \brief Tighten bounds using relaxation loaded to engine.
     *
     * \param[in,out] Modified node
     * \return true if node has been modified
     * 
     * Method assumes that relaxation rel_ is created and loaded to engine
     *
     */
    bool strongBoundsTighten_(NodePtr node);

  };

  typedef boost::shared_ptr <NodeFullRelaxer> NodeFullRelaxerPtr;
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
