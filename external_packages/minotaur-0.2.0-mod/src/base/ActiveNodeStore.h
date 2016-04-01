// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file ActiveNodeStore.h
 * \brief Declare the base class ActiveNodeStore for storing active nodes of
 * the branch-and-bound tree.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURACTIVENODESTORE_H
#define MINOTAURACTIVENODESTORE_H

#include "Types.h"

namespace Minotaur {

  class Node;
  typedef boost::shared_ptr<Node> NodePtr;

  /** 
   * \brief Save and retrieve active nodes.
   *
   * Nodes that have not been processed yet by the branch-and-bound algorithm
   * are called active nodes. This base class is used to store them in a
   * suitable data structure (stack or heap).
   */
   class ActiveNodeStore {

     public:
       /// Default Constructor.
       ActiveNodeStore () {}

       /// Destroy.
       virtual ~ActiveNodeStore() {}

       /// Find the minimum lower bound of all the active nodes.
       virtual double getBestLB() const = 0;

       /// Find the maximum depth of all active nodes.
       virtual UInt getDeepestLevel() const = 0;

       /// Get the number of active nodes.
       virtual UInt getSize() const = 0;

       /**
        * \brief Check if there are any active nodes left.
        *
        * \return true if there are no active nodes in the heap, otherwise
        * return false.
        */
       virtual bool isEmpty() const = 0;

       /// Remove the best node from the store.
       virtual void pop() = 0;

       /**
        * \brief Add a node to the set of active nodes.
        *
        * \param[in] n Pointer to node to be added.
        */
       virtual void push(NodePtr n) = 0;

       /** 
        * \brief Access to the best candidate for evaluating next.
        *
        * Does not remove the best candidate from the storage. It just returns
        * a pointer to it. Use pop() to remove the candidate.
        * \return pointer to the best candidate node.
        */
       virtual NodePtr top() const = 0;

       /**
        * \brief Display the active nodes.
        *
        * Write in order the node IDs and the criteria used to order the
        * nodes, e.g. bound value and depth. 
        * \param[in] out output stream to write to.
        */
       virtual void write(std::ostream &out) const = 0;
   };
   typedef boost::shared_ptr<ActiveNodeStore> ActiveNodeStorePtr;
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
