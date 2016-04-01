// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file NodeHeap.h
 * \brief Define the class NodeHeap for storing active nodes of the
 * branch-and-bound tree on a heap.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURNODEHEAP_H
#define MINOTAURNODEHEAP_H

#include "ActiveNodeStore.h"

namespace Minotaur {

  /** 
   * When the active nodes of the branch-and-bound tree are not explored in a
   * last-in-first-out (LIFO) order, we store them in a heap. A heap is a
   * binary tree with the properties:
   * -# Every node has a score higher than each of its children.
   * -# If a node is the only child  of its parent, then it is the last node
   *    in the heap.
   * .
   * In order to create a heap, we need to have criteria for comparing nodes.
   * These criteria are determined by the parameter for
   * node-selection-strategy: best bound, best estimate, best estimate+guided
   * dive etc.
   */
   class NodeHeap : public ActiveNodeStore {

      public:
        /// Types of ordering.
        enum Type { Value, Depth };

        /// Constructor.
        NodeHeap(Type type) : type_(type) {}

        /// Destroy.
        virtual ~NodeHeap();

        /**
         * Return true if there are no active nodes in the heap, otherwise
         * return false.
         */
        virtual bool isEmpty() const { return nodes_.empty(); }

        /**
         * Find the minimum lower bound of all the active nodes in the heap.
         * If the heap is ordered by best bound, then the root has the
         * minimum value.
         */
        virtual double getBestLB() const;

        /// Find the maximum depth of all active nodes.
        virtual UInt getDeepestLevel() const;

        /// Remove the best node from the heap.
        virtual void pop();

        /**
         * Write in order the node ID and the criteria used to order the
         * heap, e.g. bound value and depth. 
         */
        virtual void write(std::ostream &out) const;

        /// Add a node to the set of active nodes.
        virtual void push(NodePtr n);

        /// Set the type of ordering of this heap.
        virtual void setType(Type type);

        /// Get access to the best node in this heap.
        virtual NodePtr top() const { return (nodes_.front()); }

        /// Get the number of active nodes in the heap.
        virtual UInt getSize() const { return (nodes_.size()); }

        /// Get iterator to the first node in the heap.
        NodePtrIterator nodesBegin();

        /// Get iterator to the last node in the heap.
        NodePtrIterator nodesEnd();

      private:
        /// Vector of active nodes.
        NodePtrVector nodes_;

        /// The type of criteria used to order the heap.
        Type type_;
   };
   typedef boost::shared_ptr<NodeHeap> NodeHeapPtr;
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
