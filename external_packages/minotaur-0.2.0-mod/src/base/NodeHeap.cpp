// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file NodeHeap.cpp
 * \brief Define class NodeHeap for storing active nodes.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cassert>
#include <cmath>

#include "MinotaurConfig.h"
#include "Node.h"
#include "NodeHeap.h"

namespace Minotaur {
   bool valueGreaterThan(ConstNodePtr n1, ConstNodePtr n2) 
   {
      if (n1->getLb() > n2->getLb() + 1e-6) { 
        return true;
      } else if (n1->getLb() < n2->getLb() - 1e-6) {
        return false;
      } 

      if (n1->getTbScore() > n2->getTbScore() + 1e-6) {
        return true;
      } else if (n1->getTbScore() < n2->getTbScore() - 1e-6) {
        return false;
      } 

      if (n1->getDepth() < n2->getDepth()) {
        return false;
      } else if (n1->getDepth() > n2->getDepth()) {
        return true;
      }

      return (n1->getId() < n2->getId());
   }

   bool depthGreaterThan(ConstNodePtr n1, ConstNodePtr n2)
   {
      return (n1->getDepth() > n2->getDepth());
   }
}

using namespace Minotaur;


NodeHeap::~NodeHeap()
{
  nodes_.clear();
}


void NodeHeap::pop()
{
  // std::cout << "popping out node " << nodes_.front()->getId() << std::endl;
  switch(type_) {
  case (Value):
    pop_heap(nodes_.begin(), nodes_.end(), valueGreaterThan);
    break;
  case (Depth):
    pop_heap(nodes_.begin(), nodes_.end(), depthGreaterThan);
    break;
  default:
    assert(0);
  }
  nodes_.pop_back();
}


double NodeHeap::getBestLB() const
{
   double retval = INFINITY;
   if (nodes_.size() > 0) {
      if (type_ == Value) {
         retval = nodes_.front()->getLb();
      } else {
        std::vector<NodePtr>::const_iterator it = min_element(nodes_.begin(), 
             nodes_.end(), valueGreaterThan);
         retval = (*it)->getLb();
      }
   }
   return retval;
}


UInt NodeHeap::getDeepestLevel() const
{
   UInt retval = 0;
   if (nodes_.size() > 0) {
      if (type_ == Depth) {
         retval = nodes_.front()->getDepth();
      }
      else {
        std::vector<NodePtr>::const_iterator it = max_element(nodes_.begin(), 
             nodes_.end(), depthGreaterThan);
         retval = (*it)->getDepth();
      }
   }
   return retval;
}


void NodeHeap::write(std::ostream &) const
{
   //for(std::vector<NodePtr>::const_iterator it = nodes_.begin();
   //    it != nodes_.end(); it++) {
   //   ConstNodePtr n = *it;
   //   //o << (*n) << endl;
   //}
}


void NodeHeap::push(NodePtr n)
{
  // std::cout << "inserting node " << n->getId() << std::endl;
  nodes_.push_back(n);
  switch(type_) {
  case (Value):
    push_heap(nodes_.begin(), nodes_.end(), valueGreaterThan);
    break;
  case (Depth):
    push_heap(nodes_.begin(), nodes_.end(), depthGreaterThan);
    break;
  default:
    assert(0);
  }
}


void NodeHeap::setType(Type type)
{
  if (type == type_) return;
  //XXX Could probably do this fancier

  switch(type_) {
  case (Value):
    make_heap(nodes_.begin(), nodes_.end(), valueGreaterThan);
    break;
  case (Depth):
    make_heap(nodes_.begin(), nodes_.end(), depthGreaterThan);
    break;
  default:
    assert(0);
  }
}


NodePtrIterator NodeHeap::nodesBegin() 
{
  return nodes_.begin();
}


NodePtrIterator NodeHeap::nodesEnd() 
{
  return nodes_.end();
}

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
