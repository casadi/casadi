// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file NodeStack.cpp
 * \brief Implement the class NodeStack for storing active nodes of the
 * branch-and-bound tree on a stack data structure.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iostream>

#include "MinotaurConfig.h"
#include "Node.h"
#include "NodeStack.h"

using namespace Minotaur;

NodeStack::NodeStack()
{
  // nodes_ is initialized automatically.
}


NodeStack::~NodeStack()
{
  nodes_.clear();
}
  

bool NodeStack::isEmpty() const
{ 
  return nodes_.empty(); 
}


// This function is expensive and must be avoided for large trees.
double NodeStack::getBestLB() const
{
  NodeStackConstIter iter;
  double node_lb;
  double best_lb = INFINITY;

  for (iter = nodes_.begin(); iter != nodes_.end(); ++iter) {
    node_lb = (*iter)->getLb();
    if (node_lb < best_lb) {
      best_lb = node_lb;
    }
  }
  return best_lb;
}


UInt NodeStack::getDeepestLevel() const
{
  NodeStackConstIter iter = nodes_.begin();
  return (*iter)->getDepth();
}


void NodeStack::pop() 
{
  nodes_.pop_front();
}


void NodeStack::push(NodePtr n) 
{
  nodes_.push_front(n);
}


/// Write in order the node ID and the depth of each active node.
void NodeStack::write(std::ostream &out) const 
{
  out << "Nodes in NodeStack:\n";
  for (NodeStackConstIter iter = nodes_.begin(); iter != nodes_.end();
      ++iter) {
    out << "node " << (*iter)->getId() << "\t\tdepth " << 
      (*iter)->getDepth() << std::endl;
  }
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
