// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file NodeIncRelaxer.cpp
 * \brief Define the class NodeIncRelaxer.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <iostream>
#include <stack>

#include "MinotaurConfig.h"
#include "Engine.h"
#include "Environment.h"
#include "Node.h"
#include "NodeIncRelaxer.h"
#include "Option.h"
#include "Relaxation.h"


using namespace Minotaur;


NodeIncRelaxer::NodeIncRelaxer (EnvPtr env, HandlerVector handlers) 
  : engine_(EnginePtr()),  // NULL
    env_(env),
    handlers_(handlers),
    modProb_(true),
    rel_(RelaxationPtr()) // NULL
{
}


NodeIncRelaxer::~NodeIncRelaxer ()
{
  rel_.reset();
  handlers_.clear();
  env_.reset();
  engine_.reset();
}


RelaxationPtr NodeIncRelaxer::createRootRelaxation(NodePtr, bool &prune)
{
  prune = false;
  rel_ = (RelaxationPtr) new Relaxation();
  for (HandlerIterator h = handlers_.begin(); h != handlers_.end() && !prune; 
      ++h) {
    (*h)->relaxInitInc(rel_, &prune);
    if (prune) {
      break;
    }
  }
  if (!prune) {
    rel_->calculateSize();
  }

  if (engine_) {
    engine_->load(rel_);
  }

  return rel_;
}


bool NodeIncRelaxer::getModFlag()
{
  return modProb_;
}


// Set and load the engine that is used to solve the relaxation.
void NodeIncRelaxer::setEngine(EnginePtr e)
{
  engine_ = e;
  if (rel_) {
    engine_->load(rel_);
  }
}


void NodeIncRelaxer::setModFlag(bool mod_prob)
{
  modProb_ = mod_prob;
}


RelaxationPtr NodeIncRelaxer::createNodeRelaxation(NodePtr node, bool dived, 
                                                   bool &prune)
{
  NodePtr t_node; // temporary
  WarmStartPtr ws;
  prune = false;

  if (!dived) {
    // traceback to root and put in all modifications that need to go into the
    // relaxation and the engine.
    std::stack<NodePtr> predecessors;
    t_node = node->getParent();

    while (t_node) {
      predecessors.push(t_node);
      t_node = t_node->getParent();
    }

    // starting from the top, put in modifications made at each node to the
    // engine
    if (modProb_) {
      while (!predecessors.empty()) {
        t_node = predecessors.top();
        t_node->applyMods(rel_, p_);
        predecessors.pop();
      }
    } else {
      while (!predecessors.empty()) {
        t_node = predecessors.top();
        t_node->applyRMods(rel_);
        predecessors.pop();
      }
    }
  } 

  // put in the modifications that were used to create this node from
  // its parent.
  if (modProb_) {
    node->applyMods(rel_, p_);
  } else {
    node->applyRMods(rel_);
  }

  for (HandlerIterator h = handlers_.begin(); h != handlers_.end() && !prune; 
      ++h) {
    (*h)->relaxNodeInc(node, rel_, &prune);
    if (prune) {
      break;
    }
  }

  // give the warm start info to the engine.
  if (!prune) {
    ws = node->getWarmStart();
    if (ws) {
      engine_->loadFromWarmStart(ws);
    }
  } else if (ws) {
    node->removeWarmStart();
  }
  return rel_;
}


void NodeIncRelaxer::reset(NodePtr node, bool diving)
{
  NodePtr t_node = node;

  if (!diving) {
    if (modProb_) {
      while (t_node) {
        t_node->undoMods(rel_, p_);
        t_node = t_node->getParent();
      }
    } else {
      while (t_node) {
        t_node->undoRMods(rel_);
        t_node = t_node->getParent();
      }
    }
  }
}


RelaxationPtr NodeIncRelaxer::getRelaxation()
{
  return rel_;
}


void NodeIncRelaxer::setRelaxation(RelaxationPtr rel)
{
  rel_ = rel;
}


void NodeIncRelaxer::setProblem(ProblemPtr p)
{
  p_ = p;
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
