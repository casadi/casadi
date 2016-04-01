// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

// /**
// \file WarmStart.h
// \brief Define the base class WarmStart for storing the warm starting
// information for different types of engines.
// \author Ashutosh Mahajan, Argonne National Laboratory
// */


#ifndef MINOTAURWARMSTART_H
#define MINOTAURWARMSTART_H

#include "Types.h"

namespace Minotaur {

  class WarmStart;
  typedef boost::shared_ptr<WarmStart> WarmStartPtr;
  typedef boost::shared_ptr<const WarmStart> ConstWarmStartPtr;

  // /** 
  // Warm starting information enables an engine to quickly resolve a problem
  // after that has been modified. 
  //
  // We save warm start information on all those
  // nodes of the branch-and-bound tree whose at least one child is saved for
  // later processing. If all the children of a node are pruned or if the node
  // has only one child and we decide to process it next, then we don't need
  // to save warm start information for that node.
  //
  // For now we save complete warm-start information on each active node. A more
  // memory efficient method is to save warm-start information for each node
  // by just storing the differences from the parent.
  // However, the benefits of saving the complete information are:
  //  -# Ease of coding.
  //  -# We can delete warm-start information of a node when it is processed.
  //  Thus, we have at most `A' nodes that have warm-start information saved
  //  on them, where `A' is the total number of active nodes in the tree.
  // */
  class WarmStart {
    public:
      /// Default constructor
      WarmStart() {}

      /// Destroy
      virtual ~WarmStart() {}
      
      /// Return true if warm start information is initialized, false
      /// otherwise.
      virtual bool hasInfo() = 0;

      /// Write to an output stream
      virtual void write(std::ostream &out) const = 0;

      /// \todo
      //virtual void setEngine(EnginePtr engine) = 0;
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
