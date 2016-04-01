// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file LPEngine.h
 * \author Ashutosh Mahajan, Argonne National Laboratory.
 * \brief Declare the class LPEngine for solving LPs and getting solution.
 */


#ifndef MINOTAURLPENGINE_H
#define MINOTAURLPENGINE_H

#include "Engine.h"

namespace Minotaur {

  /**
   * The LPEengine class is an abstract class for interfacing LP solvers (like
   * OsiLPEngine). A derived class must implement calls to the LP solver for 
   * the methods described here.
   * 
   * \todo  add more methods for accessing simplex tableaux.
   * \todo  add more methods for accessing dual rays etc.
   */

  class LPEngine : public Engine {

    public:
      friend class Problem;

      /// Constructor. May set default parameters/options here.
      LPEngine() {};

      /// Destructor must be implemented if memory needs to be freed
      virtual ~LPEngine() {};
  };
  typedef boost::shared_ptr<LPEngine> LPEnginePtr;
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
