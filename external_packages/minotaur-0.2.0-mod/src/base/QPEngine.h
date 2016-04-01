// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 


#ifndef MINOTAURQPENGINE_H
#define MINOTAURQPENGINE_H

#include "Engine.h"

namespace Minotaur {

  /**
   * The QPEengine class is an abstract class for interfacing QP solvers (like
   * bqpd). A derived class must implement calls to the QP solver for 
   * the methods described here.
   * 
   * \todo  add more methods for modifying the problem 
   * \todo  add more methods for setting solver options and parameters 
   */

  class QPEngine : public Engine {

    public:
      friend class Problem;

      /// Constructor. May set default parameters/options here.
      QPEngine() {};

      virtual ~QPEngine() {};
    protected:
      /// Destructor must be implemented.
  };
  typedef boost::shared_ptr<QPEngine> QPEnginePtr;
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
