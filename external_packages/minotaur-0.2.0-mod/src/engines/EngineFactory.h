//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//


/**
 * \file EngineFactory.h
 * \brief Declare EngineFactory class for creating different engines.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURENGINEFACTORY_H
#define MINOTAURENGINEFACTORY_H

#include "Types.h"

namespace Minotaur {

  class Engine;
  class LPEngine;
  class NLPEngine;
  class QPEngine;
  typedef boost::shared_ptr<Engine> EnginePtr;
  typedef boost::shared_ptr<LPEngine> LPEnginePtr;
  typedef boost::shared_ptr<NLPEngine> NLPEnginePtr;
  typedef boost::shared_ptr<QPEngine> QPEnginePtr;

  class EngineFactory {
    public:
      /// Default constructor
      EngineFactory();

      /// Constructor with minotaur environment.
      EngineFactory(EnvPtr env);

      /// Destroy.
      ~EngineFactory();

      /// Get an engine. Returns NULL if no engine is available.
      EnginePtr getEngine();

      /// Get an LP Engine
      LPEnginePtr getLPEngine();

      /// Get a QP Engine. Returns NULL if none available.
      QPEnginePtr getQPEngine();

      /// Get an NLP Engine. Returns NULL if none available.
      NLPEnginePtr getNLPEngine();

    private:
      /// Minotaur Environment.
      EnvPtr env_;
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
