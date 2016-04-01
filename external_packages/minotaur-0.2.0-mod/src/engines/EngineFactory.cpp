//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file EngineFactory.cpp
 * \brief EngineFactory class finds a suitable engine.
 * \author Sven Leyffer, Argonne National Laboratory
 */

//#define SPEW 1
#include "MinotaurConfig.h"
#include "EngineFactory.h"
#include "Environment.h"
#include "LPEngine.h"
#include "NLPEngine.h"
#include "Option.h"
#include "QPEngine.h"

#ifdef USE_IPOPT
#include "IpoptEngine.h"
#endif

#ifdef USE_OSILP
#include "OsiLPEngine.h"
#endif

#ifdef USE_FILTERSQP
#include "FilterSQPEngine.h"
#endif

#ifdef USE_BQPD
#include "BqpdEngine.h"
#endif

#ifdef USE_QPOASES
#include "qpOASESEngine.h"
#endif

using namespace Minotaur;

EngineFactory::EngineFactory()
  : env_(EnvPtr())
{

}


EngineFactory::EngineFactory(EnvPtr env)
  : env_(env)
{

}


EngineFactory::~EngineFactory()
{
  env_.reset();
}


LPEnginePtr EngineFactory::getLPEngine()
{
#ifdef USE_OSILP
  if (env_->getOptions()->findString("lp_engine")->getValue()!="None") {
    return ((OsiLPEnginePtr) new OsiLPEngine(env_));
  }
#endif
  return (LPEnginePtr());
}


QPEnginePtr EngineFactory::getQPEngine()
{
#ifdef USE_BQPD
  if (env_->getOptions()->findString("qp_engine")->getValue()=="bqpd") {
    return ((BqpdEnginePtr) new BqpdEngine(env_));
  }
#endif
#ifdef USE_QPOASES
  if (env_->getOptions()->findString("qp_engine")->getValue()=="qpOASES") {
    return ((qpOASESEnginePtr) new qpOASESEngine(env_));
  }
#endif
  return (QPEnginePtr());
}


NLPEnginePtr EngineFactory::getNLPEngine()
{
  if (env_->getOptions()->findString("nlp_engine")->getValue()=="Filter-SQP") {
#ifdef USE_FILTERSQP
    return ((FilterSQPEnginePtr) new FilterSQPEngine(env_));
#endif
#ifdef USE_IPOPT
    return ((IpoptEnginePtr) new IpoptEngine(env_));
#endif
  }

  if (env_->getOptions()->findString("nlp_engine")->getValue()=="IPOPT") {
#ifdef USE_IPOPT
    return ((IpoptEnginePtr) new IpoptEngine(env_));
#endif
#ifdef USE_FILTERSQP
    return ((FilterSQPEnginePtr) new FilterSQPEngine(env_));
#endif
  }
  return (NLPEnginePtr()); // NULL
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
