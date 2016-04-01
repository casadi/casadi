// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 


#include "MinotaurConfig.h"
#include "Engine.h"
#include "Logger.h"

using namespace Minotaur;

Engine::Engine()
  : logger_(LoggerPtr()) // NULL
{

}


Engine::~Engine()
{
  logger_.reset();
}


std::string Engine::getStatusString()
{
  switch (status_) {
   case(ProvenOptimal):
     return "ProvenOptimal";
   case(ProvenLocalOptimal):
     return "ProvenLocalOptimal";
   case(ProvenInfeasible):
     return "ProvenInfeasible";
   case(ProvenLocalInfeasible):
     return "ProvenLocalInfeasible ";
   case(ProvenUnbounded):
     return "ProvenUnbounded";
   case(ProvenObjectiveCutOff):
     return "ProvenObjectiveCutOff";
   case(EngineIterationLimit):
     return "EngineIterationLimit";
   case(ProvenFailedCQFeas):
    return "ProvenFailedCQFeas";
   case(ProvenFailedCQInfeas):
    return "ProvenFailedCQInfeas";
   case(FailedFeas):
    return "FailedFeas";
   case(FailedInfeas):
    return "FailedInfeas";
   case(EngineError):
     return "EngineError";
   case(EngineUnknownStatus):
     return "EngineUnknownStatus";
   default:
     return "Engine status unknown type";
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
