//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

#include <iostream>

#include "MinotaurConfig.h"
#include "Logger.h"

using namespace Minotaur;

Logger::Logger(LogLevel max_level) 
  : maxLevel_(max_level),  nb_(), nout_(&nb_)
{

}


Logger::~Logger() 
{
  // Destructor 
}


std::ostream& Logger::msgStream(LogLevel level) const 
{
  if (level <= maxLevel_) { 
    return std::cout; 
  }
  return nout_;
}


std::ostream& Logger::errStream() const 
{ 
  if (maxLevel_ > LogNone) { 
    return std::cerr; 
  } else {
    return nout_;
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
