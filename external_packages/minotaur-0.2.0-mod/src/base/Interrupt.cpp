//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

#include <csignal>
#include <cstdlib>

#include "MinotaurConfig.h"
#include "Interrupt.h"

using namespace Minotaur;

volatile UInt Interrupt::interrupts_ = 0;
UInt Interrupt::limit_ = 0;

void Interrupt::Handler(int) {
  signal(SIGINT, SIG_IGN);
  ++interrupts_;

  if (interrupts_ >= limit_) {
    //(logger_)->ErrStream() << "Interrupt limit reached." << std::endl;
    exit(-1);
  }
  signal(SIGINT, Interrupt::Handler);
  return;
}

void Interrupt::Start() {
  originalHandler_ = signal(SIGINT, Interrupt::Handler);
  return;
}

void Interrupt::Stop() {
  signal(SIGINT, originalHandler_);
  return;
}

void Interrupt::Check() const throw(InterruptException) {
  if (interrupts_ > 0) {
    throw InterruptException();
  }
  return;
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
