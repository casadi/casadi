//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//


#ifndef MINOTAUREXCEPTION_H
#define MINOTAUREXCEPTION_H

namespace Minotaur {
  class Exception {
  public:
    Exception() { }
    Exception(const Exception&) { }		// Copy constructor
    virtual ~Exception() { }			// Destructor

  private:
    Exception& operator=(const Exception&);	// Copy assignment
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
