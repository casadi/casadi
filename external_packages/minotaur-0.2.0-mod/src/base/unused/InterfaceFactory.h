// Base method for generating information from the user interface
// Concrete instances are generated from GAMS/AMPL/CUTEr

#ifndef INTERFACEFACTORY_H
#define INTERFACEFACTORY_H

namespace Minotaur {
  class Environment;	// Logger, interrupt handlers, timers, interface data
  class Instance;	// The optimization problem

  class InterfaceFactory {
    InterfaceFactory();

    virtual void create(Environment *e, Instance *i) = 0;
  };
}

#endif

