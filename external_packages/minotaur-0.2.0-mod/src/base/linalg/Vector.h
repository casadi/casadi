
//     Minotaur -- It's only 1/2 bull
//
//     (C)opyright 2009 -- The Minotaur Team

#ifndef MINOTAURVECTOR_H
#define MINOTAURVECTOR_H

#include "Types.h"

namespace Minotaur {
  namespace LinearAlgebra {
    template <class T> class Vector {
    protected:
      UInt len_;

    public:
      Vector(UInt len) : len_(len) { }
      virtual ~Vector() { }

      UInt getNumNz();
    };

    typedef boost::shared_ptr< Vector<UInt> > UIntVectorPtr;
    typedef boost::shared_ptr< Vector<const UInt> > ConstUIntVectorPtr;

    typedef boost::shared_ptr< Vector<double> > DoubleVectorPtr;
    typedef boost::shared_ptr< Vector<const double> > ConstDoubleVectorPtr;
  }
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
