// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 -- 2010 The MINOTAUR Team.
// 

// /**
// \file SVector.h
// \brief Define sparse vectors
// \author Ashutosh Mahajan, Argonne National Laboratory
//
// This file contains definitions for methods of sparse vectors.
// */

#ifndef MINOTAURSVECTOR_H
#define MINOTAURSVECTOR_H

#include "Types.h"
#include "Vector.h"

namespace Minotaur {
  namespace LinearAlgebra {
    template <class T> class SVector : public Vector<T>  {
    private:
      UInt  max_;	// Maximum number of nonzeros
      UInt  nnz_;	// Current number of nonzeros
      UInt *idx_;	// Index of entries
      T    *dat_;	// Coefficient values

    public:
      SVector(UInt len, UInt nnz);
      SVector(UInt len, UInt nnz, UInt max);
      virtual ~SVector();
    };

    typedef SVector<UInt> IndexSet;
  };
};

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
